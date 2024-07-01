#!/bin/bash

echo -e "PREPROCESSING STARTED\n"

# Arguments and Variables

baseDir=$(dirname "$0") # folder containing all scripts for preprocessing

inputDir=$1
r1_file=($(ls --ignore=*trimmed* ${inputDir} | grep "\\.R1\\.fastq"))
r1_file=${r1_file[0]}
r2_file=($(ls --ignore=*trimmed* ${inputDir} | grep "\\.R2\\.fastq"))
r2_file=${r2_file[0]}
sample=${r1_file%.R1*}
outputDir=$2 # a $sample folder will be created inside
strand=$3 # R1, R2, or merge
primersInfo=$4 # 6 columns: primer name, seq forward, seq reverse, seq reverse complement, ref size, with CG boolean
umisSize=$5 # 16
threads=$6 # 12

if [ "${#r1_file[@]}" != "1" -o "${#r2_file[@]}" != "1" ]; then
    exit "Invalid number of R1/R2 files."
fi

echo -e "VARIABLES EXTRACTED FROM ARGUMENTS ARE:"
echo -e "- input directory: ${inputDir}"
echo -e "- r1 file: ${r1_file}"
echo -e "- r2 file: ${r2_file}"
echo -e "- sample name: ${sample}"
echo -e "- output directory: ${outputDir}"
echo -e "- strand used: ${strand}"
echo -e "- primers info file path: ${primersInfo}"
echo -e "- UMIs size: ${umisSize}"
echo -e "- number of threads: ${threads}\n"

# Filter reads by quality using fastp and fasta conversion (with merge if needed)
echo "FILTERING READS"

mkdir -p ${outputDir}/${sample}
echo "sample,reads,quality filtered reads,merged reads,reads post cutting forward,inserts,length filtered inserts,deduplicated length filtered inserts,UMIs,length filtered UMIs" > ${outputDir}/${sample}/${sample}.summary.csv

if [ ${strand} = merge ]; then

    ## Filter
    fastp -i ${inputDir}/${r1_file} -I ${inputDir}/${r2_file} -o ${outputDir}/${sample}/${sample}.R1.fastq -O ${outputDir}/${sample}/${sample}.R2.fastq -h ${outputDir}/${sample}/${sample}.fastp.html -j ${outputDir}/${sample}/${sample}.fastp.json
    
    ## add info to header
    illuminaTag=$(head -n 1 ${outputDir}/${sample}/${sample}.R1.fastq | grep -P -o "@(.+?):" | tr -d '@:')
    sed -i "s/${illuminaTag}/${sample}/g" ${outputDir}/${sample}/${sample}.R1.fastq
    sed -i "s/${illuminaTag}/${sample}/g" ${outputDir}/${sample}/${sample}.R2.fastq
    
    ## merge R1/R2
    cmd="mothur '#make.contigs(ffastq=${outputDir}/${sample}/${sample}.R1.fastq, rfastq=${outputDir}/${sample}/${sample}.R2.fastq, trimoverlap=F, processors=${threads})'"
    eval ${cmd}
    sed -i "s/\tee=.*$//g" ${outputDir}/${sample}/${sample}.R1.trim.contigs.fasta # Strip Mothur's error rate output
    mv ${outputDir}/${sample}/${sample}.R1.trim.contigs.fasta ${outputDir}/${sample}/${sample}.fasta
    cat ${outputDir}/${sample}/${sample}.fasta | grep -v ">" | awk '{print length}' > ${outputDir}/${sample}/reads_length.${sample}.csv # get reads length after merge
    
    ## infos for summary
    startReadsCount=$(zcat ${inputDir}/${r1_file} | grep "^+$" | wc -l)
	qualityFilteredReadsCount=$(grep "^+$" ${outputDir}/${sample}/${sample}.R1.fastq | wc -l)
	mergedReadsCount=$(grep ">" ${outputDir}/${sample}/${sample}.fasta | wc -l)

elif [ ${strand} = R1 ]; then

    ## Filter
    fastp -i ${inputDir}/${r1_file} -o ${outputDir}/${sample}/${sample}.fastq -h ${outputDir}/${sample}/${sample}.fastp.html -j ${outputDir}/${sample}/${sample}.fastp.json
    
    ## add info to header
    illuminaTag=$(head -n 1 ${outputDir}/${sample}/${sample}.fastq | grep -P -o "@(.+?):" | tr -d '@:')
    sed -i "s/${illuminaTag}/${sample}/g" ${outputDir}/${sample}/${sample}.fastq
    
    ## conversion to fasta
    seqkit fq2fa ${outputDir}/${sample}/${sample}.fastq -o ${outputDir}/${sample}/${sample}.fasta \
    && rm ${outputDir}/${sample}/${sample}.fastq
    
    ## infos for summary
    startReadsCount=$(zcat ${inputDir}/${r1_file} | grep "^+$" | wc -l)
	qualityFilteredReadsCount=$(grep ">" ${outputDir}/${sample}/${sample}.fasta | wc -l)
	mergedReadsCount="NA"

elif [ ${strand} = R2 ]; then

    ## Reverse complement
    vsearch --fastx_revcomp ${inputDir}/${r2_file} -fastqout ${outputDir}/${sample}/${sample}.reversed.fastq

    ## Filter
    fastp -i ${outputDir}/${sample}/${sample}.reversed.fastq -o ${outputDir}/${sample}/${sample}.fastq -h ${outputDir}/${sample}/${sample}.fastp.html -j ${outputDir}/${sample}/${sample}.fastp.json \
    && rm ${outputDir}/${sample}/${sample}.reversed.fastq

    ## add info to header
    illuminaTag=$(head -n 1 ${outputDir}/${sample}/${sample}.fastq | grep -P -o "@(.+?):" | tr -d '@:')
    sed -i "s/${illuminaTag}/${sample}/g" ${outputDir}/${sample}/${sample}.fastq

    ## conversion to fasta
    #paste - - - - < ${outputDir}/${sample}/${sample}.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${outputDir}/${sample}/${sample}.fasta \
    seqkit fq2fa ${outputDir}/${sample}/${sample}.fastq -o ${outputDir}/${sample}/${sample}.fasta \
    && rm ${outputDir}/${sample}/${sample}.fastq

    ## infos for summary
    startReadsCount=$(zcat ${inputDir}/${r2_file} | grep "^+$" | wc -l)
	qualityFilteredReadsCount=$(grep ">" ${outputDir}/${sample}/${sample}.fasta | wc -l)
	mergedReadsCount="NA"

fi

# For each primer
echo "PROCESSING EACH PRIMER"

awk -F',' '(NR!=1) { split($0, line); print line[1],line[2],line[3],line[4],line[5],line[6] }' ${primersInfo} | 
while read -r primer forward reverse reverse_revcomp refsize withCG; do 

    # Trimming

    ## trim forward...
    atropos trim -g ${forward} -se ${outputDir}/${sample}/${sample}.fasta -o ${outputDir}/${sample}/${sample}.${primer}.forward_trim.fasta --threads ${threads} -y ":${primer}" --discard-untrimmed --no-default-adapters --overlap 7 --info-file=${outputDir}/${sample}/${sample}.${primer}.forward_trim.info.csv \
    && awk -F'\t' '$2!=-1' ${outputDir}/${sample}/${sample}.${primer}.forward_trim.info.csv > ${outputDir}/${sample}/${sample}.${primer}.forward_trim.info.minus_one_filtered.csv \
    && rm -f ${outputDir}/${sample}/${sample}.${primer}.forward_trim.info.csv
    ## + isolate forward primers in a separate fasta file (not needed - already in the previous csv)
    #python ${baseDir}/info_file_to_fasta.py ${outputDir}/${sample}/${sample}.${primer}.forward_trim.info.minus_one_filtered.csv --primer ${primer} --end forward_primer > ${outputDir}/${sample}/forward_primers.${sample}.${primer}.fasta
    ## ... then trim reverse
    atropos trim -a ${reverse_revcomp} -se ${outputDir}/${sample}/${sample}.${primer}.forward_trim.fasta -o ${outputDir}/${sample}/${sample}.${primer}.trimmed.fasta --threads ${threads} --discard-untrimmed --no-default-adapters --overlap 7 --info-file=${outputDir}/${sample}/${sample}.${primer}.trimmed.info.csv \
    && awk -F'\t' '$2!=-1' ${outputDir}/${sample}/${sample}.${primer}.trimmed.info.csv > ${outputDir}/${sample}/${sample}.${primer}.trimmed.info.minus_one_filtered.csv \
    && rm -f ${outputDir}/${sample}/${sample}.${primer}.trimmed.info.csv
    ## + isolate reverse primers in a separate fasta file (not needed - already in the previous csv)
    #python ${baseDir}/info_file_to_fasta.py ${outputDir}/${sample}/${sample}.${primer}.trimmed.info.minus_one_filtered.csv --primer ${primer} --end reverse_primer > ${outputDir}/${sample}/reverse_primers.${sample}.${primer}.fasta
    ## plot reads length after complete trim
    grep -v ">" ${outputDir}/${sample}/${sample}.${primer}.trimmed.fasta | awk '{print length}' > ${outputDir}/${sample}/reads_length.${sample}.${primer}.csv
    Rscript ${baseDir}/plot_read_or_uid_length.R -t "Reads length ${sample}, primer ${primer} - exp. size = ${refsize}" -o ${outputDir}/${sample}/reads_length.${sample}.${primer}.png ${outputDir}/${sample}/reads_length.${sample}.${primer}.csv

    # UMIs

    ## Isolate UMIs
    python ${baseDir}/info_file_to_fasta.py ${outputDir}/${sample}/${sample}.${primer}.trimmed.info.minus_one_filtered.csv --primer ${primer} --end reverse_umi > ${outputDir}/${sample}/UMIs.${sample}.${primer}.fasta
    ## Filter UMIs by size
    vsearch --fastx_filter ${outputDir}/${sample}/UMIs.${sample}.${primer}.fasta --notrunclabels --fastaout ${outputDir}/${sample}/UMIs.length_filtered.${sample}.${primer}.fasta --fastq_minlen ${umisSize} --fastq_maxlen ${umisSize}
    touch ${outputDir}/${sample}/UMIs.length_filtered.${sample}.${primer}.fasta # Make sure a file exists for all primers, even empty
    ## Extract length filtered UMIs to a csv file 
    cat ${outputDir}/${sample}/UMIs.length_filtered.${sample}.${primer}.fasta | grep -v ">" > ${outputDir}/${sample}/${sample}.${primer}.UMIs.length_filtered.csv
    ## Extract UMIs frequencies
    Rscript ${baseDir}/adapter_stats.R -t "UMIs frequencies count ${sample} - ${primer}" -o ${outputDir}/${sample}/UMIs_frequencies_count_plots.${sample}.${primer}.png ${outputDir}/${sample}/${sample}.${primer}.UMIs.length_filtered.csv
    touch ${outputDir}/${sample}/${sample}.${primer}.UMIs.length_filtered.frequency.csv # Make sure a file exists for all primers, even empty

    # Inserts

    ## Filter inserts by size
    vsearch --fastx_filter ${outputDir}/${sample}/${sample}.${primer}.trimmed.fasta --notrunclabels --fastaout ${outputDir}/${sample}/inserts.length_filtered.${sample}.${primer}.fasta --fastq_minlen $(echo "${refsize} - 5" | bc) --fastq_maxlen $(echo "${refsize} + 5" | bc)
    touch ${outputDir}/${sample}/inserts.length_filtered.${sample}.${primer}.fasta # Make sure a file exists for all primers, even empty

    # Deduplication based on UMIs

    ## Concatenate length filtered inserts with length filtered UMIs
    python ${baseDir}/concatenate_UMIs_targets.py ${outputDir}/${sample}/UMIs.length_filtered.${sample}.${primer}.fasta ${outputDir}/${sample}/inserts.length_filtered.${sample}.${primer}.fasta > ${outputDir}/${sample}/inserts_UMIs_concat.length_filtered.${sample}.${primer}.fasta
    ## Deduplicate
    vsearch --notrunclabels --derep_fulllength ${outputDir}/${sample}/inserts_UMIs_concat.length_filtered.${sample}.${primer}.fasta --sizein --output ${outputDir}/${sample}/inserts_UMIs_concat_dereplicated.length_filtered.${sample}.${primer}.fasta \
    && rm -f ${outputDir}/${sample}/inserts_UMIs_concat.length_filtered.${sample}.${primer}.fasta
    ## Trim UMIs after deduplication (UMIs are appended to reads ids)
    python ${baseDir}/trim_UMIs.py ${outputDir}/${sample}/inserts_UMIs_concat_dereplicated.length_filtered.${sample}.${primer}.fasta > ${outputDir}/${sample}/inserts_dereplicated.length_filtered.${sample}.${primer}.fasta \
    && rm -f ${outputDir}/${sample}/inserts_UMIs_concat_dereplicated.length_filtered.${sample}.${primer}.fasta
    touch ${outputDir}/${sample}/inserts_dereplicated.length_filtered.${sample}.${primer}.fasta # Make sure a file exists for all primers, even empty

done
## concatenate each primer specific file into one
cat ${outputDir}/${sample}/${sample}.*.UMIs.length_filtered.csv >> ${outputDir}/${sample}/${sample}.merged.UMIs.length_filtered.csv
## Extract UMIs frequencies for merged file
Rscript ${baseDir}/adapter_stats.R -t "UMIs frequencies count ${sample} - all primers" -o ${outputDir}/${sample}/UMIs_frequencies_count_plots.${sample}.merged.png ${outputDir}/${sample}/${sample}.merged.UMIs.length_filtered.csv

# infos for summary
echo "GETTING REMAINING INFOS"

forwardCutReadsCount=$(grep ">" ${outputDir}/${sample}/${sample}.*.forward_trim.fasta | wc -l) \
&& rm -f ${outputDir}/${sample}/${sample}.*.forward_trim.fasta # sequences are already in the csv file
insertsReadsCount=$(grep ">" ${outputDir}/${sample}/${sample}.*.trimmed.fasta | wc -l) \
lengthFilteredInsertsReadsCount=$(grep ">" ${outputDir}/${sample}/inserts.length_filtered.${sample}.*.fasta | wc -l)
deduplicatedLengthFilteredInsertsReadsCount=$(grep ">" ${outputDir}/${sample}/inserts_dereplicated.length_filtered.${sample}.*.fasta | wc -l)
umisCount=$(grep ">" ${outputDir}/${sample}/UMIs.${sample}.*.fasta | wc -l) \
&& rm -f ${outputDir}/${sample}/UMIs.${sample}.*.fasta ## keeping only the length_filtered one
lengthFilteredUmisCount=$(grep ">" ${outputDir}/${sample}/UMIs.length_filtered.${sample}.*.fasta | wc -l)
echo "${sample},${startReadsCount// /},${qualityFilteredReadsCount// /},${mergedReadsCount// /},${forwardCutReadsCount// /},${insertsReadsCount// /},${lengthFilteredInsertsReadsCount// /},${deduplicatedLengthFilteredInsertsReadsCount// /},${umisCount// /},${lengthFilteredUmisCount// /}" >> ${outputDir}/${sample}/${sample}.summary.csv

# Count reads for each primer set + UMI stats

primers=$(awk -F',' '(NR!=1) { print $1 }' ${primersInfo} )

echo "sample,$(echo ${primers} | sed 's/ /,/g'),merged" > ${outputDir}/${sample}/${sample}.summary.reads.csv
echo "sample,$(echo ${primers} | sed 's/ /,/g'),merged" > ${outputDir}/${sample}/${sample}.summary.UMIs.csv
echo "sample,$(echo ${primers} | sed 's/ /,/g'),merged" > ${outputDir}/${sample}/${sample}.summary.UMIs.unique.csv
primersReads=()
primersUmis=()
primersUmisUniq=()
for primer in ${primers}; do
    primersReads+=$(echo "$(grep ">" ${outputDir}/${sample}/${sample}.${primer}.trimmed.fasta | wc -l | sed 's/ //g'),")
    primersUmis+=$(echo "$(cat ${outputDir}/${sample}/${sample}.${primer}.UMIs.length_filtered.csv | wc -l | sed 's/ //g'),")
    primersUmisUniq+=$(echo "$(cat ${outputDir}/${sample}/${sample}.${primer}.UMIs.length_filtered.frequency.csv | wc -l | sed 's/ //g'),")
done
mergedReads=$(grep ">" ${outputDir}/${sample}/${sample}.*.trimmed.fasta | wc -l) \
&& rm -f ${outputDir}/${sample}/${sample}.*.trimmed.fasta
mergedUmis=$(cat ${outputDir}/${sample}/${sample}.merged.UMIs.length_filtered.csv | wc -l)
mergedUmisUniq=$(cat ${outputDir}/${sample}/${sample}.merged.UMIs.length_filtered.frequency.csv | wc -l)
echo "${sample},${primersReads}${mergedReads// /}" >> ${outputDir}/${sample}/${sample}.summary.reads.csv
echo "${sample},${primersUmis}${mergedUmis// /}" >> ${outputDir}/${sample}/${sample}.summary.UMIs.csv
echo "${sample},${primersUmisUniq}${mergedUmisUniq// /}" >> ${outputDir}/${sample}/${sample}.summary.UMIs.unique.csv

# compress sequencing fasta file
echo "COMPRESSING FASTA FILE"

gzip ${outputDir}/${sample}/${sample}.fasta

# temporary lines, one-time use: move to /groups/
#echo "MOVING TO /GROUPS/"
#mv ${outputDir}/${sample} /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2
#chgrp -R prj_proudhon_lab /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2/${sample}
#chmod -R g+rws /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2/${sample}

echo "PREPROCESSING DONE"