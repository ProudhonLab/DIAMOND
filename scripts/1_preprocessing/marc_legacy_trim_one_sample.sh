#!/bin/bash
#SBATCH --job-name=marc_legacy_trim_not_nextflow

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/l1pa_meth/Conf/marc_legacy


cd /groups/proudhon_lab/l1pa_meth

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/l1pa_meth/Conf/marc_legacy

dataDir="/groups/proudhon_lab/l1pa_meth/Cna/raw/seq"
rootDir="/groups/proudhon_lab/l1pa_meth/klausTest"
scriptsDir="/groups/proudhon_lab/l1pa_meth/scripts"
samples= "A1146S20"
## Samples pre- and suffixes will be determined automatically, variables here are for reference:
# samplesPrefix="D541S"
# samplesSuffixes=$(echo {1..57})
# samplesSuffixes=$(seq -f "%03g" 1 164)
umisSize=16
# Strand can be either: "R1" (R1 only), "R2" (R2 only) or "merge" (merge R1 & R2)
# strand="R2"
strand="R2"
## Illumina tag will be determined automatically, variables here are for reference:
# illuminaTag="SN7001339"
# illuminaTag="A00514"
sample="A1146S20"
threads=12

primersList="/groups/proudhon_lab/l1pa_meth/data/primers_MHOVC_dataset_2010012.csv"
primersRefSizes="/groups/proudhon_lab/l1pa_meth/data/primers_MHOVC_dataset_2010012_refsizes.csv"
primers=$(sed 1d ${primersList} | cut -f 1)
# List of primer sets that contain CG in their sequences
primersWithCG=( L1HS_2 L1HS_4 L1HS_14 L1HS_15 )

mkdir -p ${rootDir}
echo ${strand} > ${rootDir}/strand_used_for_analysis.txt

## Filter reads by quality using fastp
mkdir -p ${rootDir}/${sample}
if [ ${strand} = merge ]; then
    fastp -i ${dataDir}/${sample}/${sample}.R1.fastq.gz -I ${dataDir}/${sample}/${sample}.R2.fastq.gz -o ${rootDir}/${sample}/${sample}.R1.fastq -O ${rootDir}/${sample}/${sample}.R2.fastq -h ${rootDir}/${sample}/fastp_${sample}.html -j ${rootDir}/${sample}/fastp_${sample}.json
elif [ ${strand} = R1 ]; then
    fastp -i ${dataDir}/${sample}/${sample}.R1.fastq.gz -o ${rootDir}/${sample}/${sample}.fastq -h ${rootDir}/${sample}/fastp_${sample}.html -j ${rootDir}/${sample}/fastp_${sample}.json
elif [ ${strand} = R2 ]; then
    vsearch --fastx_revcomp ${dataDir}/${sample}/${sample}.R2.fastq.gz -fastqout ${rootDir}/${sample}/${sample}.reversed.fastq
    fastp -i ${rootDir}/${sample}/${sample}.reversed.fastq -o ${rootDir}/${sample}/${sample}.fastq -h ${rootDir}/${sample}/fastp_${sample}.html -j ${rootDir}/${sample}/fastp_${sample}.json
fi

## Replace illumina device tag by sample name
if [ ${strand} = merge ]; then
    illuminaTag=$(head -n 1 ${rootDir}/${sample}/${sample}.R1.fastq | grep -P -o "@(.+?):" | tr -d '@:')
    sed -i "s/${illuminaTag}/${sample}/g" ${rootDir}/${sample}/${sample}.R1.fastq
    sed -i "s/${illuminaTag}/${sample}/g" ${rootDir}/${sample}/${sample}.R2.fastq
else
    illuminaTag=$(head -n 1 ${rootDir}/${sample}/${sample}.fastq | grep -P -o "@(.+?):" | tr -d '@:')
    sed -i "s/${illuminaTag}/${sample}/g" ${rootDir}/${sample}/${sample}.fastq
fi

## Merge if needed
if [ ${strand} = merge ]; then
    cd ${rootDir}/${sample}
    cmd="mothur '#make.contigs(ffastq=${sample}.R1.fastq, rfastq=${sample}.R2.fastq, trimoverlap=T, processors=${threads})'"
    eval ${cmd}
    # Strip Mothur's error rate output
    sed -i "s/\tee=.*$//g" ${sample}.R1.trim.contigs.fasta
    mv ${sample}.R1.trim.contigs.fasta ${sample}.fasta
fi

cd ${rootDir}

if [ ${strand} != merge ]; then
	# Make fasta out of fastq
    paste - - - - < ${sample}/${sample}.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${sample}/${sample}.fasta
fi

if [ ${strand} = merge ]; then
	# Plot reads length after merge
    cat ${sample}/${sample}.fasta | grep -v ">" | awk '{print length}' > ${sample}/${sample}.reads_length.csv
fi

## Cut adapters
mkdir -p ${sample}/forward_adapter_trimmed ${sample}/reverse_adapter_trimmed
if [ ${strand} = R1 ]; then
    mkdir -p ${sample}/illumina_adapter_trimmed
fi

## Trim adapters
while read Primer Forward Reverse Reverse_revcomp; do
    atropos trim -g ${Forward} -se ${sample}/${sample}.fasta -o ${sample}/forward_adapter_trimmed/${sample}.${Primer}.fasta --threads ${threads} -y ":${Primer}" --discard-untrimmed --no-default-adapters --overlap 7 --info-file=${sample}/forward_adapter_trimmed/${sample}.${Primer}.info.csv && grep -v "	-1	" ${sample}/forward_adapter_trimmed/${sample}.${Primer}.info.csv > ${sample}/forward_adapter_trimmed/${sample}.${Primer}.info.minus_one_filtered.csv && rm -f ${sample}/forward_adapter_trimmed/${sample}.${Primer}.info.csv
done < <(sed 1d ${primersList})

while read Primer Forward Reverse Reverse_revcomp; do
    atropos trim -a ${Reverse_revcomp} -se ${sample}/forward_adapter_trimmed/${sample}.${Primer}.fasta -o ${sample}/reverse_adapter_trimmed/${sample}.${Primer}.fasta --threads ${threads} -y ":${Primer}" --discard-untrimmed --no-default-adapters --overlap 7 --info-file=${sample}/reverse_adapter_trimmed/${sample}.${Primer}.info.csv && grep -v "	-1	" ${sample}/reverse_adapter_trimmed/${sample}.${Primer}.info.csv > ${sample}/reverse_adapter_trimmed/${sample}.${Primer}.info.minus_one_filtered.csv && rm -f ${sample}/reverse_adapter_trimmed/${sample}.${Primer}.info.csv
done < <(sed 1d ${primersList})

# Not needed at this time, stays here just in case it might be needed later
if [ ${strand} = NORUN ]; then
	# Cut illumina CS1/2 adapters to recover UMIs
	# Create fasta file from info file: whatever follows reverse primer, e.g. UMIs + partial CS2
    for primer in ${primers}; do
        python ${scriptsDir}/info_file_to_fasta.py ${sample}/reverse_adapter_trimmed/${sample}.${primer}.info.minus_one_filtered.csv --primer ${primer} --end reverse_umi > ${sample}/reverse_adapter_trimmed/UMIs_CS2.${sample}.${primer}.fasta
    done

    while read Primer Forward Reverse Reverse_revcomp; do
        # This is part of illumina CS1/2 adapters sequence
        atropos trim -a AGACCAAGTCTCTGCTACCGTA -se ${sample}/reverse_adapter_trimmed/UMIs_CS2.${sample}.${Primer}.fasta -o ${sample}/reverse_adapter_trimmed/UMIs.${sample}.${Primer}.fasta --threads ${threads} --discard-untrimmed --no-default-adapters --overlap 5 --info-file=${sample}/illumina_adapter_trimmed/${sample}.${Primer}.info.csv && grep -v "	-1	" ${sample}/illumina_adapter_trimmed/${sample}.${Primer}.info.csv > ${sample}/illumina_adapter_trimmed/${sample}.${Primer}.info.minus_one_filtered.csv && rm -f ${sample}/illumina_adapter_trimmed/${sample}.${Primer}.info.csv
    done < <(sed 1d ${primersList})
fi

## Isolate UMIs from merged reads
if [ ${strand} = merge ] || [ ${strand} = R1 ]; then
    for primer in ${primers}; do
        python ${scriptsDir}/info_file_to_fasta.py ${sample}/reverse_adapter_trimmed/${sample}.${primer}.info.minus_one_filtered.csv --primer ${primer} --end reverse_umi > ${sample}/reverse_adapter_trimmed/UMIs.${sample}.${primer}.fasta
    done
fi

## Plot inserts size
mkdir -p reads_length_plots

for primer in ${primers}; do
    grep -v ">" ${sample}/reverse_adapter_trimmed/${sample}.${primer}.fasta | awk '{print length}' > ${sample}/reverse_adapter_trimmed/reads_length.${sample}.${primer}.csv
done

while read primer refsize; do
    Rscript ${scriptsDir}/plot_read_or_uid_length.R -t "Reads length ${sample}, primer ${primer} - exp. size = ${refsize}" -o reads_length_plots/reads_length.${sample}.${primer}.png ${sample}/reverse_adapter_trimmed/reads_length.${sample}.${primer}.csv
done < <(sed 1d ${primersRefSizes})

## Filter inserts by size
while read primer refsize; do
    vsearch --fastx_filter ${sample}/reverse_adapter_trimmed/${sample}.${primer}.fasta --notrunclabels --fastaout ${sample}/reverse_adapter_trimmed/inserts.length_filtered.${sample}.${primer}.fasta --fastq_minlen $(echo "${refsize} - 5" | bc) --fastq_maxlen $(echo "${refsize} + 5" | bc)
done < <(sed 1d ${primersRefSizes})

if [ ${strand} = R2 ]; then
# Recreate fasta file from info file: reads id + UMIs
    for primer in ${primers}; do
        python ${scriptsDir}/info_file_to_fasta.py ${sample}/reverse_adapter_trimmed/${sample}.${primer}.info.minus_one_filtered.csv -p ${primer} --end reverse_umi > ${sample}/reverse_adapter_trimmed/UMIs.${sample}.${primer}.fasta
    done
fi

## Filter UMIs by size
for primer in ${primers}; do
    vsearch --fastx_filter ${sample}/reverse_adapter_trimmed/UMIs.${sample}.${primer}.fasta --notrunclabels --fastaout ${sample}/reverse_adapter_trimmed/UMIs.length_filtered.${sample}.${primer}.fasta --fastq_minlen ${umisSize} --fastq_maxlen ${umisSize}
done

## Make sure a file exists for all primers, even empty
for primer in ${primers}; do
    touch ${sample}/reverse_adapter_trimmed/UMIs.length_filtered.${sample}.${primer}.fasta
done

## Extract length filtered UMIs to a csv file
for primer in ${primers}; do
    cat ${sample}/reverse_adapter_trimmed/UMIs.length_filtered.${sample}.${primer}.fasta | grep -v ">" > ${sample}/reverse_adapter_trimmed/${sample}.${primer}.UMIs.length_filtered.csv
done

# Remove merge file from potential previous runs
rm -f ${sample}/reverse_adapter_trimmed/${sample}.merged.UMIs.length_filtered.csv
for primer in ${primers}; do
    cat ${sample}/reverse_adapter_trimmed/${sample}.${primer}.UMIs.length_filtered.csv >> ${sample}/reverse_adapter_trimmed/${sample}.merged.UMIs.length_filtered.csv
done

mkdir -p UMIs_frequencies_count_plots/${sample}

Rscript ${scriptsDir}/adapter_stats.R -t "UMIs frequencies count ${sample} - all primers" -o UMIs_frequencies_count_plots/${sample}.merged.png ${sample}/reverse_adapter_trimmed/${sample}.merged.UMIs.length_filtered.csv

for primer in ${primers}; do
    Rscript ${scriptsDir}/adapter_stats.R -t "UMIs frequencies count ${sample} - ${primer}" -o UMIs_frequencies_count_plots/${sample}/${sample}.${primer}.png ${sample}/reverse_adapter_trimmed/${sample}.${primer}.UMIs.length_filtered.csv
done


## Make sure a file exists for all primers, even empty
for primer in ${primers}; do
    touch ${sample}/reverse_adapter_trimmed/${sample}.${primer}.UMIs.length_filtered.frequency.csv
done

## Deduplication based on UMIs
# Concatenate length filtered inserts with length filtered UMIs
for primer in ${primers}; do
    python ${scriptsDir}/concatenate_UMIs_targets.py ${rootDir}/${sample}/reverse_adapter_trimmed/UMIs.length_filtered.${sample}.${primer}.fasta ${rootDir}/${sample}/reverse_adapter_trimmed/inserts.length_filtered.${sample}.${primer}.fasta > ${rootDir}/${sample}/reverse_adapter_trimmed/inserts_UMIs_concat.length_filtered.${sample}.${primer}.fasta
done

# Deduplicate
for primer in ${primers}; do
    vsearch --notrunclabels --derep_fulllength ${sample}/reverse_adapter_trimmed/inserts_UMIs_concat.length_filtered.${sample}.${primer}.fasta --sizein --output ${sample}/reverse_adapter_trimmed/inserts_UMIs_concat_dereplicated.length_filtered.${sample}.${primer}.fasta
done

# Trim UMIs after deduplication (UMIs are appended to reads ids)
for primer in ${primers}; do
    python ${scriptsDir}/trim_UMIs.py ${rootDir}/${sample}/reverse_adapter_trimmed/inserts_UMIs_concat_dereplicated.length_filtered.${sample}.${primer}.fasta > ${rootDir}/${sample}/reverse_adapter_trimmed/inserts_UMIs_concat_dereplicated.length_filtered.UMIs_trimmed.${sample}.${primer}.fasta
done

## Recover CG in primer sequences
for primer in ${primersWithCG[@]}; do
    # First isolate forward primers in a separate fasta file
    python ${scriptsDir}/info_file_to_fasta.py ${sample}/forward_adapter_trimmed/${sample}.${primer}.info.minus_one_filtered.csv --primer ${primer} --end forward_primer > ${sample}/forward_adapter_trimmed/forward_primers.${sample}.${primer}.fasta
    # Second isolate reverse primers in a separate fasta file
    python ${scriptsDir}/info_file_to_fasta.py ${sample}/reverse_adapter_trimmed/${sample}.${primer}.info.minus_one_filtered.csv --primer ${primer} --end reverse_primer > ${sample}/reverse_adapter_trimmed/reverse_primers.${sample}.${primer}.fasta
    # Third concatenate forward primers with target sequences, in 5'
    python ${scriptsDir}/concatenate_UMIs_targets.py --end 3prime ${rootDir}/${sample}/forward_adapter_trimmed/forward_primers.${sample}.${primer}.fasta ${rootDir}/${sample}/reverse_adapter_trimmed/inserts_UMIs_concat_dereplicated.length_filtered.UMIs_trimmed.${sample}.${primer}.fasta > ${rootDir}/${sample}/reverse_adapter_trimmed/forward_primers.inserts_UMIs_concat_dereplicated.length_filtered.UMIs_trimmed.${sample}.${primer}.fasta
    # Lastly concatenate reverse primers with forward primers + target sequences, in 3'
    python ${scriptsDir}/concatenate_UMIs_targets.py --end 5prime ${rootDir}/${sample}/reverse_adapter_trimmed/reverse_primers.${sample}.${primer}.fasta ${rootDir}/${sample}/reverse_adapter_trimmed/forward_primers.inserts_UMIs_concat_dereplicated.length_filtered.UMIs_trimmed.${sample}.${primer}.fasta > ${rootDir}/${sample}/reverse_adapter_trimmed/inserts_UMIs_concat_dereplicated.length_filtered.UMIs_trimmed.including_primers.${sample}.${primer}.fasta
done

