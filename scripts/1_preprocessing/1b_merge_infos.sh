#!/bin/bash

outputDir=$1
primersInfo=$2
samples=($(ls -d ${outputDir}/*/ | grep -v "logs/" | sed 's|[/]$||g'))
datasetName=$3

# Merge all summaries when all jobs from array are finished

i=0
for s in ${samples[@]}/; do
    sample=$(basename ${s})
    if [[ $i -eq 0 ]] ; then
        head -1  "${outputDir}/${sample}/${sample}.summary.csv" > "${outputDir}/${datasetName}.summary.csv"
        head -1  "${outputDir}/${sample}/${sample}.summary.reads.csv" > "${outputDir}/${datasetName}.summary.reads.csv"
        head -1  "${outputDir}/${sample}/${sample}.summary.UMIs.csv" > "${outputDir}/${datasetName}.summary.UMIs.csv"
        head -1  "${outputDir}/${sample}/${sample}.summary.UMIs.unique.csv" > "${outputDir}/${datasetName}.summary.UMIs.unique.csv"
    fi
	tail -n +2  "${outputDir}/${sample}/${sample}.summary.csv" >>  "${outputDir}/${datasetName}.summary.csv"
    tail -n +2  "${outputDir}/${sample}/${sample}.summary.reads.csv" >>  "${outputDir}/${datasetName}.summary.reads.csv"
    tail -n +2  "${outputDir}/${sample}/${sample}.summary.UMIs.csv" >>  "${outputDir}/${datasetName}.summary.UMIs.csv"
    tail -n +2  "${outputDir}/${sample}/${sample}.summary.UMIs.unique.csv" >>  "${outputDir}/${datasetName}.summary.UMIs.unique.csv"
	i=$(( $i + 1 ))
done

# Merge all inserts fasta

#awk -F',' '(NR!=1) { split($0, line); print line[1],line[2],line[3],line[4],line[5],line[6] }' ${primersInfo} | 
#while read -r primer forward reverse reverse_revcomp refsize withCG; do 
#    cat $(find ${outputDir} -name inserts_dereplicated.length_filtered.*.${primer}.fasta) > ${outputDir}/inserts_dereplicated.length_filtered.${datasetName}.${primer}.fasta
#done