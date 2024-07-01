#!/bin/bash

# Env

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/extract_methyl2

# Arguments

inputList=$1
outputDir=$2
samples_info=$3

sample_array=($(cat ${inputList} | sed 's/\.L1HS_.*$//' | uniq))
inserts_prefix=${sample_array[${SLURM_ARRAY_TASK_ID}]}
sample_id=${inserts_prefix##*.}


shDir=($(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | tr " " "\n"))
baseDir=$(dirname ${shDir[0]})

echo "ARGUMENTS USED ARE:"
echo "- inputList: ${inputList}"
echo "- outputDir: ${outputDir}"
echo "- samples_info: ${samples_info}"
echo -e "Also, the sample used for this job is: ${sample_id}.\n"

# Script
listPrimers=(L1HS_2 L1HS_14 L1HS_15 L1HS_4 L1HS_8 L1HS_9 L1HS_7 L1HS_6)
for primer in "${listPrimers[@]}"; do 
    selectedPos=/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/selected_pos_${primer}.pickle
    pypy ${baseDir}/4b_extract_methyl_haplo.py -i ${inserts_prefix}.${primer}.align -o ${outputDir} -s ${samples_info} -c ${selectedPos} -m all
    if [ "${primer}" == "L1HS_2" ]; then
        cat ${outputDir}/cg_methyl.${sample_id}.${primer}.csv > ${outputDir}/cg_methyl.${sample_id}.csv
        cat ${outputDir}/haplotypes.${sample_id}.${primer}.csv > ${outputDir}/haplotypes.${sample_id}.csv
    else
        fname=$(mktemp) && cut --complement -d, -f1,2 ${outputDir}/cg_methyl.${sample_id}.${primer}.csv | paste -d, ${outputDir}/cg_methyl.${sample_id}.csv - >> "$fname" && mv "$fname" ${outputDir}/cg_methyl.${sample_id}.csv
        fname=$(mktemp) && cut --complement -d, -f1,2 ${outputDir}/haplotypes.${sample_id}.${primer}.csv | paste -d, ${outputDir}/haplotypes.${sample_id}.csv - >> "$fname" && mv "$fname" ${outputDir}/haplotypes.${sample_id}.csv
    fi
    rm ${outputDir}/cg_methyl.${sample_id}.${primer}.csv
    rm ${outputDir}/haplotypes.${sample_id}.${primer}.csv
done

