#!/bin/bash

# Arguments

inputList=$1
outputDir=$2
samples_info=$3

mkdir -p ${outputDir}/logs

sample_array=($(cat ${inputList} | sed 's/\.L1HS_.*$//' | uniq))
nbJobs=$(echo "${#sample_array[@]} - 1" | bc -l)

# run sbatch, primers in parallel
sbatch --mem=10G -o ${outputDir}/logs/extractMethyl-%A_%a.out --array=0-${nbJobs}%24 scripts/4_extract_methylation_data/slurm_extract_methyl_haplo.sh ${inputList} ${outputDir} ${samples_info}

#chgrp -R prj_proudhon_lab ${outputDir}
#chmod -R g+rws ${outputDir}