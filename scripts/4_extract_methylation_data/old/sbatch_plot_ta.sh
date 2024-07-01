#!/bin/bash

# Arguments

inputDir=$1
outputDir=$2

# list align files for each primer to be processed
allAligns=($(ls ${inputDir}/inserts_dereplicated.length_filtered.*.align))

# run sbatch for the dataset
# Arguments for extraction of methyl:
# $1 input directory
# $2 output directory
nbJobs=$(echo "${#allAligns[@]} - 1" | bc -l)
sbatch -W --cpus-per-task=10 --mem=50G -o ${outputDir}/logs/extract_methyl-%A_%a.out --array=0-${nbJobs}%20 scripts/3_extract_methylation_data/slurm_plot_ta.sh ${inputDir} ${outputDir}

chgrp -R prj_proudhon_lab ${outputDir}
chmod -R g+rws ${outputDir}