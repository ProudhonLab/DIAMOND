#!/bin/bash

# Arguments

inputDir=$1
outputDir=$2
mkdir -p ${outputDir}/logs

allInputs=($(ls ${inputDir}/*))
nbJobs=$(echo "${#allInputs[@]} - 1" | bc -l)

# run sbatch
sbatch --mem=20G -o ${outputDir}/logs/classif-%A_%a.out --array=0-${nbJobs}%24 scripts/5_classifier/old/healthy_vs_each_class/slurm_classif_subsampling.sh ${inputDir} ${outputDir}
# sbatch --mem=20G -o ${outputDir}/logs/classif-%A_%a.out --array=0-${nbJobs}%24 scripts/5_classifier/old/healthy_vs_each_class/slurm_classif_no_subsampling.sh ${inputDir} ${outputDir}