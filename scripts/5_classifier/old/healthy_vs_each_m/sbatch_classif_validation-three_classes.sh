#!/bin/bash

# Arguments

inputDir=$1
outputDir=$2
mkdir -p ${outputDir}/logs

allInputs=($(ls ${inputDir}/*))
nbJobs=$(echo "${#allInputs[@]}" | bc -l)

# run sbatch
sbatch --mem=20G -o ${outputDir}/logs/classif-%A_%a.out --array=0-${nbJobs}%24 scripts/5_classifier/slurm_classif_validation_subsampling-three_classes.sh ${inputDir} ${outputDir}
sbatch --mem=20G -o ${outputDir}/logs/classif-%A_%a.out --array=0-${nbJobs}%24 scripts/5_classifier/slurm_classif_validation_no_subsampling-three_classes.sh ${inputDir} ${outputDir}