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
sbatch -W -o ${outputDir}/logs/extract_methyl-%A_%a.out --array=0-${nbJobs}%20 scripts/3_extract_methylation_data/slurm_methylation_haplotypes_cluster_ver.sh ${inputDir} ${outputDir}