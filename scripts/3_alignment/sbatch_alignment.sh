#!/bin/bash

# Arguments

inputList=$1
clust_template_folder=$2
outputDir=$3

mkdir -p $outputDir
mkdir -p ${outputDir}/logs

# list fasta files for each primer to be processed
readarray -t fastas < $inputList

# run sbatch for the dataset
# Arguments for clustering and alignment:
# $1 inputFasta
# $2 output directory
# $3 custom mafft
nbJobs=$(echo "${#fastas[@]} - 1" | bc -l)
sbatch --cpus-per-task=8 --mem=20G -o ${outputDir}/logs/alignment-%A_%a.out --array=0-${nbJobs}%24 scripts/3_alignment/slurm_alignment.sh ${inputList} ${clust_template_folder} ${outputDir}

#chgrp -R prj_proudhon_lab ${outputDir}
#chmod -R g+rws ${outputDir}