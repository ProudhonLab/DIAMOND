#!/bin/bash

# Arguments

inputList=$1
outputDir=$2
samples_info=$3
threshold=$4

mkdir -p ${outputDir}/logs

# run sbatch, primers in parallel
sbatch -W --mem=50G -o ${outputDir}/logs/selectPositions-%A_%a.out --array=0-7%8 scripts/4_extract_methylation_data/slurm_select_positions.sh ${inputList} ${outputDir} ${samples_info} ${threshold}

chgrp -R prj_proudhon_lab ${outputDir}
chmod -R g+rws ${outputDir}