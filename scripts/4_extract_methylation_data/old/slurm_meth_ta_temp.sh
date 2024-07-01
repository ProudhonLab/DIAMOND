#!/bin/bash

# Env

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/extract_methyl

# Arguments

baseDir=$0

inputAlign=$1
if [[ -d $inputAlign ]]; then
	allAligns=($(ls ${inputAlign}/inserts_dereplicated.length_filtered.*.align))
    inputAlign=${allAligns[${SLURM_ARRAY_TASK_ID}]}
	shDir=($(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | tr " " "\n"))
    baseDir=$(dirname ${shDir[0]})
else
    exit "Invalid argument. Expecting a directory."
fi
outputDir=$2
primer=$(echo ${inputAlign} | sed -nE 's/.*inserts_dereplicated.length_filtered\..*\.(.*)\.align$/\1/p')

echo "Checking arguments:"
echo "baseDir = ${baseDir}"
echo "inputAlign = ${inputAlign}"
echo "outputDir = ${outputDir}"
echo "primer = ${primer}"

# Script
python ${baseDir}/methylation_haplotypes_cluster_ver.py -a -p "${primer}" --cg_methyl --include_ta ${inputAlign} > ${outputDir}/allSamples.${primer}.methylation_vs_ta.csv