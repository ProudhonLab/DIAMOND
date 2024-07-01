#!/bin/bash

inputPrefix=$1
inputFastaBasename=$2
outputDir=$3
threads=$4

# run sbatch for the dataset
mkdir -p ${outputDir}/logs
sbatch -W --cpus-per-task=${threads} --mem=50G -o ${outputDir}/logs/clustering-%A_%a.out --array=0-7%8 scripts/2_construct_consensus/slurm_clustering.sh ${inputPrefix} ${inputFastaBasename} ${outputDir} ${threads}
