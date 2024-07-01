#!/bin/bash

inputSamples=$1 # can be a directory or a list of sample directories
outputDir=$2 # usually the scratch folder. A $sample folder will be created inside.
strand=$3 # R1, R2, or merge
primersInfo=$4 # 6 columns: primer name, seq forward, seq reverse, seq reverse complement, ref size, with CG boolean
umisSize=$5 # 16
threads=$6 # 12

# Check
if [ $# -ne 6 ]; then
    exit "No arguments or not enough supplied."
fi

# If directory, get list of sub-directories.
if [[ -d $inputSamples ]]; then
    echo "The input samples is a directory."
    samples=($(ls -d ${inputSamples}/*/ | sed 's|[/]$||g'))
elif [[ -f $inputSamples ]]; then
    echo "The input samples is a txt file."
    readarray -t samples < $inputSamples
else
    exit "Invalid argument. Expecting a directory or a text file listing directories."
fi

# run sbatch for the dataset
# Arguments for preprocessing:
# $1 inputDir
# $2 outputDir
# $3 strand: R1, R2, or merge
# $4 primersInfo csv (6 columns: primer name, seq forward, seq reverse, seq reverse complement, ref size, with CG boolean)
# $5 umisSize (16)
# $6 threads (12)
mkdir -p ${outputDir}/logs
nbJobs=$(echo "${#samples[@]} - 1" | bc -l)
sbatch -W --cpus-per-task=${threads} --mem=50G -o ${outputDir}/logs/preproc-%A_%a.out --array=0-${nbJobs}%24 scripts/1_preprocessing/slurm_sample_preprocessing.sh ${inputSamples} ${outputDir} ${strand} ${primersInfo} ${umisSize} ${threads}

#scripts/1_sequencing_batch_pipeline/1b_merge_infos.sh ${outputDir} ${primersInfo}