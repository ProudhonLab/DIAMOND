#!/bin/bash

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/preproc

# Arguments and Variables

inputDir=$1
if [[ -d $inputDir ]]; then
    samples=($(ls -d ${inputDir}/*/ | sed 's|[/]$||g'))
    inputDir=${samples[${SLURM_ARRAY_TASK_ID}]}
    shDir=($(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | tr " " "\n"))
    baseDir=$(dirname ${shDir[0]})
elif [[ -f $inputDir ]]; then
    readarray -t samples < $inputDir
    inputDir=${samples[${SLURM_ARRAY_TASK_ID}]}
    shDir=($(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | tr " " "\n"))
    baseDir=$(dirname ${shDir[0]})
else
    exit "Invalid argument. Expecting a directory or a text file listing directories."
fi

${baseDir}/1a_sample_preprocessing.sh ${inputDir} $2 $3 $4 $5 $6