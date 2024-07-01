#!/bin/bash

# Env

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/extract_methyl2

# Arguments

inputList=$1
outputDir=$2
samples_info=$3
threshold=$4

listPrimers=(L1HS_2 L1HS_4 L1HS_6 L1HS_7 L1HS_8 L1HS_9 L1HS_14 L1HS_15)
primer=${listPrimers[${SLURM_ARRAY_TASK_ID}]}

grep ".*${primer}\.align$" ${inputList} > ${outputDir}/tmp_list_${primer}.txt
inputList=${outputDir}/tmp_list_${primer}.txt

shDir=($(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | tr " " "\n"))
baseDir=$(dirname ${shDir[0]})

echo "ARGUMENTS USED ARE:"
echo "- inputList: ${inputList}"
echo "- outputDir: ${outputDir}"
echo "- samples_info: ${samples_info}"
echo "- threshold: ${threshold}"
echo "Also, the primer used for this job is: ${primer}."

# Script 

pypy ${baseDir}/4a_select_positions.py -i ${inputList} -o ${outputDir} -s ${samples_info} -p ${primer} -t ${threshold}

rm ${outputDir}/tmp_list_${primer}.txt