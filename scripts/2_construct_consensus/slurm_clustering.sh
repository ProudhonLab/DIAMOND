#!/bin/bash

# Env

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/clust_align

# Arguments

inputPrefix=$1
listPrimers=(L1HS_2 L1HS_4 L1HS_6 L1HS_7 L1HS_8 L1HS_9 L1HS_14 L1HS_15)
primer=${listPrimers[${SLURM_ARRAY_TASK_ID}]}
fastasToCluster=${inputPrefix}_${primer}.txt

inputFastaBasename=$2
inputFastaBasename=${inputFastaBasename}_${primer}

outputDir=$3

threads=$4

shDir=($(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | tr " " "\n"))
baseDir=$(dirname ${shDir[0]})

# Script A = clustering 

${baseDir}/2a_clustering.sh ${fastasToCluster} ${inputFastaBasename} ${outputDir} ${threads}