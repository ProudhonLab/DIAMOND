#!/bin/bash

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif

inputDir=$1
outputDir=$2
nRuns=5000
sklearnScript="/home/genouest/cnrs_umr6074/kdasilva/proudhon_lab/psl1_meth/scripts/5_classifier/old/healthy_vs_each_class/methyl_random_forest_sklearn_colors1.py"

allInputs=($(ls ${inputDir}/*))
input=${allInputs[${SLURM_ARRAY_TASK_ID}]}
input=$(basename ${input})
input=${input%.csv}

echo $input

time python ${sklearnScript} -r ${nRuns} -s -o ${outputDir}/${input}.with_subsampling.${nRuns}-runs_60-40 ${inputDir}/${input}.csv