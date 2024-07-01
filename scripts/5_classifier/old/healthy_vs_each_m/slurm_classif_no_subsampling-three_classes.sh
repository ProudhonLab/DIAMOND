#!/bin/bash

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif

inputDir=$1
outputDir=$2
nRuns=5000
sklearnScript="/home/genouest/cnrs_umr6074/kdasilva/proudhon_lab/psl1_meth/scripts/5_classifier/methyl_random_forest_sklearn_three_classes.py"

allInputs=($(ls ${inputDir}/*))
input=${allInputs[${SLURM_ARRAY_TASK_ID}]}
input=$(basename ${input})
input=${input%.csv}

echo $input

time python ${sklearnScript} -r ${nRuns} -o ${outputDir}/${input}.no_subsampling.${nRuns}-runs_60-40 ${inputDir}/${input}.csv