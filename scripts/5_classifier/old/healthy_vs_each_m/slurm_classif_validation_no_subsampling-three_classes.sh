#!/bin/bash

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif

inputDir=$1
outputDir=$2
nRuns=5000
sklearnScript="/home/genouest/cnrs_umr6074/kdasilva/proudhon_lab/psl1_meth/scripts/5_classifier/methyl_random_forest_sklearn_three_classes.py"

allInputs=($(ls ${inputDir}/*))
input=${allInputs[${SLURM_ARRAY_TASK_ID}]}
training="${input/validation/discovery}"
training="${training/validation_cohort/discovery_cohort}"
inputBase=$(basename ${input})
inputBase=${inputBase%.csv}

echo $training
echo $input
echo $inputBase

time python ${sklearnScript} -r ${nRuns} ${training} --test ${input} -o ${outputDir}/${inputBase}.no_subsampling.${nRuns}-runs