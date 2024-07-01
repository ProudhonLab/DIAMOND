#!/bin/bash

# Env

. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/clust_align

# Arguments

inputList=$1
clust_template_folder=$2
outputDir=$3

readarray -t fastas < $inputList
fasta_to_process=${fastas[${SLURM_ARRAY_TASK_ID}]}
primer=$(echo ${fasta_to_process} | awk -F. '{print $(NF-1)}')
clust_template=$(ls ${clust_template_folder}/*${primer}*msa.fasta) # ls should give only one file !

# Script

mothur "#set.dir(output=${outputDir});align.seqs(candidate=${fasta_to_process}, template=${clust_template}, align=needleman, match=1, mismatch=-1, gapopen=-1, gapextend=0)"

# make sure an align file exist, even if empty
align_name=$(basename ${fasta_to_process%.fasta})
touch ${outputDir}/${align_name}.align