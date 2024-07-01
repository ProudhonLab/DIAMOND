#!/bin/bash
#SBATCH --job-name=launcher


. /local/env/envconda.sh
conda activate /scratch/kvongrafenst/Nextflow


echo $mail
nextflow run  $2 $3 $4 $5 $6 $7 $8 $9 -N ${1} -with-timeline -with-report -resume
#sbatch baseLaunchingNextflow.sh klaus.von-grafenstein@univ-rennes.fr workflowPCRUnique.nf
