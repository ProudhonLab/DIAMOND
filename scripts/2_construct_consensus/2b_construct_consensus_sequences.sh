#!/bin/bash

consensusFasta=$1
outputDir=$2
mafft_matrix=$3
nRefSeqs=$4

mkdir -p ${outputDir}/${nRefSeqs}_largest
outputDir=${outputDir}/${nRefSeqs}_largest
inputFastaBasename=$(basename ${consensusFasta%".consensus.fasta"})
nLargestClustersFasta=${outputDir}/${inputFastaBasename}.consensus.${nRefSeqs}_largest.fasta
multipleAlignFasta=${outputDir}/${inputFastaBasename}.consensus.${nRefSeqs}_largest.msa.fasta

# Isolate consensus sequences from <nRefSeqs> largest clusters
echo -e "INFO - Get consensus sequences\n"
cmd="awk \"/^>/ {n++} n>${nRefSeqs} {exit} 1\" ${consensusFasta} > ${nLargestClustersFasta}"
echo "${cmd}" 1>&2
eval "${cmd} || exit 1"

# Align consensus sequences using a methylation-focused score matrix
echo -e "INFO - Aligning consensus sequences\n"
if [ "${mafft_matrix}" == "None" ]; then
    #cmd="mafft --retree 2 ${nLargestClustersFasta} > ${multipleAlignFasta}"
    cmd="mafft --globalpair --maxiterate 1000 ${nLargestClustersFasta} > ${multipleAlignFasta}"
else
    cmd="mafft --textmatrix ${mafft_matrix} --retree 2 ${nLargestClustersFasta} > ${multipleAlignFasta}"
fi
echo "${cmd}" 1>&2
eval "${cmd} || exit 1"