#!/bin/bash

ListOfFastas=$1
inputFastaBasename=$2
outputDir=$3
threads=$4
nSeqsSample=20000000

mkdir -p ${outputDir}/clustering
outputDir=${outputDir}/clustering
cat ${ListOfFastas} | xargs cat > ${outputDir}/${inputFastaBasename}.fasta
inputFasta=${outputDir}/${inputFastaBasename}.fasta
subsampledInputFasta=${outputDir}/${inputFastaBasename}.subsampled.fasta
consensusFasta=${outputDir}/${inputFastaBasename}.consensus.fasta

# Check that there is enough reads for subsampling
nSeq=$(grep '>' ${inputFasta} | wc -l)

if [ "${nSeq}" -lt "${nSeqsSample}" ]; then
	# Use all sequences from input fasta
	echo -e "INFO - No subsampling\n"
	subsampledInputFasta=${inputFasta}
else
	# Subsample input fasta for clustering
	echo -e "INFO - Subsampling\n"
	cmd="vsearch --fastx_subsample ${inputFasta} --sample_size ${nSeqsSample} --fastaout ${subsampledInputFasta}"
	echo "${cmd}" 1>&2
	eval "${cmd} || exit 1"
fi

# Cluster
echo -e "INFO - Clustering\n"
cmd="vsearch --cluster_fast ${subsampledInputFasta} --notrunclabels --threads ${threads} --fasta_width 0 --iddef 4 --id 0 --qmask none --clusterout_sort --consout ${consensusFasta}"
#cmd="vsearch --cluster_size ${subsampledInputFasta} --notrunclabels --threads ${threads} --fasta_width 0 --iddef 4 --id 0 --qmask none --clusterout_sort --consout ${consensusFasta}"
echo "${cmd}" 1>&2
eval "${cmd} || exit 1"

# remove temporary file
rm ${outputDir}/${inputFastaBasename}.fasta
rm ${subsampledInputFasta}