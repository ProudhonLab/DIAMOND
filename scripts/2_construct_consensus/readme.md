# Clustering scripts

## Main script (2a)

### Principle

Reads are clustered using vsearch.
Undersampling if there are too much reads.

### Arguments

Positional arguments:
1. A text file listing the paths of the fastas to cluster
2. Path and prefix of the output file
3. Output directory
4. Number of threads
5. Limit for undersampling [20000000]

## slurm script

The slurm script is only used to call the main script when using the sbatch script (see below). The slurm script will process only reads from one primer based on the task ID (job array).

## sbatch script

The sbatch script allows to run the main script for multiple primers in parallel. The arguments are the same as the main script except the limit for undersampling stays by default.

# Consensus construction script (2b)

## Principle

Isolate representative sequences of the clusters.
Multiple alignment of those representative sequences.

## Arguments

Positional arguments:
1. consensus.fasta file previously generated
2. output directory
3. Mafft matrix used for the alignment scores [None or custom score matrix]
4. Number of representative sequences [currently 10 for DIAMOND]

# Other scripts

- extract_clusters_nb_seqs.sh
- extract_clusters_nb_seqs.R
are scripts to analyze the composition of the clusters. During development, we were wondering if we should select the number of representative sequences according to the number of sequences in the clusters.

- generate_list_of_fastas_to_cluster.R
generates the text file to give in input to the 2a clustering script.

# Example (user=kdasilva)

Clustering using C1/C2 samples (only plasma and cleaned from duplicates).

```
mamba create -p /groups/proudhon_lab/projects/l1pa_meth/conf/clust_align -c conda-forge -c bioconda boost mothur vsearch mafft
```

## Prepare list of fastas to cluster

See R script generate_list_of_fastas_to_cluster.R

## Run clustering

All primers in parallel.
```
cd ~/proudhon_lab/psl1_meth/

time scripts/2_construct_consensus/sbatch_clustering.sh /groups/proudhon_lab/projects/l1pa_meth/0_data/fastas_to_cluster_c1c2 c1c2 /scratch/kdasilva/l1pameth_2023/consensus_sequences 12
```

Time = 27m22.277s

## Extract consensus sequences and MSA

Test with/without custom score matrix.
Test 1, 2, 3, 4, 5, 10 largest clusters.

```
for scores in 0 1;do
	if [ ${scores} -eq 0 ]; then
		mafft_matrix=None
		mkdir /scratch/kdasilva/l1pameth_2023/consensus_sequences/default_scores
		outputDir=/scratch/kdasilva/l1pameth_2023/consensus_sequences/default_scores
	else
		mafft_matrix=scripts/2_construct_consensus/mafft_custom_score_matrix.txt
		mkdir /scratch/kdasilva/l1pameth_2023/consensus_sequences/custom_scores
		outputDir=/scratch/kdasilva/l1pameth_2023/consensus_sequences/custom_scores
	fi
	for nseqs in 1 2 3 4 5 10; do
		for primers in 2 4 6 7 8 9 14 15; do
			scripts/2_construct_consensus/2b_construct_consensus_sequences.sh /scratch/kdasilva/l1pameth_2023/consensus_sequences/clustering/c1c2_L1HS_${primers}.consensus.fasta ${outputDir} ${mafft_matrix} ${nseqs}
		done
	done
done
```

One line version:
```
. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/clust_align

for scores in 0 1;do if [ ${scores} -eq 0 ]; then mafft_matrix=None;mkdir /scratch/kdasilva/l1pameth_2023/consensus_sequences/default_scores; outputDir=/scratch/kdasilva/l1pameth_2023/consensus_sequences/default_scores;	else mafft_matrix=scripts/2_construct_consensus/mafft_custom_score_matrix.txt; mkdir /scratch/kdasilva/l1pameth_2023/consensus_sequences/custom_scores; outputDir=/scratch/kdasilva/l1pameth_2023/consensus_sequences/custom_scores;fi; for nseqs in 1 2 3 4 5 10; do for primers in 2 4 6 7 8 9 14 15; do scripts/2_construct_consensus/2b_construct_consensus_sequences.sh /scratch/kdasilva/l1pameth_2023/consensus_sequences/clustering/c1c2_L1HS_${primers}.consensus.fasta ${outputDir} ${mafft_matrix} ${nseqs}; done; done; done
```