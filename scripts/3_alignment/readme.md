# Alignment scripts

## No main script

Contrary to the other steps, the alignment does not require a main script.
mothur is used directly in the slurm and sbatch scripts.

## slurm script

The slurm script is only used to call mothur when using the sbatch script (see below). The slurm script will process only one fasta file based on the task ID (job array).

## sbatch script

The sbatch script allows to run the alignment for multiple fastas in parallel. The arguments are:
1. A text file listing the paths of the fasta files to align.
2. The path to the folder containing the reference/template to align to (generated during step 2).
3. The output directory.

# Example (user=kdasilva)

Align all available samples.
Using default_scores and 10_largest clusters.

## Prepare list of fastas to align

```
find /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2/ -name inserts_dereplicated.length_filtered.*.fasta > /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastas_to_align.txt; find /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_merge/ -name inserts.length_filtered.*.fasta >> /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastas_to_align.txt
```

## Run Alignment

All fastas in parallel.
```
cd ~/proudhon_lab/psl1_meth/

time scripts/3_alignment/sbatch_alignment.sh /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastas_to_align.txt /scratch/kdasilva/l1pameth_2023/consensus_sequences/default_scores/10_largest /scratch/kdasilva/l1pameth_2023/alignments/default_scores/10_largest
```

Time = 119m42.055s - 334m

## Check for empty align

```
find /scratch/kdasilva/l1pameth_2023/alignments/default_scores/10_largest -type f -empty -name "*.align" > /scratch/kdasilva/l1pameth_2023/align_failed_summary.txt
```

No empty alignments.
Total: 1136 samples for each primer, except L1HS-15, only 1126 samples.
