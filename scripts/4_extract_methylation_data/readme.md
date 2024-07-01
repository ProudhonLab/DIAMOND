# Position selection

## Main script (4a)

### Principle

Alignments are read to determine the CG positions in the sequences, based on a threshold of CG+TG ratio.
Usually, only healthy samples are used for the selection.

### Arguments

Named arguments:
- -i: A text file listing the paths of alignment files to process (usually healthy donors data)
- -o: Output directory
- -s: The path to the samples_infos.csv (generated at step 0)
- -p: Name of the primer processed
- -t: Threshold of (CG+TG)/NN

## slurm script

The slurm script is only used to call the main script when using the sbatch script (see below). The slurm script will process only alignment files from one primer based on the task ID (job array).

## sbatch script

The sbatch script allows to run the main script for multiple primers in parallel. The arguments are the same as the main script.

# Extract methylation data

## Main script (4b)

### Principle

For each sequence in the alignment, count the number of each dinucleotide at the selected positions (from script 4a). Also count the number of haplotypes.
Counts are converted into proportion and written in csv files.

### Arguments

Named arguments:
- -i: A text file listing the paths of alignment files to process (usually all data)
- -o: Output directory
- -s: The path to the samples_infos.csv (generated at step 0)
- -c: The pickle file storing the CG selected positions
- -mode: The method for CG proportion computation [all or cgtg]

## slurm script

The slurm script is used to call the main script when using the sbatch script (see below). The slurm script will process only alignment files from one primer based on the task ID (job array).
It also merges the methylation data of the samples processed.

Caution: The mode is currently set to "all". Edit the slurm script to change the mode.
(would be good to improve this)

## sbatch script

The sbatch script allows to run the main script for multiple primers in parallel. The arguments are positional:
1. A text file listing the paths of alignment files to process (usually all data)
2. Output directory
3. The path to the samples_infos.csv (generated at step 0)

# Other scripts

- generate_alignment_plot.py
is used for further analyses if we want to generate an alignment plot for selected aligments files (not only healthy donors like in 4a).

# Example (user=kdasilva)

```
conda create -p /groups/proudhon_lab/projects/l1pa_meth/conf/extract_methyl -c conda-forge -c anaconda python matplotlib numpy pandas pickle
```

## Prepare list of aligns

```
ls /scratch/kdasilva/l1pameth_2023/alignments/default_scores/10_largest/*.align > /groups/proudhon_lab/projects/l1pa_meth/0_data/all_aligns_to_select_pos.txt 
```

## Position selection and alignment plots

All primers in parallel
```
cd ~/proudhon_lab/psl1_meth/

time scripts/4_extract_methylation_data/sbatch_select_positions.sh /groups/proudhon_lab/projects/l1pa_meth/0_data/all_aligns_to_select_pos.txt /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest /groups/proudhon_lab/projects/l1pa_meth/0_data/samples_infos.csv 20
```

## Run methylation data extraction

All samples in parallel
```
cd ~/proudhon_lab/psl1_meth/

time scripts/4_extract_methylation_data/sbatch_extract_methyl_haplo.sh /groups/proudhon_lab/projects/l1pa_meth/0_data/all_aligns_to_extract_methyl.txt /scratch/kdasilva/l1pameth_2023/methylation_data/test /groups/proudhon_lab/projects/l1pa_meth/0_data/samples_infos.csv
```

Time = 15m41.909s

## Merge infos

```
i=0
for s in $(ls /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest/vs_all/cg_methyl*); do
    if [[ $i -eq 0 ]] ; then
        head -1  "${s}" > /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest/cg_methyl.all_samples.csv
    fi
	tail -n +2  "${s}" >>  /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest/cg_methyl.all_samples.csv
	i=$(( $i + 1 ))
done

i=0
for s in $(ls /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest/vs_all/haplotypes*); do
    if [[ $i -eq 0 ]] ; then
        head -1  "${s}" > /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest/haplotypes.all_samples.csv
    fi
	tail -n +2  "${s}" >>  /scratch/kdasilva/l1pameth_2023/methylation_data/default_scores/10_largest/haplotypes.all_samples.csv
	i=$(( $i + 1 ))
done
```
