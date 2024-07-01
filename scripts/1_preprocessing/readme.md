# Preprocessing scripts 

## Main script (1a)

### Principle


The final output is fasta files of good quality reads, demultiplexed, and deduplicated.

Reads are first filtered using fastp.
Then, for each primer, reads are trimmed of their primers and UMIs.
Remaining inserts are filtered by size and deduplicated based on the UMIs.

Various statistics are stored in the process.


### Arguments

Positional arguments:
1. The path of the folder containing the fastq file(s). The script automatically detects the R1 and/or R2 files, but to work it needs to have R1/R2 in its name. The sample name will be the string before "R1".
2. The path to the output folder. Inside it, a new folder named after the sample will be created.
3. The strand to process [R1, R2, or merge]
4. The path to a csv file containing the infos of the primers. The csv should have 6 columns: primer name, seq forward, seq reverse, seq reverse complement, ref size, with CG boolean
5. The size of the UMIs [16 for DIAMOND]
6. The number of threads

## slurm script

The slurm script is only used to call the main script when using the sbatch script (see below). The input is a folder containing multiple sequencing data folders or a text file listing the samples paths. The slurm script will process only one sample based on the task ID (job array).

## sbatch script

The sbatch script allows to run the main script for multiple samples in parallel. The arguments are the same as the preprocessing script except that the input directory is a folder containing multiple sequencing data folders or a text file listing the samples paths.

# Merge script (1b)

## Principle

The outputs for all processed samples can be merged in single files with this script.

## Arguments

Positional arguments:
1. The path of the output folder
2. The path to a csv file containing the infos of the primers. The csv should have 6 columns: primer name, seq forward, seq reverse, seq reverse complement, ref size, with CG boolean
3. A prefix used for the output files generated.

# Other scripts

- adapter_stats.R
- concatenate_UMIs_targets.py
- info_file_to_fasta.py
- plot_read_or_uid_length.R
- trim_UMIs.py
are scripts used inside the 1a preprocessing script.

Written by Marc Michel and untouched since.


# Example (user=kdasilva)

## Create env

```
ssh kdasilva@genossh.genouest.org
tmux
srun --mem 20G --pty bash
. /local/env/envconda.sh
conda create -p /groups/proudhon_lab/projects/l1pa_meth/conf/preproc -c bioconda -c conda-forge -c anaconda boost python biopython r-base rtools r-tidyverse r-optparse r-plotly r-plyr mothur vsearch fastp atropos seqkit 
```

## Prepare list of samples to preprocess

MHOVC 1, 2 et 15 ont été demandé ensuite pour les rajouter aux analyses déjà faites auparavant.
Attention, MHOVC 1 et 2 n'ont pas d'UMIs, ils n'utilisaient pas non plus les mêmes méthodes, d'où l'option merge plutôt que de récupérer uniquement R2.

```
ls -d /groups/proudhon_lab/sequencing_data/2000014/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc1_mergeMode.txt
ls -d /groups/proudhon_lab/sequencing_data/2001929/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc2_mergeMode.txt
ls -d /groups/proudhon_lab/sequencing_data/2010164/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc15_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2010131/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc13_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2009579/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc10_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2008909/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc8_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2008910/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc9_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2012336/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc16_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2010130/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc14_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2010012/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc12_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2012811/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc17_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2015506/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc20_r2Mode.txt
ls -d /groups/proudhon_lab/sequencing_data/2016240/*/ | sed 's|[/]$||g' > /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_mhovc22_r2Mode.txt

cat /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_*_r2Mode.txt > /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastq_to_preproc_r2Mode.txt
cat /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_*_mergeMode.txt > /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastq_to_preproc_mergeMode.txt

rm /groups/proudhon_lab/projects/l1pa_meth/0_data/fastq_to_preproc_*
```

## Generate samples infos

Script R generate_mega_sample_infos.R pour générer une matrice contenant tous les échantillons préprocessés, et récupérer leurs informations cliniques et autres via les fichiers du One Drive.

## Run preprocessing

```
srun -p tiny --pty bash
cd ~/proudhon_lab/psl1_meth/
mkdir /scratch/kdasilva/l1pameth_2023/

time scripts/1_preprocessing/sbatch_dataset_preprocessing.sh /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastq_to_preproc_mergeMode.txt /scratch/kdasilva/l1pameth_2023/preprocessing_merge merge /groups/proudhon_lab/projects/l1pa_meth/0_data/primers_MHOVC_dataset_2010012.csv 16 12

time scripts/1_preprocessing/sbatch_dataset_preprocessing.sh /groups/proudhon_lab/projects/l1pa_meth/0_data/all_fastq_to_preproc_r2Mode.txt /scratch/kdasilva/l1pameth_2023/preprocessing_r2 R2 /groups/proudhon_lab/projects/l1pa_meth/0_data/primers_MHOVC_dataset_2010012.csv 16 12
```

1074 samples processed with R2 mode.
Time r2Mode = 1052m1.295s (17h30min)

62 samples processed with merge mode.
Time mergeMode = 12m58.438s

Total: 1136 samples

## Check for empty inserts

```
find /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2/ -type f -empty -name "inserts_dereplicated.*.fasta" > /scratch/kdasilva/l1pameth_2023/preproc_failed_summary.txt
find /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_merge/ -type f -empty -name "inserts.*.fasta" >> /scratch/kdasilva/l1pameth_2023/preproc_failed_summary.txt
```

No inserts L1HS15 for 10 samples :
* L53S1
* L53S4
* L53S5
* L56S6
* L53S9
* L53S10
* L53S11
* L53S12
* L53S13
* L53S15

## Merge preprocessing infos

Can be done with a script for all samples preprocessed using R2 mode.
```
time scripts/1_preprocessing/1b_merge_infos.sh /groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2 /groups/proudhon_lab/projects/l1pa_meth/0_data/primers_MHOVC_dataset_2010012.csv allSamplesR2
```

Has to be done manually for MHOVC 1 and 2 because there is no UMIs, so the files to use are inserts.length_filtered rather than inserts_dereplicated.
```
outputDir="/groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_merge"
primersInfo="/groups/proudhon_lab/projects/l1pa_meth/0_data/primers_MHOVC_dataset_2010012.csv"
samples=($(ls -d ${outputDir}/*/ | grep -v "logs/" | sed 's|[/]$||g'))
datasetName="allSamplesMerge"
i=0
for s in ${samples[@]}/; do
    sample=$(basename ${s})
    if [[ $i -eq 0 ]] ; then
        head -1  "${outputDir}/${sample}/${sample}.summary.csv" > "${outputDir}/${datasetName}.summary.csv"
        head -1  "${outputDir}/${sample}/${sample}.summary.reads.csv" > "${outputDir}/${datasetName}.summary.reads.csv"
    fi
	tail -n +2  "${outputDir}/${sample}/${sample}.summary.csv" >>  "${outputDir}/${datasetName}.summary.csv"
    tail -n +2  "${outputDir}/${sample}/${sample}.summary.reads.csv" >>  "${outputDir}/${datasetName}.summary.reads.csv"
	i=$(( $i + 1 ))
done
```