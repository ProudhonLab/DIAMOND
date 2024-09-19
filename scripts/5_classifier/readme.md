# Classification

## Prepare data

Before running the classification scripts, prepare the data to process.
[TODO]

## Main script (5)

### Principle

Run random forests for healthy vs cancer classifications.

### Arguments

Named arguments:
- -d: csv file with methylation data for testing set
- -v: csv file with methylation data for 
- -r: Number of runs [5000]
- -b: balancing method for imbalanced classes [undersampling, oversampling, class_weight, None]
- -o: output directory
- -c: csv file with a color table (colors need to be in the second column)

## sbatch script

Contrary to the previous step, slurm and sbatch scripts are merged.

The sbatch script allows to run the main script for multiple sets of training/test data in parallel. The arguments are positional:
1. A text file with tab-separated paths to data used for classification. First column is the training set, second column is the testing set (if training only, put "None").
2. Output directory
3. Balance method [undersampling, oversampling, class_weigth, None]
4. csv file with a color table (colors need to be in the second column)

# Extract features

[TODO]


# Example (user=kdasilva)

Refaire les mêmes classifications du papier avec mes scripts.

```
mkdir -p /scratch/kdasilva/20230612_classic_classifications/discovery/inputs
mkdir -p /scratch/kdasilva/20230612_classic_classifications/discovery/results
mkdir -p /scratch/kdasilva/20230612_classic_classifications/validation/inputs
mkdir -p /scratch/kdasilva/20230612_classic_classifications/validation/results
```

## Prepare data
See prepare_data_for_classif.R

## List for discovery
```
cd /scratch/kdasilva/20230612_classic_classifications/
touch list_discovery.txt
for file in /scratch/kdasilva/20230612_classic_classifications/discovery/inputs/*; do echo -e "$file\tNone" >> list_discovery.txt; done
```

## List for validation
```
cd /scratch/kdasilva/20230612_classic_classifications/
touch list_validation.txt
for file in /scratch/kdasilva/20230612_classic_classifications/validation/inputs/*; do discovery=$(echo $file | sed 's/validation/discovery/g'); echo -e "$discovery\t$file" >> list_validation.txt; done
```

## Discovery
```
bash ~/proudhon_lab/psl1_meth/scripts/5_classifier/sbatch_classifications.sh /scratch/kdasilva/20230612_classic_classifications/list_discovery.txt /scratch/kdasilva/20230612_classic_classifications/discovery/results undersampling /groups/proudhon_lab/projects/l1pa_meth/5_classifications/color_table_c1.csv

bash ~/proudhon_lab/psl1_meth/scripts/5_classifier/sbatch_classifications.sh /scratch/kdasilva/20230612_classic_classifications/list_discovery.txt /scratch/kdasilva/20230612_classic_classifications/discovery/results None /groups/proudhon_lab/projects/l1pa_meth/5_classifications/color_table_c1.csv
```

## Validation
```
bash ~/proudhon_lab/psl1_meth/scripts/5_classifier/sbatch_classifications.sh /scratch/kdasilva/20230612_classic_classifications/list_validation.txt /scratch/kdasilva/20230612_classic_classifications/validation/results undersampling /groups/proudhon_lab/projects/l1pa_meth/5_classifications/color_table_c2.csv

bash ~/proudhon_lab/psl1_meth/scripts/5_classifier/sbatch_classifications.sh /scratch/kdasilva/20230612_classic_classifications/list_validation.txt /scratch/kdasilva/20230612_classic_classifications/validation/results None /groups/proudhon_lab/projects/l1pa_meth/5_classifications/color_table_c2.csv
```

## Extract scores
```
srun --mem=1G --pty bash
. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif2

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_scores.py -i  /scratch/kdasilva/20230612_classic_classifications/discovery/results/ -o /scratch/kdasilva/20230612_classic_classifications/discovery.scores_summary.csv

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_scores.py -i  /scratch/kdasilva/20230612_classic_classifications/validation/results/ -o /scratch/kdasilva/20230612_classic_classifications/validation.scores_summary.csv
```

## Extract predictions
```
srun --mem=1G --pty bash
. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif2

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_predictions.py -i  /scratch/kdasilva/20230612_classic_classifications/discovery/results/ -o /scratch/kdasilva/20230612_classic_classifications/discovery.predictions_summary.csv
```

## Extract performances
```
srun --mem=1G --pty bash
. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif2

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_sensib_at_99perc_spec.py -i  /scratch/kdasilva/20230612_classic_classifications/discovery/results/ -o /scratch/kdasilva/20230612_classic_classifications/discovery.sensibilities_summary.csv

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_sensib_at_99perc_spec.py -i  /scratch/kdasilva/20230612_classic_classifications/validation/results/ -o /scratch/kdasilva/20230612_classic_classifications/validation.sensibilities_summary.csv
```

## Extract features
```
srun --mem=1G --pty bash
. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif2

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_features_importance.py -i  /scratch/kdasilva/20230612_classic_classifications/discovery/results/ -o /scratch/kdasilva/20230612_classic_classifications/discovery.features_importance_summary.csv

python3 ~/proudhon_lab/psl1_meth/scripts/5_classifier/extract_features_importance.py -i  /scratch/kdasilva/20230612_classic_classifications/validation/results/ -o /scratch/kdasilva/20230612_classic_classifications/validation.features_importance_summary.csv
```

## Stack and blind models
Stacked and blind models ( a cancer subtype or type is removed each run to the testing set) are in the stack_blind brach on this repo.
