#!/bin/bash

# Arguments

script_dir=$( dirname "$0" )

input_file=$1
num_samples=$(wc -l < "$input_file")

output_dir=$2
mkdir -p "${output_dir}/logs"

balance_method=$3 # [undersampling, oversampling, class_weigth, None]

color_table_path=$4

sbatch << EOF
#!/bin/bash
#SBATCH --job-name=classifications
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --output="${output_dir}/logs/classif-%A_%a.out"
#SBATCH --array=1-${num_samples}

# Environment
. /local/env/envconda.sh
conda activate /groups/proudhon_lab/projects/l1pa_meth/conf/classif2

# Extract the corresponding line from the input file for this job
line_id="\${SLURM_ARRAY_TASK_ID}p"
line=\$(sed -n "\${line_id}" "$input_file")

# Split the line into two variables based on the tab separator
IFS=$'\t' read -r training_set testing_set <<< "\$line"

# Run the script

if [ "\$testing_set" == "None" ]; then
    output_dir="${output_dir}/training.\$(basename \${training_set%.*}).${balance_method}"
    python3 "$script_dir"/5_random_forests_classification.py -d "\${training_set}" -b "$balance_method" -o "\${output_dir}" -c "$color_table_path"
else
    output_dir="${output_dir}/testing.\$(basename \${testing_set%.*}).${balance_method}"
    python3 "$script_dir"/5_random_forests_classification.py -d "\${training_set}" -v "\${testing_set}" -b "$balance_method" -o "\${output_dir}" -c "$color_table_path"
fi

echo "THE END"
EOF