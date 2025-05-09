#!/bin/bash
#SBATCH --job-name=alpha_diversity_class_EA
#SBATCH --output=output_%j.out
#SBATCH --error=error_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --partition=batch,guest
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=snakamura@unomaha.edu
#SBATCH --licenses=common

# Load required modules
module load krakentools/py310/2.1
module load python/3.9

# Set paths
INPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bracken2_class/16S_GM_EA"
OUTPUT_FILE="/common/biocore/coffeebean/senior_project_analysis/alpha_diversity_class/alpha_diversity_EA.csv"

# Ensure the CSV file has a header
echo "Sample,Shannon_Diversity" > $OUTPUT_FILE

# Loop through all txt files
for bracken2_file in ${INPUT_DIR}/*_kraken2_bracken.txt; do
    # Get sample prefix (remove "_bracken.txt" to get sample name)
    sample_prefix=$(basename "$bracken2_file" | sed 's/_kraken2_bracken.txt//')

    # Run alpha_diversity.py and capture the Shannon diversity index
    shannon_value=$(python /common/lymphomics/coffeebean/murine_fecal_analysis/alpha_diversity.py -f ${bracken2_file} | awk '{print $NF}')

    # Append results to CSV
    echo "$sample_prefix,$shannon_value" >> $OUTPUT_FILE

    echo "Processed ${bracken2_file}"
done


echo "Alpha diversity calculation complete!"
