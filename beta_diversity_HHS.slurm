#!/bin/bash
#SBATCH --job-name=beta_diversity_HHS_fecal # Job name
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
INPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bracken2/HHS_fecal"
OUTPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/beta_diversity_plots"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Run beta diversity
echo "Running beta_diversity.py on all bracken files..."

python /common/lymphomics/coffeebean/murine_fecal_analysis/beta_diversity.py \
    -i ${INPUT_DIR}/*_kraken2_bracken.txt \
    --type bracken > ${OUTPUT_DIR}/beta_diversity_HHS.csv
echo "Beta diversity calculation complete!"
