#!/bin/bash
#SBATCH --job-name=bracken2_standard_EA_class_level
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

# Load Bracken module (adjust as needed)
module load bracken/2.6
module load kraken2/2.1

# Set paths
KRAKEN_DB="/work/lymphomics/coffeebean/kraken2_db"
READ_LENGTH=100  # Adjust based on sequencing reads
INPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/kraken2_results/16S_GM_EA"
OUTPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bracken2_class/16S_GM_EA"
mkdir -p ${OUTPUT_DIR}

# Set parameters
THRESHOLD=10

# Loop through all .kreport files
for kreport_file in ${INPUT_DIR}/*_report.txt; do
    # Get sample prefix
    sample_prefix=$(basename ${kreport_file} _report.txt)

    # Run Bracken at species level
    bracken -d ${KRAKEN_DB} \
            -i ${kreport_file} \
            -o ${OUTPUT_DIR}/${sample_prefix}_bracken.txt \
            -t ${THRESHOLD} \
            -r ${READ_LENGTH} \
            -l C


    echo "Processed ${kreport_file}"
done

echo "Bracken processing for taxonomic classification completed."
