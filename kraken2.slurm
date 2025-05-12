#!/bin/bash
#SBATCH --job-name=kraken2_taxonomic_classification # Job name
#SBATCH --output=outout_%j.log # Output file (%j will be replaced by the job ID)
#SBATCH --error=error_%j.log # Error log file (%j wil be replaced by the job ID)
#SBATCH --time=96:00:00 # Wall time (total time to run the job): 8 hours
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks=1 # Number of tasks
#SBATCH --cpus-per-task=8  # Adjust based on your job
#SBATCH --mem-per-cpu=25G  # This ensures 8 CPUs get a total of 200GB
#SBATCH --partition=batch,guest
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=snakamura@unomaha.edu
#SBATCH --licenses=common

# Load necessary modules
echo "Loading Kraken2..."
module load kraken2/2.1

# Define input and output directories
INPUT_DIRS=(
    "/common/biocore/coffeebean/senior_project_analysis/bowtie2/16S_GM_AO"
    "/common/biocore/coffeebean/senior_project_analysis/bowtie2/16S_GM_EA"
    "/common/biocore/coffeebean/senior_project_analysis/bowtie2/HHS_fecal"
)
OUTPUT_PARENT_DIR="/common/biocore/coffeebean/senior_project_analysis/kraken2_results"
OUTPUT_DIRS=(
    "${OUTPUT_PARENT_DIR}/16S_GM_AO"
    "${OUTPUT_PARENT_DIR}/16S_GM_EA"
    "${OUTPUT_PARENT_DIR}/HHS_fecal"
)

# Ensure output directories exist
for DIR in "${OUTPUT_DIRS[@]}"; do
    mkdir -p "$DIR"
done

# Function to process samples
process_sample() {
    INPUT_DIR=$1
    OUTPUT_DIR=$2
    SAMPLE=$3

    UNMAPPED_R1="${INPUT_DIR}/${SAMPLE}_unmapped.fastq.1.gz"
    UNMAPPED_R2="${INPUT_DIR}/${SAMPLE}_unmapped.fastq.2.gz"
    OUTPUT_REPORT="${OUTPUT_DIR}/${SAMPLE}_kraken2_report.txt"
    OUTPUT_CLASSIFIED="${OUTPUT_DIR}/${SAMPLE}_kraken2_classified.txt"

    echo "Running Kraken2 for $SAMPLE..."
    kraken2 --db "$KRAKEN2_DB" \
            --threads 8 \
            --paired "$UNMAPPED_R1" "$UNMAPPED_R2" \
            --report "$OUTPUT_REPORT" \
            --output "$OUTPUT_CLASSIFIED"

    echo "Completed Kraken2 for $SAMPLE!"
}

# Loop through input directories and process all unmapped read files sequentially
for i in "${!INPUT_DIRS[@]}"; do
    INPUT_DIR="${INPUT_DIRS[$i]}"
    OUTPUT_DIR="${OUTPUT_DIRS[$i]}"
    for FILE in "$INPUT_DIR"/*_unmapped.fastq.1.gz; do
        SAMPLE=$(basename "$FILE" _unmapped.fastq.1.gz)  # Extract sample name
        process_sample "$INPUT_DIR" "$OUTPUT_DIR" "$SAMPLE"
    done
done

echo "All samples processed!"
