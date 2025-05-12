#!/bin/bash
#SBATCH --job-name=bowtie2_HHS # Job name
#SBATCH --output=outout_%j.log # Output file (%j will be replaced by the job ID)
#SBATCH --error=error_%j.log # Error log file (%j wil be replaced by the job ID)
#SBATCH --time=96:00:00 # Wall time (total time to run the job): 8 hours
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks=1 # Number of tasks
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --mem=200GB
#SBATCH --partition=batch,guest
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=snakamura@unomaha.edu
#SBATCH --licenses=common

# Load necessary modules
module load bowtie/2.5

# Directories
INPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/cutadapt/HHS_fecal"
OUTPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bowtie2/HHS_fecal"
HUMAN_DIR="/common/claytonlab/coffeebean/GRCh38_noalt_decoy_as/"

# Run the job
# Loop through all the trimmed fastq files in the directory
echo "Starting murine fecal Bowtie2 (host contamination removal)"
for file in ${INPUT_DIR}/*.fastq; do
    echo "Processing ${file}"

    # Get the sample prefix (e.g., 'SRS260331') by stripping off the sample-specific parts
    sample_prefix=$(basename ${file} | cut -d'.' -f1)

    # Define the corresponding paired-end files for L001 and L002
    r1="${INPUT_DIR}/${sample_prefix}.denovo_duplicates_marked.trimmed.1_trimmed.fastq"
    r2="${INPUT_DIR}/${sample_prefix}.denovo_duplicates_marked.trimmed.2_trimmed.fastq"

    # Output files for aligned reads
    output_l_1="${OUTPUT_DIR}/${sample_prefix}_1_unmapped.fastq"
    output_l_2="${OUTPUT_DIR}/${sample_prefix}_2_unmapped.fastq"

    # Run Bowtie2 for paired-end alignment
    bowtie2 -p 8 -x ${HUMAN_DIR}GRCh38_noalt_decoy_as -1 ${r1} -2 ${r2} --un-conc-gz ${OUTPUT_DIR}/${sample_prefix}_unmapped.fastq.gz -S ${OUTPUT_DIR}/${sample_prefix}.sam

    echo "Processed ${sample_prefix} for both L001 and L002"
done

echo "Done with murine fecal Bowtie2 (host contamination removal)"
