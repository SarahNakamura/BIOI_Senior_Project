#!/bin/bash
#SBATCH --job-name=cutadapt_HHS # Job name
#SBATCH --output=output_%j.log # Output log filr (%j will be replaced by the job ID)
#SBATCH --error=error_%j.log # Error log file (%j will be replaced by the job ID) 
#SBATCH --time=24:00:00 # Wall time (total time to eint the job): 8 hours
#SBATCH --nodes=1 # Number of nodes
#SBATCH --ntasks=1 # Number of tasks
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH --mem=200GB
#SBATCH --partition=batch,guest
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=snakamura@unomaha.edu


# Load necessary modules
module load cutadapt/2.9

# Directories
INPUT_DIR="/common/biocore/coffeebean/senior_project_data/HHS_fecal"
OUTPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/cutadapt/HHS_fecal"

# Run the job
# Loop through all the fastq files in the directory
echo "Starting murine fecal cutadapt..."

for file in ${INPUT_DIR}/*fastq; do
        echo ${file}
        # Extract the base name of the file (without path and extension)
        base_name=$(basename ${file} .fastq)

	# Run cutadapt over file
	cutadapt -u -15 -o ${OUTPUT_DIR}/${base_name}_trimmed.fastq  ${file}
        
        echo "Processed: ${file}"
done

echo "Done with murine fecal cutadapt..."
