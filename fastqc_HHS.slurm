#!/bin/bash
#SBATCH --job-name=fastqc_HHS       # Job name
#SBATCH --output=output_%j.log  # Output log file (%j will be replaced by the job ID)
#SBATCH --error=error_%j.log    # Error log file (%j will be replaced by the job ID)
#SBATCH --time=24:00:00         # Wall time (total time to run the job): 24 hours
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=1              # Number of tasks
#SBATCH --cpus-per-task=4       # Number of CPU cores per task
#SBATCH --mem=200GB
#SBATCH --partition=batch,guest
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=snakamura@unomaha.edu
#SBATCH --licenses=common

# Load necessary modules
module load fastqc/0.12

# Run the job


echo "Done with HHS_fecal..."
fastqc /common/biocore/coffeebean/senior_project_data/HHS_fecal/*.fastq -o /common/biocore/coffeebean/senior_project_analysis/fastqc/HHS_fecal
echo "Done with HHS_fecal..." 
~                                                 
