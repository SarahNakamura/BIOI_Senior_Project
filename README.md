# BIOI_Senior_Project
This is the code repository for my bioinformatics senior project. Most steps within the pipeline were completed through job submissions using slurm files. The slurm files can be run via command line as follows:
```
sbatch path/to/slurm_file.slurm
```
## Data Download
The data for this project were produced from the NIH Human Microbiome Project and was retrieved from their designated database. With advanced search, I have narrowed down the data of interest in the FASTQ format. The overview of the data retrieved is below:
| Study | Number of files | Type | File format|
|:-----:|:---------------:|:----:|:----------:|
|Obesity Study (16S-GM-AO)|12|16S sequence|FASTQ|
|Esophageal Adenocarcinoma Study (16S-GM-EA)|28|16S sequence|FASTQ|
|Human Microbiome Project (control/HHS)|24|16S sequence|FASTQ|
|(total)|64|||
- The summary of the 16S shotgun sequencing data files retrieved from NIH Human Mocrobime Project 
Links for the advance search in the Human Microbiome Project Database is the following:
- [16S-GM-AO & 16S-GM-EA](https://portal.hmpdacc.org/search/s?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.body_site%22,%22value%22:%5B%22gastrointestinal%20tract%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22subject.project_name%22,%22value%22:%5B%22Human%20Microbiome%20Project%20(HMP)%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.study_name%22,%22value%22:%5B%2216S-GM-AO%22,%2216S-GM-EA%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22file.format%22,%22value%22:%5B%22FASTQ%22%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:100,%22sort%22:%22file_id:asc%22%7D%7D&facetTab=cases)
- [Control](https://portal.hmpdacc.org/search/s?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.visit_visit_number%22,%22value%22:%5B1%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22subject.project_name%22,%22value%22:%5B%22Human%20Microbiome%20Project%20(HMP)%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.study_name%22,%22value%22:%5B%22HHS%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22sample.body_site%22,%22value%22:%5B%22feces%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22file.format%22,%22value%22:%5B%22FASTQ%22%5D%7D%7D%5D%7D&facetTab=cases)

The files were downloaded to Holland Computing Center (HCC) under my account using the manifest file generated from the database and Portal-client on HCC. Portal-client is a python-based client used for downloading large datafiles. It will read the manifest files as an input and would download the files to the specified destination. HCC has a pre-installed version of Portal-client (igs-portal-client/1.4) and will be available once installed. (If a newer version if available, I suggest to use the updated version of the software.)
Reference: Igs. (n.d.). IGS/portal_client: Python-based client for downloading data made available through portals powered by the GDC-based portal system.. GitHub. https://github.com/IGS/portal_client 
```
# Load portal-client to your system
module load igs-portal-client/1.4
# Run portal-client and download all files
portal_client -m path/to/MANIFEST_FILE -d path/to/destination
```

## Quality Control
### 1.FastQC
FastQC was used to analyze the quality of the raw sequence file. Reference:Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Bash scripts (Since the sequence files were very large, I wrote three separate bash scripts to run them separately):
- [FastQC_AO](fastqc_AO.slurm)
- [FastQC_EA](fastqc_EA.slurm)
- [FastQC_HHS](fastqc_HHS.slurm)

The link will take you to the designated location where the file resides, but below is FastQC slurm file written for AO.
```
#!/bin/bash
#SBATCH --job-name=fastqc_AO       # Job name
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
echo "Starting 16S_GM_AO..."
fastqc /common/biocore/coffeebean/senior_project_data/16S_GM_AO/*.fastq -o /common/biocore/coffeebean/senior_project_analysis/fastqc/16S_GM_AO
echo "Done with 16S_GM_AO..."

~                                                 
```
### 2.MultiQC
After running FastQC, we could run MiltiQC to visualize the sequence quality. Reference:P. Ewels, M. Magnusson, S. Lundin, and M. Käller. 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 32(19): 3047-3048.

For this, I ran
```
module load multiqc/py37/1.8
multiqc -o path/to/output_dir path/to/input_dir
```
## Adapter Removal
From the MultiQC report, there seem to be conatmination around the last 15 bases in the control dataset. Therefore, I ran cutadapt to remove the low sequnece quality region from all seuquence data in the control group. Cutadapt documentation can be found [here](https://cutadapt.readthedocs.io/en/stable/).
Reference: Martin, M. (2011). CUTADAPT removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1), 10. https://doi.org/10.14806/ej.17.1.200 

Slurm script can be found here.
```
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
```
## Contamination Removal
For the host contamination removal step, "GRCh38_noalt_decoy_as" was retrieved from the Bowtie2 documentation. The documentation can be found [here](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).
Reference:  Wood DE, Lu J, Langmead B. Improved metagenomic analysis with Kraken 2 (2019). Genome Biology. 2019 Nov;p. 76230

The slurm files can be found below:
- [Bowtie2_AO](bowtie2_AO.slurm)
- [Bowtie2_EA](bowtie2_EA.slurm)
- [Bowtie2_HHS](bowtie2_HHS.slurm)

Below is the slurm file for AO datasets.
```
#!/bin/bash
#SBATCH --job-name=bowtie2_AO # Job name
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
INPUT_DIR="/common/biocore/coffeebean/senior_project_data/16S_GM_AO"
OUTPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bowtie2/16S_GM_AO_test"
HUMAN_DIR="/common/claytonlab/coffeebean/GRCh38_noalt_decoy_as/"

# Run the job
# Loop through all the trimmed fastq files in the directory
echo "Starting murine fecal Bowtie2 (host contamination removal)"
for file in ${INPUT_DIR}/*.fastq; do
    echo "Processing ${file}"

    # Get the sample prefix (e.g., 'SRS260331') by stripping off the sample-specific parts
    sample_prefix=$(basename ${file} | cut -d'.' -f1)

    # Define the corresponding paired-end files for L001 and L002
    r1="${INPUT_DIR}/${sample_prefix}.denovo_duplicates_marked.trimmed.1.fastq"
    r2="${INPUT_DIR}/${sample_prefix}.denovo_duplicates_marked.trimmed.2.fastq"

    # Output files for aligned reads
    output_l_1="${OUTPUT_DIR}/${sample_prefix}_1_unmapped.fastq"
    output_l_2="${OUTPUT_DIR}/${sample_prefix}_2_unmapped.fastq"

    # Run Bowtie2 for paired-end alignment
    bowtie2 -p 8 -x ${HUMAN_DIR}GRCh38_noalt_decoy_as -1 ${r1} -2 ${r2} --un-conc-gz ${OUTPUT_DIR}/${sample_prefix}_unmapped.fastq.gz -S ${OUTPUT_DIR}/${sample_prefix}.sam

    echo "Processed ${sample_prefix} for both L001 and L002"
done

echo "Done with murine fecal Bowtie2 (host contamination removal)"
```
## Taxonomic Classification
Kraken2 was used for the taxonomic abundance analysis. The slurm script can be found [here](kraken2.slurm) and below. The pre-intalled Kraken2 standard database in HCC was used in this step. Reference: Lu, J., Rincon, N., Wood, D.E. et al. Metagenome analysis using the Kraken software suite. Nat Protoc 17, 2815–2839 (2022). https://doi.org/10.1038/s41596-022-00738-y![image](https://github.com/user-attachments/assets/a0518553-2f7f-42d0-a095-7f643f620b86)

Kraken2 documentation can be found [here](https://ccb.jhu.edu/software/kraken2/).

The following is the slurm script for Kraken2.
```
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
```
## Taxonomic Abunance Analysis
The abundance analysis was performed using bracken in the class level. Reference: Lu J, Breitwieser FP, Thielen P, Salzberg SL. 2017. Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science 3:e104 https://doi.org/10.7717/peerj-cs.104

The slurm scripts are as follows:
- [Bracken_AO](bracken2_bacteria_AO.slurm)
- [Bracken_EA](bracken2_bacteria_EA.slurm)
- [Bracke_HHS](bracken2_bacteria_HHS.slurm)

Below is the slurm file for AO dataset.
```
#!/bin/bash
#SBATCH --job-name=bracken2_standard_AO_class_level
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
INPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/kraken2_results/16S_GM_AO"
OUTPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bracken2_class/16S_GM_AO"
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
```
The taxonomic classification data are separated into files by sample. To combine the abundance data from all samples into one for each project group, KBKen was used. The "combine_bracken_outputs.py" was used. Reference: SumeetTiwari07. (n.d.). Sumeettiwari07/KBKen: Converting Kraken2 and Bracken output into phyloseq object. GitHub. https://github.com/SumeetTiwari07/KBKen

To run this:
```
python3 combine_bracken_outputs.py --files path/to/*bracken_outputfiles -o path/to/combined_bracken_output
```

## Abundance analysis Top10 abundant taxa
In order to identify the top10 abundant taxa across three sample groups, the following R script was written:
```
# taxon_boxplot_with_stats.R

# Load libraries
library(tidyverse)
library(ggpubr)  # For stat_compare_means

# ----------------- FUNCTION TO LOAD BRACKEN FRACTION FILES -----------------
load_bracken_frac_long <- function(file_path, project_label) {
  df <- read.delim(file_path, check.names = FALSE)
  frac_cols <- grep("_frac$", colnames(df), value = TRUE)
  
  abundance <- df[, c("name", frac_cols)]
  abundance_long <- pivot_longer(abundance, cols = -name,
                                 names_to = "SampleID",
                                 values_to = "Abundance")
  
  abundance_long$Project <- project_label
  return(abundance_long)
}

# ----------------- LOAD ALL THREE FILES -----------------
ao <- load_bracken_frac_long("combined_abundance_class_AO.txt", "AO")
ea <- load_bracken_frac_long("combined_abundance_class_EA.txt", "EA")
hhs <- load_bracken_frac_long("combined_abundance_class_HHS.txt", "HHS")

# Combine all
combined <- bind_rows(ao, ea, hhs)

# ----------------- FILTER TO TOP 10 MOST ABUNDANT TAXA -----------------
top_taxa <- combined %>%
  group_by(name) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>%
  pull(name)

filtered <- combined %>%
  filter(name %in% top_taxa)

# ----------------- RUN KRUSKAL-WALLIS TEST FOR EACH TAXON -----------------
kw_results <- filtered %>%
  group_by(name) %>%
  summarise(p_value = kruskal.test(Abundance ~ Project)$p.value)

# Adjust p-values (FDR correction)
kw_results$p_adj <- p.adjust(kw_results$p_value, method = "fdr")

# Merge adjusted p-values back into data
filtered_stats <- left_join(filtered, kw_results, by = "name")

# ----------------- PLOT WITH SIGNIFICANCE -----------------
# Function to generate significance label from p-value
get_sig_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

kw_results$sig_label <- sapply(kw_results$p_adj, get_sig_label)

# Base plot
p <- ggplot(filtered, aes(x = Project, y = Abundance, fill = Project)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ name, scales = "free_y") +
  theme_minimal() +
  labs(title = "Top 10 Taxa Abundance Across Projects (with Kruskal–Wallis p-values)",
       y = "Relative Abundance",
       x = "Project") +
  theme(legend.position = "none")

# Add p-values to each facet
p + stat_compare_means(method = "kruskal.test", label = "p.format")
# (Displays raw p-values on plot)

# ----------------- SAVE PLOT -----------------
ggsave("taxon_boxplot_with_stats.pdf", width = 14, height = 8, dpi = 300)

# ----------------- OPTIONAL: Export stats table -----------------
write.csv(kw_results, "kruskal_results_top10_taxa.csv", row.names = FALSE)
```

## Diversity Analysis (apha diversity & beta diversity)
Python script from KrakenTools was used in alpha diversity. The GitHub repository for KrakenTools can be found [here](https://github.com/jenniferlu717/KrakenTools)

Reference: [Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)] (https://www.nature.com/articles/s41596-022-00738-y)

The slurm scripts for alpha diversity analysis is as follows:
- [Alpha_diversity_AO](alpha_diversity_AO.slurm)
- [Alpha_diversity_EA](alpha_diversity_EA.slurm)
- [Alpha_diversity_HHS](alpha_diversity_HHS.slurm)

Below is the slurm file for the AO dataset.
```
#!/bin/bash
#SBATCH --job-name=alpha_diversity_class_AO
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
INPUT_DIR="/common/biocore/coffeebean/senior_project_analysis/bracken2_class/16S_GM_AO"
OUTPUT_FILE="/common/biocore/coffeebean/senior_project_analysis/alpha_diversity_class/alpha_diversity_AO.csv"

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
```

It takes the abundance results from bracken as input and outputs a csv file that contain the alpha diversity metrics.

For visualization, the following R script was ran:

```
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)

# Load the CSV files and add a 'Group' column to each
data1 <- read.csv("alpha_diversity_class/alpha_diversity_AO.csv", header = TRUE)
data1$Group <- "16S_GM_AO"

data2 <- read.csv("alpha_diversity_class/alpha_diversity_EA.csv", header = TRUE)
data2$Group <- "16S_GM_EA"

data3 <- read.csv("alpha_diversity_class/alpha_diversity_HHS.csv", header = TRUE)
data3$Group <- "HHS_fecal"

# Combine all datasets into one dataframe
data <- bind_rows(data1,data2,data3)

# Check structure of merged data
head(data)

# Use ANOVA for statistical test since there are three data groups
anova_result <- aov(Shannon_Diversity ~ Group, data = data)
summary(anova_result)
TukeyHSD(anova_result)


ggplot(data, aes(x = Group, y = Shannon_Diversity, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +  # Add individual points
  stat_compare_means(method = "anova") +  # ANOVA p-value (change to "kruskal.test" if non-parametric)
  theme_minimal() +
  labs(title = "Shannon Diversity Across Datasets",
       x = "Dataset Group",
       y = "Shannon Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("shannon_diversity_boxplot_class.pdf", width = 8, height = 6)
```

The R for beta diversity analysis is as follows using PCoA:
```
# multi_project_pcoa.R

# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)

# ---------- FUNCTION TO LOAD AND MERGE BRACKEN FILES ----------
load_bracken_frac_files <- function(file_paths, project_labels) {
  if (length(file_paths) != length(project_labels)) {
    stop("file_paths and project_labels must have the same length.")
  }
  
  all_data <- list()
  
  for (i in seq_along(file_paths)) {
    file <- file_paths[i]
    label <- project_labels[i]
    
    df <- read.delim(file, check.names = FALSE)
    frac_cols <- grep("_frac$", colnames(df), value = TRUE)
    
    # Subset and transpose
    abundance <- df[, c("name", frac_cols)]
    rownames(abundance) <- abundance$name
    abundance$name <- NULL
    transposed <- as.data.frame(t(abundance))
    
    # Add metadata
    transposed$SampleID <- rownames(transposed)
    transposed$Project <- label
    
    all_data[[i]] <- transposed
  }
  
  # Combine all samples
  combined <- bind_rows(all_data)
  rownames(combined) <- combined$SampleID
  metadata <- combined[, c("SampleID", "Project")]
  abundance_matrix <- combined[, !(colnames(combined) %in% c("SampleID", "Project"))]
  
  # Fill missing taxa with 0s
  abundance_matrix[is.na(abundance_matrix)] <- 0
  
  # Convert all to numeric
  abundance_matrix <- as.data.frame(lapply(abundance_matrix, as.numeric))
  rownames(abundance_matrix) <- metadata$SampleID
  
  return(list(
    abundance_matrix = abundance_matrix,
    metadata = metadata
  ))
}

# ---------- FILE PATHS AND LABELS ----------
file_paths <- c(
  "combined_abundance_class_AO.txt",
  "combined_abundance_class_EA.txt",
  "combined_abundance_class_HHS.txt"
)
project_labels <- c("AO", "EA", "HHS")

# ---------- LOAD AND MERGE DATA ----------
result <- load_bracken_frac_files(file_paths, project_labels)
abundance_matrix <- result$abundance_matrix
metadata <- result$metadata

# ---------- FILTER LOW-PREVALENCE TAXA ----------
taxa_prevalence <- colSums(abundance_matrix > 0, na.rm = TRUE)
if (sum(taxa_prevalence >= 3, na.rm = TRUE) > 0) {
  abundance_matrix <- abundance_matrix[, taxa_prevalence >= 3]
} else {
  warning("No taxa present in 3 or more samples. Skipping filter.")
}

# ---------- PCoA ANALYSIS ----------
dist_mat <- vegdist(abundance_matrix, method = "bray")
pcoa_res <- cmdscale(dist_mat, eig = TRUE, k = 2)

# ---------- BUILD PLOT DATA ----------
pcoa_df <- as.data.frame(pcoa_res$points)
pcoa_df$SampleID <- rownames(pcoa_df)
plot_df <- left_join(pcoa_df, metadata, by = "SampleID")

# ---------- PLOT ----------
ggplot(plot_df, aes(x = V1, y = V2, color = Project)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCoA of Microbiome Composition Across Projects",
    x = "PC1", y = "PC2", color = "Project"
  )

# ---------- OPTIONAL: SAVE PLOT ----------
ggsave("pcoa_across_projects_no_label.pdf", width = 8, height = 6, dpi = 300)
```

## Network Analysis
Taxonomic abundance networks were constructed for all three sample groups. The R script is as folows (change input combined abundance file to the project of interest):
```
# taxon_network.R

# Load libraries
library(tidyverse)
library(Hmisc)
library(igraph)
library(reshape2)
library(ggraph)
library(tidygraph)

# Read and parse input file
df <- read.delim("combined_abundance_class_HHS.txt", check.names = FALSE)

# Extract only the "_frac" columns (relative abundance)
frac_cols <- grep("_frac$", colnames(df), value = TRUE)
abundance_frac <- df[, c("name", frac_cols)]
rownames(abundance_frac) <- abundance_frac$name
abundance_frac$name <- NULL

# Transpose: rows = samples, columns = taxa
data_t <- as.data.frame(t(abundance_frac))

# Optional: filter rare taxa
keep_taxa <- colSums(data_t > 0) >= nrow(data_t) * 0.2
data_filtered <- data_t[, keep_taxa]

# Calculate Spearman correlation and p-values
cor_res <- rcorr(as.matrix(data_filtered), type = "spearman")
cor_mat <- cor_res$r
p_mat <- cor_res$P

# Filter for strong, significant correlations
cor_thresh <- 0.6
p_thresh <- 0.05
sig_mask <- (cor_mat >= cor_thresh) & (p_mat < p_thresh)
diag(sig_mask) <- FALSE

# Create edge list
edges <- which(sig_mask, arr.ind = TRUE)
edge_df <- data.frame(
  from = rownames(cor_mat)[edges[,1]],
  to = colnames(cor_mat)[edges[,2]],
  weight = cor_mat[edges]
) %>% filter(from < to)  # remove duplicates

# Build graph
g <- graph_from_data_frame(edge_df, directed = FALSE)


# Convert igraph object to tidygraph
tg <- as_tbl_graph(g)

# Set node degree as size
tg <- tg %>%
  mutate(degree = centrality_degree())

# Set edge sign
E(tg)$sign <- ifelse(E(tg)$weight > 0, "positive", "negative")

pdf("taxon_network_ggraph_HHS.pdf", width = 14, height = 10)

# Plot with ggraph
ggraph(tg, layout = "fr") +  # Try layout = "kk" or "fr" (Fruchterman-Reingold)
  geom_edge_link(aes(color = sign, width = abs(weight)), alpha = 0.6) +
  geom_node_point(aes(size = degree), color = "steelblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_edge_color_manual(values = c(positive = "forestgreen", negative = "firebrick")) +
  scale_edge_width(range = c(0.3, 2)) +
  theme_void() +
  labs(title = "Taxon Co-occurrence Network HHS_fecal",
       edge_color = "Correlation Sign") +
  theme(legend.position = "right")

dev.off()
```
