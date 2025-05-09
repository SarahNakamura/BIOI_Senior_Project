# BIOI_Senior_Project
This is the code repository for my bioinformatics senior project 
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
```
module load igs-portal-client/1.4
```
## Quality Control
### 1.FastQC
Bash scripts (Since the sequence files were very large, I wrote three separate bash scripts to run them separately):
- [FastQC_AO](fastqc_AO.slurm)
- [FastQC_EA](fastqc_EA.slurm)
- [FastQC_HHS](fastqc_HHS.slurm)
### 2.MultiQC
After running FastQC, we could run MiltiQC to visualize the sequence quality. For this, I ran
```
module load multiqc/py37/1.8
multiqc -o path/to/output_dir path/to/input_dir
```
## Adapter Removal
From the MultiQC report, there seem to be conatmination around 10 bases in the 3' end of most sequence data. Therefore, I ran cutadapt to remove the low sequnece quality region from all seuquence data.
Bash script can be found here.
## Host Contamination Removal
## Contamination Removal
## Taxonomic Classification
## Taxonomic Abunance Analysis
The abundance analysis was performed using bracken. The slurm scripts are as follows:
- [Bracken_AO](bracken2_bacteria_AO.slurm)
- [Bracken_EA](bracken2_bacteria_EA.slurm)
- [Bracke_HHS](bracken2_bacteria_HHS.slurm)
## Diversity Analysis (apha diversity & beta diversity)
Python scripts from KrakenTools were used in both analyses, alpha diversity and beta diversity. The GitHub repository for KrakenTools can be found [here](https://github.com/jenniferlu717/KrakenTools)
The slurm scripts for alpha diversity analysis is as follows:
- [Alpha_diversity_AO](alpha_diversity_AO.slurm)
- [Alpha_diversity_EA](alpha_diversity_EA.slurm)
- [Alpha_diversity_HHS](alpha_diversity_HHS.slurm)
It takes the abundance results from bracken as input and outputs a csv file that contain the alpha diversity metrics.

The slurm scripts for beta diversity analysis is as follows:
- [Beta_diversity_AO](beta_diversity_AO.slurm)
- [Beta_diversity_EA](beta_diversity_EA.slurm)
- [Beta_diversity_HHS](beta_diversity_HHS.slurm)
The input files are the taxonomic abundance data from bracken. The output 
## Network Analysis
