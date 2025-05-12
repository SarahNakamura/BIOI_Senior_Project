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

