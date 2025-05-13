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

