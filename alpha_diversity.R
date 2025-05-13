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
