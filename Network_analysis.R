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
