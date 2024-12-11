
# Load the data
ir_data <- read.csv("ir_gene_family_ranked.csv", row.names = 1)
obp_data <- read.csv("obp_gene_family_ranked.csv", row.names = 1)
or_data <- read.csv("or_gene_family_ranked.csv", row.names = 1)
# Install and load pvclust for hierarchical clustering with bootstrap resampling
library(pvclust)

# Function to plot dendrogram
plot_dendrogram <- function(data, title) {
  # Perform hierarchical clustering with complete linkage and Euclidean distance
  fit <- pvclust(data, method.hclust = "complete", method.dist = "euclidean", nboot = 1000)

  # Plot the dendrogram with bootstrap values
  plot(fit, main = title)
}

# Plot dendrograms for each gene family
plot_dendrogram(ir_data, "IR Gene Family")
plot_dendrogram(obp_data, "OBP Gene Family")
plot_dendrogram(or_data, "OR Gene Family")
