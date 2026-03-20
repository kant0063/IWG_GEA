#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)

# Read in the data
gebvs <- read.csv("/mnt/user-data/uploads/GEBVs_GBLUP-Gaussian_all.csv", row.names = 1)
pca_clusters <- read.csv("/mnt/user-data/uploads/IWG_pca_clusters_k6.csv")
lat_long <- read.csv("/mnt/user-data/uploads/IWG_lat_long_V2.csv")

cat("Loading and processing data...\n")

# Prepare GEBV matrix
gebvs_mat <- as.matrix(gebvs)

# Match accessions between datasets
common_accessions <- intersect(rownames(gebvs_mat), pca_clusters$Accession)
gebvs_mat <- gebvs_mat[common_accessions, ]
pca_clusters_matched <- pca_clusters[match(common_accessions, pca_clusters$Accession), ]

# Create annotation for PCA clusters
annotation_pca <- data.frame(
  Accession = common_accessions,
  Cluster = as.factor(pca_clusters_matched$Cluster),
  row.names = common_accessions
)

# Process lat/long data
lat_long_clean <- lat_long[!is.na(lat_long$Latitude) & !is.na(lat_long$Longitude), ]

# Function to assign country/region based on coordinates
assign_country <- function(lat, lon) {
  if (is.na(lat) || is.na(lon)) return(NA)
  if (lat >= 36 && lat <= 42 && lon >= 26 && lon <= 45) return("Turkey")
  if (lat >= 39 && lat <= 43 && lon >= 43 && lon <= 47) return("Armenia/Georgia")
  if (lat >= 35 && lat <= 43 && lon >= 58 && lon <= 75) return("Central Asia")
  if (lat >= 42 && lat <= 56 && lon >= 48 && lon <= 90) return("Russia/Kazakhstan")
  if (lat >= 29 && lat <= 40 && lon >= 44 && lon <= 63) return("Iran/Iraq")
  return("Other")
}

# Apply country assignment
lat_long_clean$Country <- mapply(assign_country, lat_long_clean$Latitude, lat_long_clean$Longitude)

# Match with GEBVs data for country analysis
accessions_with_coords <- intersect(common_accessions, lat_long_clean$X...ID)
gebvs_mat_coords <- gebvs_mat[accessions_with_coords, ]
lat_long_matched <- lat_long_clean[match(accessions_with_coords, lat_long_clean$X...ID), ]

# Create annotation for countries
annotation_country <- data.frame(
  Accession = accessions_with_coords,
  Country = as.factor(lat_long_matched$Country),
  row.names = accessions_with_coords
)

# Scale the data (z-score normalization)
gebvs_scaled <- t(scale(t(gebvs_mat)))
gebvs_scaled_coords <- t(scale(t(gebvs_mat_coords)))

# Handle any remaining non-finite values
gebvs_scaled[!is.finite(gebvs_scaled)] <- 0
gebvs_scaled_coords[!is.finite(gebvs_scaled_coords)] <- 0

cat("Data processing complete.\n")
cat("Creating enhanced ggplot2 heatmaps...\n\n")

# =====================================
# HEATMAP 1: Grouped by PCA Cluster
# =====================================

# Order accessions by cluster
cluster_order <- order(annotation_pca$Cluster)
gebvs_ordered_cluster <- gebvs_scaled[cluster_order, ]
annotation_ordered_cluster <- annotation_pca[cluster_order, ]

# Convert to long format for ggplot
gebvs_long_cluster <- melt(gebvs_ordered_cluster)
colnames(gebvs_long_cluster) <- c("Accession", "Trait", "GEBV")

# Add cluster information and position
gebvs_long_cluster$Cluster <- annotation_ordered_cluster$Cluster[
  match(gebvs_long_cluster$Accession, annotation_ordered_cluster$Accession)
]

# Create ordered factor for plotting
gebvs_long_cluster$Accession <- factor(gebvs_long_cluster$Accession, 
                                       levels = rownames(gebvs_ordered_cluster))
gebvs_long_cluster$AccessionNum <- as.numeric(gebvs_long_cluster$Accession)

# Define cluster colors
cluster_colors <- setNames(brewer.pal(6, "Set2"), levels(annotation_pca$Cluster))

# Create annotation dataframe for sidebar
cluster_annotation <- data.frame(
  AccessionNum = 1:nrow(gebvs_ordered_cluster),
  Cluster = annotation_ordered_cluster$Cluster
)

# Create the main heatmap
p1_main <- ggplot(gebvs_long_cluster, aes(x = Trait, y = AccessionNum, fill = GEBV)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdBu")),
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(
    title = "GEBVs Heatmap Grouped by PCA Cluster",
    subtitle = sprintf("%d accessions × %d traits", nrow(gebvs_ordered_cluster), ncol(gebvs_ordered_cluster)),
    x = "Trait",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 5)
  ) +
  # Add cluster boundaries
  geom_hline(yintercept = cumsum(table(annotation_ordered_cluster$Cluster)) + 0.5, 
             color = "white", linewidth = 1.5)

# Create the annotation sidebar
p1_anno <- ggplot(cluster_annotation, aes(x = 1, y = AccessionNum, fill = Cluster)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = cluster_colors, name = "PCA\nCluster") +
  labs(x = "", y = "Accessions") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 0, 10, 10),
    legend.position = "left"
  ) +
  geom_hline(yintercept = cumsum(table(annotation_ordered_cluster$Cluster)) + 0.5, 
             color = "white", linewidth = 1.5)

# Combine plots
library(gridExtra)
p1_combined <- grid.arrange(p1_anno, p1_main, ncol = 2, widths = c(0.5, 6))

# Save the plot
ggsave("/mnt/user-data/outputs/heatmap_cluster_ggplot2.pdf", 
       plot = p1_combined, width = 18, height = 12, dpi = 300)
cat("✓ Created: heatmap_cluster_ggplot2.pdf\n")

# =====================================
# HEATMAP 2: Grouped by Country
# =====================================

# Order accessions by country
country_order <- order(annotation_country$Country)
gebvs_ordered_country <- gebvs_scaled_coords[country_order, ]
annotation_ordered_country <- annotation_country[country_order, ]

# Convert to long format for ggplot
gebvs_long_country <- melt(gebvs_ordered_country)
colnames(gebvs_long_country) <- c("Accession", "Trait", "GEBV")

# Add country information
gebvs_long_country$Country <- annotation_ordered_country$Country[
  match(gebvs_long_country$Accession, annotation_ordered_country$Accession)
]

# Create ordered factor for plotting
gebvs_long_country$Accession <- factor(gebvs_long_country$Accession, 
                                       levels = rownames(gebvs_ordered_country))
gebvs_long_country$AccessionNum <- as.numeric(gebvs_long_country$Accession)

# Define country colors
n_countries <- length(unique(annotation_country$Country))
country_colors <- setNames(brewer.pal(max(3, n_countries), "Set1")[1:n_countries], 
                           levels(annotation_country$Country))

# Create annotation dataframe for sidebar
country_annotation <- data.frame(
  AccessionNum = 1:nrow(gebvs_ordered_country),
  Country = annotation_ordered_country$Country
)

# Create the main heatmap
p2_main <- ggplot(gebvs_long_country, aes(x = Trait, y = AccessionNum, fill = GEBV)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "RdBu")),
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Z-score"
  ) +
  labs(
    title = "GEBVs Heatmap Grouped by Geographic Region",
    subtitle = sprintf("%d accessions × %d traits", nrow(gebvs_ordered_country), ncol(gebvs_ordered_country)),
    x = "Trait",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 5)
  ) +
  # Add country boundaries
  geom_hline(yintercept = cumsum(table(annotation_ordered_country$Country)) + 0.5, 
             color = "white", linewidth = 1.5)

# Create the annotation sidebar
p2_anno <- ggplot(country_annotation, aes(x = 1, y = AccessionNum, fill = Country)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = country_colors, name = "Geographic\nRegion") +
  labs(x = "", y = "Accessions") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 0, 10, 10),
    legend.position = "left",
    legend.text = element_text(size = 9)
  ) +
  geom_hline(yintercept = cumsum(table(annotation_ordered_country$Country)) + 0.5, 
             color = "white", linewidth = 1.5)

# Combine plots
p2_combined <- grid.arrange(p2_anno, p2_main, ncol = 2, widths = c(0.5, 6))

# Save the plot
ggsave("/mnt/user-data/outputs/heatmap_country_ggplot2.pdf", 
       plot = p2_combined, width = 18, height = 12, dpi = 300)
cat("✓ Created: heatmap_country_ggplot2.pdf\n")

# =====================================
# Summary Statistics
# =====================================

cat("\n==========================================\n")
cat("      ENHANCED GGPLOT2 HEATMAPS\n")
cat("==========================================\n\n")
cat("Dataset Information:\n")
cat("  • Total accessions with GEBVs:", nrow(gebvs_mat), "\n")
cat("  • Accessions with PCA data:", nrow(annotation_pca), "\n")
cat("  • Accessions with coordinates:", nrow(annotation_country), "\n")
cat("  • Number of traits analyzed:", ncol(gebvs_mat), "\n\n")

cat("PCA Cluster Distribution:\n")
cluster_table <- table(annotation_pca$Cluster)
for(i in 1:length(cluster_table)) {
  cat(sprintf("  • Cluster %s: %d accessions (%.1f%%)\n", 
              names(cluster_table)[i], 
              cluster_table[i],
              100 * cluster_table[i] / sum(cluster_table)))
}

cat("\nGeographic Region Distribution:\n")
country_table <- table(annotation_country$Country)
for(i in 1:length(country_table)) {
  cat(sprintf("  • %s: %d accessions (%.1f%%)\n", 
              names(country_table)[i], 
              country_table[i],
              100 * country_table[i] / sum(country_table)))
}

cat("\n==========================================\n")
cat("Enhanced ggplot2 heatmaps created!\n")
cat("==========================================\n")
cat("Files saved:\n")
cat("  1. heatmap_cluster_ggplot2.pdf\n")
cat("  2. heatmap_country_ggplot2.pdf\n")
cat("\nFeatures:\n")
cat("  • Color-coded annotation sidebars\n")
cat("  • White lines separating groups\n")
cat("  • Z-score normalized values\n")
cat("  • High-resolution output (300 DPI)\n")
cat("==========================================\n")
