#!/usr/bin/env Rscript
# ============================================================================
# Allele Frequency Mapping with Real Geographic Data
# Author: Michael Kantar
# Date: December 2024
# Purpose: Create categorical genotype maps with proper continental outlines
# ============================================================================

# Load required packages
library(vcfR)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

# ============================================================================
# CONFIGURATION
# ============================================================================

VCF_FILE <- "IWG_PI.vcf"
COORDS_FILE <- "IWG_lat_long_V2.csv"
MARKERS_FILE <- "top_5_IWG_markers.csv"
OUTPUT_DIR <- "IWG_final_outputs"

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ============================================================================
# READ DATA
# ============================================================================

cat("\n=== Allele Frequency Mapping with Geographic Context ===\n\n")
cat("Step 1: Reading input files...\n")

# Read VCF
vcf <- read.vcfR(VCF_FILE, verbose = FALSE)
cat(sprintf("  VCF: %d variants, %d samples\n", nrow(vcf), ncol(vcf@gt) - 1))

# Read markers
markers_df <- read_csv(MARKERS_FILE, show_col_types = FALSE) %>%
  filter(!is.na(marker))
cat(sprintf("  Markers: %d\n", nrow(markers_df)))

# Read coordinates
coords_df <- read_csv(COORDS_FILE, show_col_types = FALSE) %>%
  filter(!is.na(Latitude) & !is.na(Longitude))
cat(sprintf("  Coordinates: %d accessions\n", nrow(coords_df)))

# ============================================================================
# EXTRACT GENOTYPES AS CATEGORICAL
# ============================================================================

cat("\nStep 2: Extracting categorical genotypes...\n")

# Function to classify genotypes
classify_genotype <- function(gt_vector) {
  gt_vector <- gsub("\\|", "/", gt_vector)  # Handle both separators
  
  classified <- case_when(
    gt_vector %in% c("0/0", "0|0") ~ "Ref/Ref",
    gt_vector %in% c("0/1", "1/0", "0|1", "1|0") ~ "Ref/Alt",
    gt_vector %in% c("1/1", "1|1") ~ "Alt/Alt",
    TRUE ~ "Missing"
  )
  
  return(classified)
}

# Extract genotypes for each marker
marker_names <- markers_df$marker
genotype_list <- list()

for (marker in marker_names) {
  marker_idx <- which(vcf@fix[, "ID"] == marker)
  
  if (length(marker_idx) == 0) {
    cat(sprintf("  Warning: %s not found\n", marker))
    next
  }
  
  if (length(marker_idx) > 1) {
    marker_idx <- marker_idx[1]
  }
  
  # Extract genotypes
  gt_data <- extract.gt(vcf, element = "GT")[marker_idx, , drop = FALSE]
  sample_names <- colnames(gt_data)
  
  # Classify each genotype
  genotypes <- classify_genotype(as.vector(gt_data))
  
  # Create data frame
  marker_df <- data.frame(
    ID = sample_names,
    marker = marker,
    genotype = genotypes,
    stringsAsFactors = FALSE
  )
  
  genotype_list[[marker]] <- marker_df
  cat(sprintf("  Processed: %s\n", marker))
}

# Combine all data
genotype_data <- bind_rows(genotype_list)

# ============================================================================
# MERGE WITH COORDINATES
# ============================================================================

cat("\nStep 3: Merging with coordinates...\n")

map_data <- genotype_data %>%
  left_join(coords_df, by = "ID") %>%
  filter(!is.na(Latitude) & !is.na(Longitude))

cat(sprintf("  Complete records: %d\n", nrow(map_data)))

# Print genotype distribution
cat("\n  Genotype distribution:\n")
map_data %>%
  group_by(genotype) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  print()

# ============================================================================
# GET WORLD MAP WITH REAL GEOGRAPHY
# ============================================================================

cat("\nStep 4: Loading world map data...\n")

# Get world countries from Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")
cat("  World map loaded successfully\n")

# Calculate extent
lat_range <- range(map_data$Latitude)
lon_range <- range(map_data$Longitude)
buffer <- 15  # degrees

map_extent <- c(
  xmin = lon_range[1] - buffer,
  xmax = lon_range[2] + buffer,
  ymin = lat_range[1] - buffer,
  ymax = lat_range[2] + buffer
)

cat(sprintf("  Map extent: Lon [%.1f, %.1f], Lat [%.1f, %.1f]\n",
            map_extent[1], map_extent[2], map_extent[3], map_extent[4]))

# ============================================================================
# CREATE CATEGORICAL GENOTYPE MAPS
# ============================================================================

cat("\nStep 5: Creating categorical genotype maps...\n")

# Define categorical colors
genotype_colors <- c(
  "Ref/Ref" = "#2E5EAA",   # Dark blue
  "Ref/Alt" = "#FFD700",   # Gold
  "Alt/Alt" = "#C41E3A",   # Cardinal red
  "Missing" = "#888888"    # Gray
)

# Function to create map with proper geography
create_geographic_map <- function(marker_name, data, world_map, extent) {
  
  # Filter data
  marker_data <- data %>% filter(marker == marker_name)
  
  # Get trait info
  trait_info <- markers_df %>% filter(marker == marker_name)
  subtitle <- if (nrow(trait_info) > 0 && !is.na(trait_info$trait_name[1])) {
    sprintf("Trait: %s | p-value: %.2e", 
            trait_info$trait_name[1], 
            trait_info$p_value[1])
  } else {
    ""
  }
  
  # Count genotypes
  genotype_counts <- marker_data %>%
    count(genotype) %>%
    arrange(match(genotype, c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing")))
  
  # Create plot
  p <- ggplot() +
    # Ocean background
    theme(panel.background = element_rect(fill = "#AED6F1", color = NA)) +
    
    # World map - land masses
    geom_sf(data = world_map, 
            fill = "#E8D4A2",        # Tan/sand color for land
            color = "#8B7355",       # Brown for borders
            size = 0.3) +
    
    # Set map extent
    coord_sf(xlim = c(extent[1], extent[2]),
             ylim = c(extent[3], extent[4]),
             expand = FALSE) +
    
    # Plot genotype points
    geom_point(data = marker_data,
               aes(x = Longitude, y = Latitude, 
                   color = genotype, shape = genotype),
               size = 4,
               alpha = 0.85,
               stroke = 1.2) +
    
    # Color scale - categorical
    scale_color_manual(name = "Genotype",
                      values = genotype_colors,
                      breaks = c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing"),
                      labels = sprintf("%s (n=%d)", 
                                     c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing"),
                                     genotype_counts$n[match(c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing"), 
                                                            genotype_counts$genotype)])) +
    
    # Shape scale
    scale_shape_manual(name = "Genotype",
                      values = c("Ref/Ref" = 16, "Ref/Alt" = 17, 
                                "Alt/Alt" = 15, "Missing" = 4),
                      breaks = c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing"),
                      labels = sprintf("%s (n=%d)", 
                                     c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing"),
                                     genotype_counts$n[match(c("Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing"), 
                                                            genotype_counts$genotype)])) +
    
    # Labels
    labs(title = marker_name,
         subtitle = subtitle,
         x = "Longitude (°E)",
         y = "Latitude (°N)") +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major = element_line(color = "white", size = 0.5, linetype = "dashed"),
      panel.border = element_rect(fill = NA, color = "black", size = 1.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11)
    )
  
  return(p)
}

# Create individual maps
unique_markers <- unique(map_data$marker)
map_list <- list()

for (marker in unique_markers) {
  cat(sprintf("  Creating map: %s\n", marker))
  
  p <- create_geographic_map(marker, map_data, world, map_extent)
  map_list[[marker]] <- p
  
  # Save
  filename <- file.path(OUTPUT_DIR, 
                        sprintf("geographic_categorical_%s.png", 
                                gsub("[^A-Za-z0-9_]", "_", marker)))
  ggsave(filename, p, width = 14, height = 10, dpi = 300, bg = "white")
  cat(sprintf("    Saved: %s\n", filename))
}

# ============================================================================
# CREATE COMBINED PANEL
# ============================================================================

cat("\nStep 6: Creating combined panel...\n")

if (length(map_list) > 1) {
  combined_plot <- plot_grid(plotlist = map_list,
                              ncol = 2,
                              align = "hv",
                              labels = "AUTO",
                              label_size = 14)
  
  combined_file <- file.path(OUTPUT_DIR, "geographic_categorical_combined.png")
  ggsave(combined_file, combined_plot,
         width = 20, height = 12 * ceiling(length(map_list)/2),
         dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", combined_file))
}

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\nStep 7: Generating summary statistics...\n")

# Genotype counts by marker
genotype_summary <- map_data %>%
  group_by(marker, genotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = count, values_fill = 0)

# Add trait information
summary_stats <- genotype_summary %>%
  left_join(markers_df %>% select(marker, trait_name, p_value, effect), 
            by = "marker")

# Reorder columns
col_order <- c("marker", "trait_name", "Ref/Ref", "Ref/Alt", "Alt/Alt", "Missing", "p_value", "effect")
col_order <- col_order[col_order %in% names(summary_stats)]
summary_stats <- summary_stats[, col_order]

# Save
summary_file <- file.path(OUTPUT_DIR, "genotype_summary_categorical.csv")
write_csv(summary_stats, summary_file)

cat("\nGenotype Summary:\n")
print(summary_stats, n = Inf)

# ============================================================================
# COMPLETION
# ============================================================================

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nGenerated files:\n")
cat("  - Individual maps: geographic_categorical_*.png\n")
cat("  - Combined panel: geographic_categorical_combined.png\n")
cat("  - Summary statistics: genotype_summary_categorical.csv\n")
cat(sprintf("\nMarkers analyzed: %d\n", length(unique_markers)))
cat(sprintf("Total data points: %d\n", nrow(map_data)))
cat("\nDone!\n\n")
