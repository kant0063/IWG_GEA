# Load required libraries
library(vcfR)
library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(MetBrewer)
library(corrplot)
library(carat)

# Read VCF file and extract genotypes (assuming this is already done from previous scripts)
if (!exists("vcf")) {
  vcf <- read.vcfR("IWG_PI.vcf")
}

if (!exists("gt_numeric")) {
  gt_numeric <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
}

# Load environmental data (assuming this is already done from previous scripts)
if (!exists("env_data_clean")) {
  env_data <- read.csv("IWG_bioclim.csv", header = TRUE, stringsAsFactors = FALSE)
  rows_to_remove <- c(1, 2, 3, 4, 5, 21, 39, 40, 42, 45, 46, 47, 48, 50, 
                      56, 57, 58, 59, 60, 62, 65, 70, 71, 81, 307, 318, 319, 320, 321, 325, 327, 328,                    
                      329, 330, 331)
  env_data_clean <- env_data[-rows_to_remove, ]
}

# Match genotype and environmental data (from previous scripts)
if (!exists("genotype_data") || !exists("csv_matched")) {
  vcf_samples <- colnames(gt_numeric)
  csv_samples <- env_data_clean$PI_Accession
  common_samples <- intersect(vcf_samples, csv_samples)
  
  gt_matched <- gt_numeric[, colnames(gt_numeric) %in% common_samples]
  csv_matched <- env_data_clean[env_data_clean$PI_Accession %in% common_samples, ]
  csv_matched <- csv_matched[match(colnames(gt_matched), csv_matched$PI_Accession), ]
  
  # Impute missing values in genotype data
  genotype_data <- gt_matched
  genotype_data <- apply(genotype_data, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
  genotype_data <- t(genotype_data)  # Transpose so samples are rows, SNPs are columns
}

# Prepare environmental data for RDA
exclude_cols <- c("PI_Accession", "Latitude", "Longitude")
env_vars <- csv_matched[, !colnames(csv_matched) %in% exclude_cols, drop = FALSE]

# Remove environmental variables with no variation or too many missing values
env_vars_clean <- env_vars[, sapply(env_vars, function(x) {
  var(x, na.rm = TRUE) > 0 & sum(!is.na(x)) > (nrow(env_vars) * 0.5)
})]

# Impute missing values in environmental data with column means
env_vars_clean <- apply(env_vars_clean, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})

# Check for and remove highly correlated environmental variables (|r| > 0.9)
env_cor <- cor(env_vars_clean, use = "complete.obs")

# Function to identify highly correlated variables
identify_highly_correlated <- function(cor_matrix, cutoff = 0.9) {
  # Get upper triangle of correlation matrix
  upper_tri <- upper.tri(cor_matrix)
  
  # Find pairs with high correlation
  high_cor_pairs <- which(abs(cor_matrix) > cutoff & upper_tri, arr.ind = TRUE)
  
  if (nrow(high_cor_pairs) == 0) {
    return(integer(0))
  }
  
  # For each highly correlated pair, remove the variable with higher mean correlation
  vars_to_remove <- c()
  
  for (i in 1:nrow(high_cor_pairs)) {
    var1 <- high_cor_pairs[i, 1]
    var2 <- high_cor_pairs[i, 2]
    
    # Calculate mean absolute correlation for each variable
    mean_cor_var1 <- mean(abs(cor_matrix[var1, -var1]), na.rm = TRUE)
    mean_cor_var2 <- mean(abs(cor_matrix[var2, -var2]), na.rm = TRUE)
    
    # Remove the variable with higher mean correlation
    if (mean_cor_var1 > mean_cor_var2) {
      vars_to_remove <- c(vars_to_remove, var1)
    } else {
      vars_to_remove <- c(vars_to_remove, var2)
    }
  }
  
  # Return unique variables to remove
  return(unique(vars_to_remove))
}

highly_corr <- identify_highly_correlated(env_cor, cutoff = 0.9)

if (length(highly_corr) > 0) {
  cat("Removing", length(highly_corr), "highly correlated environmental variables:\n")
  cat(colnames(env_vars_clean)[highly_corr], "\n")
  env_vars_final <- env_vars_clean[, -highly_corr]
} else {
  cat("No highly correlated variables found (|r| > 0.9)\n")
  env_vars_final <- env_vars_clean
}

# Scale environmental variables
env_vars_scaled <- scale(env_vars_final)

cat("Final dataset dimensions:\n")
cat("Samples:", nrow(genotype_data), "\n")
cat("SNPs:", ncol(genotype_data), "\n") 
cat("Environmental variables:", ncol(env_vars_scaled), "\n")

# Perform RDA
cat("\nPerforming RDA analysis...\n")
rda_result <- rda(genotype_data ~ ., data = as.data.frame(env_vars_scaled))

# Summary of RDA
cat("\nRDA Summary:\n")
print(summary(rda_result))

# Test significance of RDA model
cat("\nTesting significance of RDA model...\n")
rda_anova <- anova.cca(rda_result, permutations = 100)
print(rda_anova)

# Test significance of each environmental variable
cat("\nTesting significance of each environmental variable...\n")
rda_anova_by_axis <- anova.cca(rda_result, by = "terms", permutations = 100)
print(rda_anova_by_axis)

# Test significance of RDA axes
cat("\nTesting significance of RDA axes...\n")
rda_anova_by_axis <- anova.cca(rda_result, by = "axis", permutations = 100)
print(rda_anova_by_axis)

# Extract variance explained
total_var <- rda_result$tot.chi
constrained_var <- rda_result$CCA$tot.chi
unconstrained_var <- rda_result$CA$tot.chi

var_explained <- constrained_var / total_var
cat("\nVariance explained by environmental variables:", round(var_explained * 100, 2), "%\n")

# Extract RDA scores
site_scores <- scores(rda_result, choices = 1:min(4, rda_result$CCA$rank), display = "sites")
species_scores <- scores(rda_result, choices = 1:min(4, rda_result$CCA$rank), display = "species")
env_scores <- scores(rda_result, choices = 1:min(4, rda_result$CCA$rank), display = "bp")

# Create data frame for plotting
rda_df <- data.frame(
  Sample = rownames(site_scores),
  RDA1 = site_scores[, 1],
  RDA2 = site_scores[, 2]
)

if (ncol(site_scores) > 2) {
  rda_df$RDA3 <- site_scores[, 3]
}

# Add environmental data for coloring points
rda_df <- cbind(rda_df, env_vars_final[match(rda_df$Sample, rownames(env_vars_final)), ])

# Environmental arrows data frame
env_arrows <- data.frame(
  Variable = rownames(env_scores),
  RDA1 = env_scores[, 1],
  RDA2 = env_scores[, 2]
)

if (ncol(env_scores) > 2) {
  env_arrows$RDA3 <- env_scores[, 3]
}

# Calculate variance explained by each axis
axis_var <- eigenvals(rda_result) / total_var * 100

# Create RDA biplot for RDA1 vs RDA2
rda_plot <- ggplot(rda_df, aes(x = RDA1, y = RDA2)) +
  geom_point(size = 2, alpha = 0.7, color = "steelblue") +
  geom_segment(data = env_arrows, 
               aes(x = 0, y = 0, xend = RDA1 * 3, yend = RDA2 * 3),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "red", size = 1, alpha = 0.8) +
  geom_text(data = env_arrows, 
            aes(x = RDA1 * 3.2, y = RDA2 * 3.2, label = Variable),
            color = "red", size = 3, fontface = "bold") +
  labs(
    title = "RDA Biplot: Genomic Variation Constrained by Environment",
    subtitle = paste0("Total variance explained by environment: ", round(var_explained * 100, 2), "%"),
    x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
    y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  coord_fixed()

print(rda_plot)

# Create RDA plot colored by first environmental variable (if available)
if (ncol(env_vars_final) > 0) {
  first_env_var <- colnames(env_vars_final)[1]
  
  rda_plot_colored <- ggplot(rda_df, aes(x = RDA1, y = RDA2)) +
    geom_point(aes_string(color = first_env_var), size = 2.5, alpha = 0.8) +
    geom_segment(data = env_arrows, 
                 aes(x = 0, y = 0, xend = RDA1 * 3, yend = RDA2 * 3),
                 arrow = arrow(length = unit(0.3, "cm")), 
                 color = "black", size = 1, alpha = 0.8) +
    geom_text(data = env_arrows, 
              aes(x = RDA1 * 3.2, y = RDA2 * 3.2, label = Variable),
              color = "black", size = 3, fontface = "bold") +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = paste("RDA Biplot Colored by", first_env_var),
      subtitle = paste0("Total variance explained by environment: ", round(var_explained * 100, 2), "%"),
      x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
      y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    ) +
    coord_fixed()
  
  print(rda_plot_colored)
}

# Identify candidate SNPs (outliers) based on loadings
# Calculate loadings for SNPs on RDA axes
snp_loadings <- species_scores

# Calculate Mahalanobis distances to identify outlier SNPs
if (ncol(snp_loadings) >= 2) {
  # Use first two RDA axes for outlier detection
  loadings_subset <- snp_loadings[, 1:2]
  
  # Calculate Mahalanobis distances
  maha_dist <- mahalanobis(loadings_subset, 
                           center = colMeans(loadings_subset), 
                           cov = cov(loadings_subset))
  
  # Identify outliers using 99.9th percentile threshold
  outlier_threshold <- quantile(maha_dist, 0.999)
  outlier_snps <- which(maha_dist > outlier_threshold)
  
  cat("\nIdentified", length(outlier_snps), "candidate SNPs (outliers) associated with environmental variation\n")
  cat("Outlier threshold (99.9th percentile):", round(outlier_threshold, 3), "\n")
  
  # Create outlier SNPs data frame
  outlier_df <- data.frame(
    SNP_index = outlier_snps,
    RDA1_loading = snp_loadings[outlier_snps, 1],
    RDA2_loading = snp_loadings[outlier_snps, 2],
    Mahalanobis_distance = maha_dist[outlier_snps]
  )
  
  # Add chromosome and position information if available
  if (exists("vcf")) {
    outlier_df$chromosome <- vcf@fix[outlier_snps, "CHROM"]
    outlier_df$position <- as.numeric(vcf@fix[outlier_snps, "POS"])
  }
  
  # Save outlier SNPs
  write.csv(outlier_df, "rda_outlier_snps.csv", row.names = FALSE)
  
  # Create loading plot highlighting outliers
  loading_df <- data.frame(
    SNP = 1:nrow(snp_loadings),
    RDA1 = snp_loadings[, 1],
    RDA2 = snp_loadings[, 2],
    Outlier = 1:nrow(snp_loadings) %in% outlier_snps
  )
  
  loading_plot <- ggplot(loading_df, aes(x = RDA1, y = RDA2)) +
    geom_point(aes(color = Outlier), size = 1, alpha = 0.6) +
    scale_color_manual(values = c("gray70", "red"), 
                       labels = c("Background SNPs", "Candidate SNPs")) +
    labs(
      title = "SNP Loadings on RDA Axes",
      subtitle = paste("Candidate SNPs (n =", length(outlier_snps), ") highlighted in red"),
      x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
      y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom"
    ) +
    coord_fixed()
  
  print(loading_plot)
}

# Create environmental correlation plot
if (ncol(env_vars_final) > 1) {
  env_cor_final <- cor(env_vars_final, use = "complete.obs")
  
  png("environmental_correlations.png", width = 8, height = 8, units = "in", res = 300)
  corrplot(env_cor_final, method = "color", type = "upper", 
           tl.col = "black", tl.srt = 45, 
           title = "Environmental Variable Correlations",
           mar = c(0, 0, 2, 0))
  dev.off()
}

# Save RDA results and plots
save(rda_result, file = "rda_analysis_results.RData")

# Save plots
ggsave("rda_biplot.png", rda_plot, width = 10, height = 8, dpi = 300)
if (exists("rda_plot_colored")) {
  ggsave("rda_biplot_colored.png", rda_plot_colored, width = 10, height = 8, dpi = 300)
}
if (exists("loading_plot")) {
  ggsave("rda_snp_loadings.png", loading_plot, width = 10, height = 8, dpi = 300)
}

# Create summary table
rda_summary_table <- data.frame(
  Metric = c("Total Variance", "Constrained Variance", "Unconstrained Variance", 
             "Proportion Explained", "Number of Environmental Variables", 
             "Number of Samples", "Number of SNPs", "Candidate SNPs"),
  Value = c(
    round(total_var, 2),
    round(constrained_var, 2), 
    round(unconstrained_var, 2),
    paste0(round(var_explained * 100, 2), "%"),
    ncol(env_vars_scaled),
    nrow(genotype_data),
    ncol(genotype_data),
    ifelse(exists("outlier_snps"), length(outlier_snps), 0)
  )
)

write.csv(rda_summary_table, "rda_summary_table.csv", row.names = FALSE)

# Print final summary
cat("\n=== RDA Analysis Complete! ===\n")
cat("Files created:\n")
cat("- rda_analysis_results.RData: Complete RDA results object\n")
cat("- rda_biplot.png: Basic RDA biplot\n")
if (exists("rda_plot_colored")) {
  cat("- rda_biplot_colored.png: RDA biplot colored by environmental variable\n")
}
if (exists("loading_plot")) {
  cat("- rda_snp_loadings.png: SNP loadings plot with candidate SNPs highlighted\n")
  cat("- rda_outlier_snps.csv: List of candidate SNPs associated with environment\n")
}
cat("- environmental_correlations.png: Environmental variable correlation matrix\n")
cat("- rda_summary_table.csv: Summary statistics table\n")

cat("\nKey Results:\n")
cat("- Environmental variables explain", round(var_explained * 100, 2), "% of genomic variation\n")
if (exists("outlier_snps")) {
  cat("- Identified", length(outlier_snps), "candidate SNPs associated with environmental variation\n")
}
cat("- RDA model significance: p-value from ANOVA test shown above\n")