# Load required libraries
library(vcfR)
library(genetics)
library(ggplot2)
library(dplyr)
library(parallel)
library(MetBrewer)

# Read VCF file
vcf <- read.vcfR("IWG_PI.vcf")

# Extract genotype matrix and convert to numeric
gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# Get chromosome and position information
chrom_info <- data.frame(
  chromosome = vcf@fix[, "CHROM"],
  position = as.numeric(vcf@fix[, "POS"]),
  snp_id = 1:nrow(vcf@fix)
)

# Function to calculate r-squared between two SNP vectors
calculate_r2 <- function(snp1, snp2) {
  # Remove individuals with missing data for either SNP
  complete_cases <- complete.cases(snp1, snp2)
  
  if (sum(complete_cases) < 10) return(NA)  # Need at least 10 individuals
  
  snp1_clean <- snp1[complete_cases]
  snp2_clean <- snp2[complete_cases]
  
  # Calculate correlation and square it for r²
  correlation <- cor(snp1_clean, snp2_clean, use = "complete.obs")
  r2 <- correlation^2
  
  return(r2)
}

# Function to estimate LD decay for a single chromosome
estimate_ld_decay_chromosome <- function(chrom_name, max_distance = 1000000, sample_size = 1000) {
  cat("Processing chromosome:", chrom_name, "\n")
  
  # Get SNPs for this chromosome
  chrom_snps <- chrom_info[chrom_info$chromosome == chrom_name, ]
  
  if (nrow(chrom_snps) < 10) {
    cat("Skipping chromosome", chrom_name, "- insufficient SNPs\n")
    return(NULL)
  }
  
  # Get genotype data for this chromosome
  chrom_gt <- gt_matrix[chrom_snps$snp_id, ]
  chrom_positions <- chrom_snps$position
  
  # Initialize results storage
  ld_results <- data.frame()
  
  # Define distance bins for analysis
  distance_bins <- c(0, 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000)
  
  # Sample SNP pairs to make analysis computationally feasible
  n_snps <- nrow(chrom_gt)
  
  # Create a sample of SNP pairs
  if (n_snps > 50) {
    # Sample SNPs for analysis
    sampled_indices <- sample(1:n_snps, min(n_snps, 200))  # Limit to 200 SNPs per chromosome
    sampled_gt <- chrom_gt[sampled_indices, ]
    sampled_positions <- chrom_positions[sampled_indices]
    n_sampled <- length(sampled_indices)
  } else {
    sampled_gt <- chrom_gt
    sampled_positions <- chrom_positions
    n_sampled <- n_snps
  }
  
  # Calculate LD for sampled SNP pairs
  pair_count <- 0
  max_pairs <- min(sample_size, (n_sampled * (n_sampled - 1)) / 2)
  
  for (i in 1:(n_sampled - 1)) {
    for (j in (i + 1):n_sampled) {
      if (pair_count >= max_pairs) break
      
      distance <- abs(sampled_positions[j] - sampled_positions[i])
      
      if (distance <= max_distance && distance > 0) {
        r2 <- calculate_r2(sampled_gt[i, ], sampled_gt[j, ])
        
        if (!is.na(r2)) {
          ld_results <- rbind(ld_results, data.frame(
            chromosome = chrom_name,
            distance = distance,
            r2 = r2,
            snp1_pos = sampled_positions[i],
            snp2_pos = sampled_positions[j]
          ))
          pair_count <- pair_count + 1
        }
      }
    }
    if (pair_count >= max_pairs) break
  }
  
  cat("Chromosome", chrom_name, ": calculated", nrow(ld_results), "LD values\n")
  return(ld_results)
}

# Get unique chromosomes
unique_chromosomes <- unique(chrom_info$chromosome)
cat("Found chromosomes:", paste(unique_chromosomes, collapse = ", "), "\n")

# Process each chromosome
all_ld_results <- data.frame()

for (chrom in unique_chromosomes) {
  chrom_results <- estimate_ld_decay_chromosome(chrom)
  if (!is.null(chrom_results)) {
    all_ld_results <- rbind(all_ld_results, chrom_results)
  }
}

# Save raw results
write.csv(all_ld_results, "ld_decay_raw_results.csv", row.names = FALSE)

# Create distance bins for summary statistics
all_ld_results$distance_bin <- cut(all_ld_results$distance, 
                                   breaks = c(0, 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000),
                                   labels = c("0-1kb", "1-5kb", "5-10kb", "10-25kb", "25-50kb", 
                                              "50-100kb", "100-250kb", "250-500kb", "500kb-1Mb"),
                                   include.lowest = TRUE)

# Calculate summary statistics by distance bin and chromosome
ld_summary <- all_ld_results %>%
  group_by(chromosome, distance_bin) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    median_r2 = median(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    n_pairs = n(),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(distance_bin))

# Save summary results
write.csv(ld_summary, "ld_decay_summary.csv", row.names = FALSE)

# Generate smooth decay curves and prepare data for plotting
create_smooth_decay_data <- function() {
  smooth_data <- data.frame()
  
  if (nrow(decay_models$parameters) > 0) {
    for (i in 1:nrow(decay_models$parameters)) {
      chrom <- decay_models$parameters$chromosome[i]
      params <- decay_models$parameters[i, ]
      
      # Create sequence of distances for smooth curve
      max_dist <- max(ld_summary$mean_distance[ld_summary$chromosome == chrom], na.rm = TRUE)
      distances <- seq(1000, max_dist, length.out = 100)
      
      # Calculate predicted r² values using fitted model
      predicted_r2 <- params$a * exp(-params$b * distances) + params$c
      
      smooth_data <- rbind(smooth_data, data.frame(
        chromosome = chrom,
        distance = distances,
        predicted_r2 = predicted_r2,
        half_decay_distance = params$half_decay_distance
      ))
    }
  }
  return(smooth_data)
}

smooth_decay_data <- create_smooth_decay_data()

# Calculate genome-wide average LD decay
genome_wide_summary <- all_ld_results %>%
  mutate(distance_bin = cut(distance, 
                            breaks = c(0, 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000),
                            labels = c("0-1kb", "1-5kb", "5-10kb", "10-25kb", "25-50kb", 
                                       "50-100kb", "100-250kb", "250-500kb", "500kb-1Mb"),
                            include.lowest = TRUE)) %>%
  filter(!is.na(distance_bin)) %>%
  group_by(distance_bin) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    median_r2 = median(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    n_pairs = n(),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(chromosome = "Genome-wide")

# Fit genome-wide decay model
genome_wide_model <- NULL
genome_wide_params <- NULL

if (nrow(genome_wide_summary) >= 4) {
  tryCatch({
    genome_wide_model <- nls(mean_r2 ~ a * exp(-b * mean_distance) + c, 
                             data = genome_wide_summary,
                             start = list(a = max(genome_wide_summary$mean_r2), b = 0.000001, c = 0.1),
                             control = nls.control(maxiter = 100, warnOnly = TRUE))
    
    params <- coef(genome_wide_model)
    genome_wide_params <- data.frame(
      chromosome = "Genome-wide",
      a = params["a"],
      b = params["b"], 
      c = params["c"],
      half_decay_distance = log(2) / params["b"]
    )
    
    cat("Genome-wide half-decay distance:", round(genome_wide_params$half_decay_distance), "bp\n")
    
  }, error = function(e) {
    cat("Could not fit genome-wide model:", e$message, "\n")
  })
}

# Create smooth decay data for genome-wide
genome_wide_smooth <- data.frame()
if (!is.null(genome_wide_params)) {
  max_dist <- max(genome_wide_summary$mean_distance, na.rm = TRUE)
  distances <- seq(1000, max_dist, length.out = 100)
  predicted_r2 <- genome_wide_params$a * exp(-genome_wide_params$b * distances) + genome_wide_params$c
  
  genome_wide_smooth <- data.frame(
    chromosome = "Genome-wide",
    distance = distances,
    predicted_r2 = predicted_r2,
    half_decay_distance = genome_wide_params$half_decay_distance
  )
}

# Create genome-wide LD decay plot
if (nrow(genome_wide_smooth) > 0) {
  genome_wide_plot <- ggplot() +
    # Add observed data points
    geom_point(data = genome_wide_summary, aes(x = mean_distance/1000, y = mean_r2), 
               color = "darkblue", size = 3, alpha = 0.7) +
    # Add smooth exponential decay curve
    geom_line(data = genome_wide_smooth, aes(x = distance/1000, y = predicted_r2), 
              color = "steelblue", size = 1.5) +
    # Add vertical line for half-decay distance
    geom_vline(xintercept = genome_wide_params$half_decay_distance/1000, 
               color = "red", linetype = "dashed", alpha = 0.8, size = 1) +
    # Add text annotation for half-decay distance
    annotate("text", 
             x = genome_wide_params$half_decay_distance/1000, 
             y = max(genome_wide_summary$mean_r2) * 0.8,
             label = paste0("Half-decay distance:\n", round(genome_wide_params$half_decay_distance), " bp"),
             hjust = -0.1, vjust = 0.5, size = 4, color = "red", fontface = "bold",
             fill = "white", alpha = 0.8) +
    scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")) +
    labs(
      title = "Genome-wide Average LD Decay",
      subtitle = paste0("Exponential decay curve with half-decay distance at ", 
                        round(genome_wide_params$half_decay_distance), " bp"),
      x = "Distance (kb, log scale)",
      y = "Mean r²"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      panel.grid.minor = element_blank()
    )
  
  print(genome_wide_plot)
  
  # Save genome-wide summary
  write.csv(genome_wide_summary, "ld_decay_genome_wide_summary.csv", row.names = FALSE)
  if (!is.null(genome_wide_params)) {
    write.csv(genome_wide_params, "ld_decay_genome_wide_parameters.csv", row.names = FALSE)
  }
}

# Create LD decay plot with smooth exponential curves
if (nrow(smooth_decay_data) > 0) {
  ld_decay_plot <- ggplot() +
    # Add observed data points
    geom_point(data = ld_summary, aes(x = mean_distance/1000, y = mean_r2, color = chromosome), 
               size = 2, alpha = 0.6) +
    # Add smooth exponential decay curves
    geom_line(data = smooth_decay_data, aes(x = distance/1000, y = predicted_r2, color = chromosome), 
              size = 1.2, alpha = 0.9) +
    # Add vertical lines for half-decay distances
    geom_vline(data = decay_models$parameters, 
               aes(xintercept = half_decay_distance/1000, color = chromosome), 
               linetype = "dashed", alpha = 0.7, size = 0.8) +
    scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")) +
    labs(
      title = "Linkage Disequilibrium Decay Across Chromosomes",
      subtitle = "Exponential decay curves with half-decay distance markers",
      x = "Distance (kb, log scale)",
      y = "Mean r²",
      color = "Chromosome"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    scale_color_manual(values = met.brewer("VanGogh1", n = length(unique(ld_summary$chromosome))))
  
  print(ld_decay_plot)
  
  # Create faceted plot with smooth curves and annotations
  faceted_ld_plot <- ggplot() +
    # Add observed data points
    geom_point(data = ld_summary, aes(x = mean_distance/1000, y = mean_r2), 
               color = "darkblue", size = 2, alpha = 0.6) +
    # Add smooth exponential decay curves
    geom_line(data = smooth_decay_data, aes(x = distance/1000, y = predicted_r2), 
              color = "steelblue", size = 1.2) +
    # Add vertical lines for half-decay distances
    geom_vline(data = decay_models$parameters, 
               aes(xintercept = half_decay_distance/1000), 
               color = "red", linetype = "dashed", alpha = 0.8, size = 0.8) +
    # Add text annotations for half-decay distances
    geom_text(data = decay_models$parameters, 
              aes(x = half_decay_distance/1000, y = Inf, 
                  label = paste0("Half-decay:\n", round(half_decay_distance), " bp")),
              vjust = 1.2, hjust = 0.5, size = 3, color = "red", fontface = "bold",
              bg.colour = "white", bg.r = 0.1) +
    facet_wrap(~ chromosome, scales = "free", ncol = 4) +
    scale_x_log10(breaks = c(1, 10, 100, 1000), labels = c("1", "10", "100", "1000")) +
    labs(
      title = "LD Decay by Chromosome with Half-Decay Distances",
      subtitle = "Red dashed line indicates half-decay distance",
      x = "Distance (kb, log scale)",
      y = "Mean r²"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
      axis.text = element_text(size = 8)
    )
  
  print(faceted_ld_plot)
  
} else {
  # Fallback to original plots if no models could be fitted
  ld_decay_plot <- ggplot(ld_summary, aes(x = mean_distance/1000, y = mean_r2, color = chromosome)) +
    geom_line(size = 1, alpha = 0.8) +
    geom_point(size = 2, alpha = 0.8) +
    scale_x_log10() +
    labs(
      title = "Linkage Disequilibrium Decay Across Chromosomes",
      x = "Distance (kb, log scale)",
      y = "Mean r²",
      color = "Chromosome"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    scale_color_manual(values = met.brewer("VanGogh1", n = length(unique(ld_summary$chromosome))))
  
  print(ld_decay_plot)
  
  faceted_ld_plot <- ggplot(ld_summary, aes(x = mean_distance/1000, y = mean_r2)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "darkblue", size = 2) +
    facet_wrap(~ chromosome, scales = "free", ncol = 4) +
    scale_x_log10() +
    labs(
      title = "LD Decay by Chromosome",
      x = "Distance (kb, log scale)",
      y = "Mean r²"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text = element_text(size = 8)
    )
  
  print(faceted_ld_plot)
}

# Fit exponential decay models for each chromosome
fit_decay_models <- function() {
  models <- list()
  decay_parameters <- data.frame()
  
  for (chrom in unique(ld_summary$chromosome)) {
    chrom_data <- ld_summary[ld_summary$chromosome == chrom & !is.na(ld_summary$mean_r2), ]
    
    if (nrow(chrom_data) >= 4) {  # Need at least 4 points to fit
      tryCatch({
        # Fit exponential decay model: r² = a * exp(-b * distance) + c
        model <- nls(mean_r2 ~ a * exp(-b * mean_distance) + c, 
                     data = chrom_data,
                     start = list(a = max(chrom_data$mean_r2), b = 0.000001, c = 0.1),
                     control = nls.control(maxiter = 100, warnOnly = TRUE))
        
        models[[chrom]] <- model
        
        # Extract parameters
        params <- coef(model)
        decay_parameters <- rbind(decay_parameters, data.frame(
          chromosome = chrom,
          a = params["a"],
          b = params["b"],
          c = params["c"],
          half_decay_distance = log(2) / params["b"]  # Distance where LD decays to half
        ))
        
      }, error = function(e) {
        cat("Could not fit model for chromosome", chrom, ":", e$message, "\n")
      })
    }
  }
  
  return(list(models = models, parameters = decay_parameters))
}

# Fit the decay models
decay_models <- fit_decay_models()

if (nrow(decay_models$parameters) > 0) {
  write.csv(decay_models$parameters, "ld_decay_parameters.csv", row.names = FALSE)
  
  # Print summary of decay parameters
  cat("\nLD Decay Parameters Summary:\n")
  print(decay_models$parameters)
  
  # Create a plot showing half-decay distances
  if (nrow(decay_models$parameters) > 1) {
    half_decay_plot <- ggplot(decay_models$parameters, aes(x = chromosome, y = half_decay_distance/1000)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      labs(
        title = "LD Half-Decay Distance by Chromosome",
        x = "Chromosome",
        y = "Half-Decay Distance (kb)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      )
    
    print(half_decay_plot)
    
    # Save plots with updated names
    ggsave("ld_decay_smooth_curves.png", ld_decay_plot, width = 14, height = 8, dpi = 300)
    ggsave("ld_decay_faceted_annotated.png", faceted_ld_plot, width = 16, height = 12, dpi = 300)
    ggsave("ld_half_decay_distances.png", half_decay_plot, width = 10, height = 6, dpi = 300)
    
    if (exists("genome_wide_plot")) {
      ggsave("ld_decay_genome_wide.png", genome_wide_plot, width = 10, height = 8, dpi = 300)
    }
    
    cat("Smooth exponential decay curves created with half-decay annotations!\n")
    cat("Genome-wide average LD decay plot also created!\n")
  }
}

# Print overall summary
cat("\nLD Decay Analysis Complete!\n")
cat("Total SNP pairs analyzed:", nrow(all_ld_results), "\n")
cat("Chromosomes processed:", length(unique(all_ld_results$chromosome)), "\n")
if (!is.null(genome_wide_params)) {
  cat("Genome-wide half-decay distance:", round(genome_wide_params$half_decay_distance), "bp\n")
}
cat("\nFiles created:\n")
cat("- ld_decay_raw_results.csv: Raw LD calculations\n")
cat("- ld_decay_summary.csv: Summary statistics by distance bins and chromosome\n")
cat("- ld_decay_parameters.csv: Exponential decay model parameters per chromosome\n")
cat("- ld_decay_genome_wide_summary.csv: Genome-wide average summary statistics\n")
if (!is.null(genome_wide_params)) {
  cat("- ld_decay_genome_wide_parameters.csv: Genome-wide decay model parameters\n")
}
cat("- PNG files: Visualization plots including genome-wide average\n")