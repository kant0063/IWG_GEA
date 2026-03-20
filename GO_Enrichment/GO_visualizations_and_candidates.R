# GO Enrichment Analysis Visualizations and Candidate Gene Extraction
# Updated with correct file paths

library(ggplot2)
library(dplyr)
library(tidyr)

# Set correct results directory (use full path from your working directory)
results_dir <- "Analysis_Outputs_and_Visualizations/GO_Analysis_Results/"

# Create output directory for figures if it doesn't exist
output_dir <- "Analysis_Outputs_and_Visualizations/GO_Figures/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load GO enrichment results
cat("Loading GO enrichment results...\n")
all_traits <- read.csv(file.path(results_dir, "GO_enrichment_results_all_traits.csv"))
temp_results <- read.csv(file.path(results_dir, "GO_enrichment_temperature_group.csv"))
precip_results <- read.csv(file.path(results_dir, "GO_enrichment_precipitation_group.csv"))

# Rename columns to match expected names
# Your data has: FDR (not p.adjust), Significant (not Count), ontology (not Ont)
names(all_traits)[names(all_traits) == "FDR"] <- "p.adjust"
names(all_traits)[names(all_traits) == "Significant"] <- "Count"
names(all_traits)[names(all_traits) == "ontology"] <- "Ont"

names(temp_results)[names(temp_results) == "FDR"] <- "p.adjust"
names(temp_results)[names(temp_results) == "Significant"] <- "Count"
names(temp_results)[names(temp_results) == "ontology"] <- "Ont"

names(precip_results)[names(precip_results) == "FDR"] <- "p.adjust"
names(precip_results)[names(precip_results) == "Significant"] <- "Count"
names(precip_results)[names(precip_results) == "ontology"] <- "Ont"

# Convert to numeric
all_traits$p.adjust <- as.numeric(all_traits$p.adjust)
all_traits$Count <- as.numeric(all_traits$Count)

temp_results$p.adjust <- as.numeric(temp_results$p.adjust)
temp_results$Count <- as.numeric(temp_results$Count)

precip_results$p.adjust <- as.numeric(precip_results$p.adjust)
precip_results$Count <- as.numeric(precip_results$Count)

cat("Loaded:\n")
cat("  All traits:", nrow(all_traits), "GO terms\n")
cat("  Temperature:", nrow(temp_results), "GO terms\n")
cat("  Precipitation:", nrow(precip_results), "GO terms\n\n")

# Function to create bar plot of top GO terms
plot_top_go_terms <- function(go_data, title, top_n = 20, output_file = NULL) {
  
  # Filter significant terms and select top N by p-value
  sig_terms <- go_data %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(Term = reorder(Term, -log10(p.adjust)))
  
  if (nrow(sig_terms) == 0) {
    cat("No significant GO terms found for", title, "\n")
    return(NULL)
  }
  
  # Create plot
  p <- ggplot(sig_terms, aes(x = Term, y = -log10(p.adjust), fill = Ont)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = title,
         x = "GO Term",
         y = "-log10(Adjusted P-value)",
         fill = "Ontology") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    cat("Saved plot to:", output_file, "\n")
  }
  
  return(p)
}

# Function to create dot plot
plot_go_dotplot <- function(go_data, title, top_n = 20, output_file = NULL) {
  
  sig_terms <- go_data %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    mutate(Term = reorder(Term, -log10(p.adjust)))
  
  if (nrow(sig_terms) == 0) {
    cat("No significant GO terms found for", title, "\n")
    return(NULL)
  }
  
  p <- ggplot(sig_terms, aes(x = -log10(p.adjust), y = Term)) +
    geom_point(aes(size = Count, color = Ont)) +
    labs(title = title,
         x = "-log10(Adjusted P-value)",
         y = "GO Term",
         size = "Gene Count",
         color = "Ontology") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)
    cat("Saved plot to:", output_file, "\n")
  }
  
  return(p)
}

# Function to extract candidate genes
extract_candidate_genes <- function(go_data, output_file = NULL) {
  
  sig_terms <- go_data %>%
    filter(p.adjust < 0.05)
  
  if (nrow(sig_terms) == 0) {
    cat("No significant GO terms found\n")
    return(NULL)
  }
  
  # Check if geneID column exists
  if (!"geneID" %in% colnames(go_data)) {
    cat("Note: geneID column not found in data.\n")
    cat("Your GO enrichment results contain", nrow(sig_terms), "significant GO terms.\n")
    cat("To extract candidate genes, you'll need to map these GO terms back to genes.\n\n")
    
    # Return the significant GO terms instead (using base R subsetting)
    sig_summary <- sig_terms[, c("GO.ID", "Term", "Ont", "Count", "p.adjust")]
    sig_summary <- sig_summary[order(sig_summary$p.adjust), ]
    
    if (!is.null(output_file)) {
      write.csv(sig_summary, output_file, row.names = FALSE)
      cat("Saved significant GO terms to:", output_file, "\n")
    }
    
    return(sig_summary)
  }
  
  # Extract unique genes from geneID column
  all_genes <- unique(unlist(strsplit(sig_terms$geneID, "/")))
  
  cat("Found", length(all_genes), "unique candidate genes\n")
  
  # Create summary table
  candidate_summary <- data.frame(
    Gene = all_genes,
    stringsAsFactors = FALSE
  )
  
  # Add count of how many GO terms each gene appears in
  candidate_summary$GO_term_count <- sapply(candidate_summary$Gene, function(gene) {
    sum(grepl(gene, sig_terms$geneID, fixed = TRUE))
  })
  
  # Sort by frequency
  candidate_summary <- candidate_summary %>%
    arrange(desc(GO_term_count))
  
  if (!is.null(output_file)) {
    write.csv(candidate_summary, output_file, row.names = FALSE)
    cat("Saved candidate genes to:", output_file, "\n")
  }
  
  return(candidate_summary)
}

# Generate visualizations and extract candidates

cat("\n=== GENERATING VISUALIZATIONS ===\n\n")

# All traits
cat("Processing all traits...\n")
plot_top_go_terms(all_traits, 
                  "Top 20 GO Terms - All Climate Traits",
                  output_file = file.path(output_dir, "GO_barplot_all_traits.png"))

plot_go_dotplot(all_traits,
                "GO Enrichment - All Climate Traits",
                output_file = file.path(output_dir, "GO_dotplot_all_traits.png"))

# Temperature
cat("\nProcessing temperature traits...\n")
plot_top_go_terms(temp_results,
                  "Top 20 GO Terms - Temperature Adaptation",
                  output_file = file.path(output_dir, "GO_barplot_temperature.png"))

plot_go_dotplot(temp_results,
                "GO Enrichment - Temperature Adaptation",
                output_file = file.path(output_dir, "GO_dotplot_temperature.png"))

# Precipitation
cat("\nProcessing precipitation traits...\n")
plot_top_go_terms(precip_results,
                  "Top 20 GO Terms - Precipitation Adaptation",
                  output_file = file.path(output_dir, "GO_barplot_precipitation.png"))

plot_go_dotplot(precip_results,
                "GO Enrichment - Precipitation Adaptation",
                output_file = file.path(output_dir, "GO_dotplot_precipitation.png"))

cat("\n=== EXTRACTING CANDIDATE GENES ===\n\n")

# Extract candidate genes for each analysis
all_genes <- extract_candidate_genes(all_traits,
                                     file.path(results_dir, "candidate_genes_all_traits.csv"))

temp_genes <- extract_candidate_genes(temp_results,
                                      file.path(results_dir, "candidate_genes_temperature.csv"))

precip_genes <- extract_candidate_genes(precip_results,
                                        file.path(results_dir, "candidate_genes_precipitation.csv"))

cat("\n=== SUMMARY ===\n")
cat("Total candidate genes (all traits):", nrow(all_genes), "\n")
cat("Temperature-specific candidates:", nrow(temp_genes), "\n")
cat("Precipitation-specific candidates:", nrow(precip_genes), "\n")
cat("\nAll outputs saved to:", output_dir, "and", results_dir, "\n")
