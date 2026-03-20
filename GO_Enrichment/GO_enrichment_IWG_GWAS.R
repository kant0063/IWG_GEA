################################################################################
# GO Enrichment Analysis for Intermediate Wheatgrass Climate Adaptation GWAS
# Author: Elizabeth (Lizzy)
# Project: IWG Chapter 2 - Climate Adaptation
################################################################################

# This script performs Gene Ontology (GO) enrichment analysis on GWAS results
# to identify biological processes involved in climate adaptation

################################################################################
# PART 1: SETUP AND PARAMETERS
################################################################################

# Install required packages (run once)
# install.packages(c("data.table", "dplyr", "ggplot2", "topGO", "stringr"))
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("topGO")

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(topGO)
library(stringr)

# Set working directory
setwd("C:/Users/Lizzy/OneDrive/Documents/UH_Manoa_Classes/IWG_Chapter_2")

# ============================================================================
# ADJUSTABLE PARAMETERS - CHANGE THESE AS NEEDED
# ============================================================================

# P-value threshold for significant SNPs
PVAL_THRESHOLD <- 0.01  # Try 0.001, 0.0001, or other values

# Window size around SNPs to search for genes (in base pairs)
WINDOW_SIZE <- 50000  # 25 kb upstream and downstream (try 25000, 50000, 100000)

# Minimum number of genes required for GO enrichment test
MIN_GENES <- 5  # Need at least this many genes for meaningful results

# FDR threshold for significant GO terms
GO_FDR_THRESHOLD <- 0.05

# Output directory
OUTPUT_DIR <- "GO_Analysis_Results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

cat("=================================================================\n")
cat("GO Enrichment Analysis Parameters:\n")
cat("=================================================================\n")
cat("P-value threshold:", PVAL_THRESHOLD, "\n")
cat("Window size:", WINDOW_SIZE, "bp (", WINDOW_SIZE/1000, "kb )\n")
cat("Minimum genes for enrichment:", MIN_GENES, "\n")
cat("GO FDR threshold:", GO_FDR_THRESHOLD, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("=================================================================\n\n")

################################################################################
# PART 2: LOAD AND PREPARE GENOME ANNOTATION FILES
################################################################################

cat("Step 1: Loading genome annotation files...\n")

# Load GFF3 file (gene locations)
cat("  Reading GFF3 file...\n")

# Read GFF3 - the key is to let fread handle it without fill=TRUE initially
# Then filter comment lines after reading
gff <- fread("Tintermedium_503_v2.1.gene.gff3", 
             sep = "\t",
             header = FALSE)

cat("  Initial read: ", ncol(gff), " columns, ", nrow(gff), " rows\n", sep = "")

# Filter out comment lines (lines starting with #)
gff <- gff[!grepl("^#", V1)]

cat("  After removing comments: ", nrow(gff), " rows\n", sep = "")

# Check we have 9 columns (standard GFF3 format)
if (ncol(gff) != 9) {
  cat("  ERROR: Expected 9 columns but got", ncol(gff), "\n")
  cat("  First few rows:\n")
  print(head(gff))
  stop("GFF3 file format issue")
}

# Set proper column names
colnames(gff) <- c("seqid", "source", "type", "start", "end", 
                   "score", "strand", "phase", "attributes")

# Filter for gene features only
genes <- gff[type == "gene"]
cat("  Found", nrow(genes), "genes in GFF3\n")

# Extract gene IDs from attributes column
genes[, gene_id := str_extract(attributes, "ID=[^;]+")]
genes[, gene_id := gsub("ID=", "", gene_id)]

# Remove version suffix to match annotation file format
# GFF3 has: Thint.01G0000100.v2.1
# Annotation has: Thint.01G0000100
genes[, gene_id := gsub("\\.v[0-9\\.]+$", "", gene_id)]

# Clean up chromosome names - they should already be in Chr## format in GFF3
genes[, chr := seqid]  # Keep original chromosome names from GFF3

# Select relevant columns
genes <- genes[, .(gene_id, chr, start, end, strand)]

cat("  Processed", nrow(genes), "gene annotations\n")
cat("  Sample gene IDs (after version removal):", paste(head(genes$gene_id, 3), collapse = ", "), "\n")
cat("  Sample gene chromosomes:", paste(head(unique(genes$chr), 5), collapse = ", "), "\n")

# Load functional annotation file (GO terms)
cat("  Reading annotation info file...\n")

# This file may be .txt or .txt.gz
annotation_file <- "Tintermedium_503_v2.1.P14.annotation_info.txt"
if (!file.exists(annotation_file)) {
  annotation_file <- "Tintermedium_503_v2.1.P14.annotation_info.txt.gz"
}

# Read annotation file
# Note: This file has tab-separated columns with GO terms
annotations <- fread(annotation_file, sep = "\t", fill = TRUE, quote = "")

cat("  Loaded annotations for", nrow(annotations), "entries\n")

# Display first few column names to understand structure
cat("  Annotation file columns:\n")
print(colnames(annotations)[1:min(10, ncol(annotations))])

# Based on the T. intermedium annotation file structure:
# Column 2: locusName (gene IDs like Tintermedium_503_v2.1.01G000100)
# Column 10: GO (space-separated GO terms like "GO:0004568 GO:0006032")

gene_col <- "locusName"
go_col <- "GO"

cat("  Using column '", gene_col, "' as gene identifier\n", sep = "")
cat("  Using column '", go_col, "' for GO terms\n", sep = "")

################################################################################
# PART 3: LOAD ALL GWAS RESULTS
################################################################################

cat("\nStep 2: Loading GWAS results for all 19 bioclimatic variables...\n")

# List of all bioclimatic variables
bio_vars <- paste0("bio", 1:19)

# Read all GWAS result files
gwas_results_list <- list()

for (bio in bio_vars) {
  file_path <- file.path("Analysis_Outputs_and_Visualizations", "Extras",
                        paste0("gwas_results_", bio, ".csv"))
  
  if (file.exists(file_path)) {
    gwas_results_list[[bio]] <- fread(file_path)
    n_sig <- sum(gwas_results_list[[bio]]$significant, na.rm = TRUE)
    cat("  ", bio, ": ", nrow(gwas_results_list[[bio]]), " markers, ", 
        n_sig, " significant\n", sep = "")
  } else {
    cat("  WARNING: File not found:", file_path, "\n")
  }
}

# Combine all results
all_gwas <- rbindlist(gwas_results_list, use.names = TRUE, fill = TRUE)
cat("\nTotal GWAS results loaded:", nrow(all_gwas), "marker-trait combinations\n")

################################################################################
# PART 4: PARSE SNP POSITIONS AND FILTER SIGNIFICANT SNPs
################################################################################

cat("\nStep 3: Parsing SNP positions and filtering significant markers...\n")

# Parse chromosome and position from marker names
# Format: SJ07_677959512 or SV03_79875884
all_gwas[, chr := str_extract(marker, "^[^_]+")]
all_gwas[, chr := gsub("^S[JV]", "", chr)]  # Remove SJ or SV prefix
all_gwas[, pos := as.numeric(str_extract(marker, "[0-9]+$"))]

# Standardize chromosome names to match GFF3 format
# GFF3 has "Chr01", "Chr02", etc.
# GWAS has "01", "02", etc. after removing prefix
all_gwas[, chr_gff := paste0("Chr", chr)]

cat("  Sample chromosome conversions:\n")
cat("    GWAS markers -> GFF3 format:\n")
sample_chrs <- head(unique(all_gwas[, .(marker, chr, chr_gff)]), 5)
for (i in 1:nrow(sample_chrs)) {
  cat("      ", sample_chrs$marker[i], " -> chr:", sample_chrs$chr_gff[i], "\n", sep = "")
}

# Filter for significant SNPs based on p-value threshold
sig_snps <- all_gwas[p_value < PVAL_THRESHOLD]

cat("  Total significant SNPs (p <", PVAL_THRESHOLD, "):", nrow(sig_snps), "\n")
cat("  Breakdown by trait:\n")
sig_summary <- sig_snps[, .(n_snps = .N, 
                            min_pval = min(p_value),
                            mean_effect = mean(abs(effect))), 
                        by = trait]
print(sig_summary)

# Save significant SNPs summary
fwrite(sig_snps, file.path(OUTPUT_DIR, "significant_SNPs_all_traits.csv"))
fwrite(sig_summary, file.path(OUTPUT_DIR, "significant_SNPs_summary.csv"))

################################################################################
# PART 5: MAP SNPs TO GENES
################################################################################

cat("\nStep 4: Mapping significant SNPs to genes...\n")
cat("  Using window size:", WINDOW_SIZE, "bp\n")

# Function to find genes near a SNP
map_snp_to_genes <- function(snp_chr, snp_pos, gene_data, window = WINDOW_SIZE) {
  # Find genes on the same chromosome within the window
  nearby_genes <- gene_data[chr == snp_chr & 
                           ((start - window) <= snp_pos) & 
                           ((end + window) >= snp_pos)]
  
  if (nrow(nearby_genes) > 0) {
    nearby_genes[, distance := pmin(abs(start - snp_pos), abs(end - snp_pos))]
    return(nearby_genes$gene_id)
  }
  return(character(0))
}

# Map each significant SNP to nearby genes
sig_snps[, genes := lapply(1:.N, function(i) {
  map_snp_to_genes(chr_gff[i], pos[i], genes)
})]

# Expand to one row per SNP-gene pair
snp_gene_map <- sig_snps[, .(gene_id = unlist(genes)), 
                         by = .(marker, chr, chr_gff, pos, p_value, effect, trait)]

cat("  Mapped", nrow(snp_gene_map), "SNP-gene associations\n")
cat("  Unique genes identified:", length(unique(snp_gene_map$gene_id)), "\n")

# Save SNP-gene mapping
fwrite(snp_gene_map, file.path(OUTPUT_DIR, "SNP_to_gene_mapping.csv"))

# Get gene lists for each trait
trait_genes <- snp_gene_map[, .(genes = list(unique(gene_id)), 
                                n_genes = length(unique(gene_id))), 
                            by = trait]

cat("\nGenes per trait:\n")
print(trait_genes[, .(trait, n_genes)])

################################################################################
# PART 6: PREPARE GO ANNOTATIONS
################################################################################

cat("\nStep 5: Preparing GO term annotations...\n")

# The T. intermedium annotation file has GO terms in column 10 "GO"
# Format: space-separated GO terms like "GO:0004568 GO:0006032 GO:0016998"

# Extract gene to GO mapping
gene2GO <- annotations[, .(gene = get(gene_col), GO = get(go_col))]

# Remove rows with missing or empty GO terms
gene2GO <- gene2GO[!is.na(GO) & GO != "" & GO != "-"]

cat("  Found", nrow(gene2GO), "genes with GO annotations\n")

# Split space-separated GO terms into individual terms
gene2GO[, GO_list := strsplit(GO, " +")]  # Split on one or more spaces
gene2GO_expanded <- gene2GO[, .(GO_term = unlist(GO_list)), by = gene]

# Clean up GO terms (remove any extra whitespace)
gene2GO_expanded[, GO_term := trimws(GO_term)]

# Keep only valid GO terms (format: GO:XXXXXXX)
gene2GO_expanded <- gene2GO_expanded[grep("^GO:[0-9]{7}$", GO_term)]

cat("  Total gene-GO associations:", nrow(gene2GO_expanded), "\n")
cat("  Unique genes with GO terms:", length(unique(gene2GO_expanded$gene)), "\n")
cat("  Unique GO terms:", length(unique(gene2GO_expanded$GO_term)), "\n")

# Convert to named list format required by topGO
# gene2GO_list <- list(
#   "gene1" = c("GO:0001", "GO:0002"),
#   "gene2" = c("GO:0003"),
#   ...
# )
gene2GO_list <- split(gene2GO_expanded$GO_term, gene2GO_expanded$gene)

cat("  Created gene2GO list for", length(gene2GO_list), "genes\n")

# Show some example genes with their GO terms
cat("\n  Example gene-GO mappings:\n")
example_genes <- head(names(gene2GO_list), 3)
for (gene in example_genes) {
  cat("    ", gene, ": ", length(gene2GO_list[[gene]]), " GO terms (", 
      paste(head(gene2GO_list[[gene]], 3), collapse = ", "), "...)\n", sep = "")
}

# Verify that our GWAS genes can be found in the GO annotation
if (exists("snp_gene_map") && nrow(snp_gene_map) > 0) {
  gwas_genes <- unique(snp_gene_map$gene_id)
  genes_with_go <- sum(gwas_genes %in% names(gene2GO_list))
  cat("\n  GWAS genes with GO annotations:", genes_with_go, "/", length(gwas_genes), 
      "(", round(100*genes_with_go/length(gwas_genes), 1), "%)\n", sep = "")
  
  if (genes_with_go == 0) {
    cat("  WARNING: None of the GWAS genes have GO annotations!\n")
    cat("  This suggests a gene ID mismatch problem.\n")
    cat("  Sample GWAS gene IDs:", paste(head(gwas_genes, 3), collapse = ", "), "\n")
    cat("  Sample GO gene IDs:", paste(head(names(gene2GO_list), 3), collapse = ", "), "\n")
  }
}

################################################################################
# PART 7: RUN GO ENRICHMENT ANALYSIS FOR EACH TRAIT
################################################################################

cat("\nStep 6: Running GO enrichment analysis...\n")

# Function to run GO enrichment for a gene list
run_GO_enrichment <- function(gene_list, background_genes, gene2GO, 
                              trait_name, ont = "BP") {
  
  cat("\n  Analyzing trait:", trait_name, "\n")
  
  # Check if we have enough genes
  if (length(gene_list) < MIN_GENES) {
    cat("    WARNING: Only", length(gene_list), "genes - need at least", 
        MIN_GENES, "for enrichment. Skipping.\n")
    return(NULL)
  }
  
  # Check how many genes have GO annotations
  genes_with_GO <- sum(gene_list %in% names(gene2GO))
  cat("    Genes with GO annotations:", genes_with_GO, "/", length(gene_list), "\n")
  
  if (genes_with_GO < MIN_GENES) {
    cat("    WARNING: Not enough genes with GO annotations. Skipping.\n")
    return(NULL)
  }
  
  # Create gene list factor (1 = interesting, 0 = background)
  gene_universe <- factor(as.integer(background_genes %in% gene_list))
  names(gene_universe) <- background_genes
  
  # Create topGO object
  GOdata <- new("topGOdata",
                ontology = ont,
                allGenes = gene_universe,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO,
                nodeSize = 5)  # Minimum 5 genes per GO term
  
  # Run Fisher's exact test
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Get results table
  results_table <- GenTable(GOdata, 
                           Fisher = resultFisher,
                           orderBy = "Fisher",
                           topNodes = 50)
  
  # Add FDR correction
  results_table$FDR <- p.adjust(results_table$Fisher, method = "BH")
  
  # Add trait name
  results_table$trait <- trait_name
  results_table$ontology <- ont
  
  # Filter for significant results
  sig_results <- results_table[results_table$FDR < GO_FDR_THRESHOLD, ]
  
  cat("    Significant GO terms (FDR <", GO_FDR_THRESHOLD, "):", 
      nrow(sig_results), "\n")
  
  return(results_table)
}

# Run enrichment for each trait with enough genes
if (length(gene2GO_list) > 0) {
  
  # Get background (all genes with GO annotations)
  background_genes <- names(gene2GO_list)
  
  # Filter trait_genes to only include traits with genes
  trait_genes <- trait_genes[n_genes >= MIN_GENES]
  
  if (nrow(trait_genes) == 0) {
    cat("\n  WARNING: No traits have enough genes (>= ", MIN_GENES, ") for enrichment analysis\n", sep = "")
    cat("  This usually means:\n")
    cat("    - P-value threshold is too strict (try increasing PVAL_THRESHOLD)\n")
    cat("    - Window size is too small (try increasing WINDOW_SIZE)\n")
    cat("    - Chromosome names don't match between GWAS and GFF3\n")
    cat("\n  Skipping GO enrichment analysis.\n")
  } else {
    cat("  Traits with sufficient genes for analysis:", nrow(trait_genes), "\n")
    
    # Run for each trait
    all_GO_results <- list()
    
    for (i in 1:nrow(trait_genes)) {
      trait_name <- trait_genes$trait[i]
      gene_list <- trait_genes$genes[[i]]
      
      # Run for Biological Process
      go_result <- tryCatch({
        run_GO_enrichment(gene_list, background_genes, gene2GO_list, 
                         trait_name, ont = "BP")
      }, error = function(e) {
        cat("    ERROR:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(go_result)) {
        all_GO_results[[trait_name]] <- go_result
      }
    }
    
    # Combine all results
    if (length(all_GO_results) > 0) {
      combined_GO_results <- rbindlist(all_GO_results, use.names = TRUE, fill = TRUE)
      
      # Save results
      fwrite(combined_GO_results, 
             file.path(OUTPUT_DIR, "GO_enrichment_results_all_traits.csv"))
      
      # Save significant results only
      sig_GO <- combined_GO_results[FDR < GO_FDR_THRESHOLD]
      fwrite(sig_GO, 
             file.path(OUTPUT_DIR, "GO_enrichment_significant.csv"))
      
      cat("\n=================================================================\n")
      cat("GO Enrichment Summary:\n")
      cat("=================================================================\n")
      cat("Total GO terms tested:", nrow(combined_GO_results), "\n")
      cat("Significant GO terms (FDR <", GO_FDR_THRESHOLD, "):", nrow(sig_GO), "\n")
      cat("\nResults saved to:", OUTPUT_DIR, "\n")
      cat("=================================================================\n")
    }
  }  # End of else block
  
} else {
  cat("\nWARNING: No GO annotations found. Cannot run enrichment analysis.\n")
  cat("Please check the annotation file structure.\n")
}

################################################################################
# PART 8: CREATE VISUALIZATIONS
################################################################################

cat("\nStep 7: Creating visualizations...\n")

if (exists("combined_GO_results") && nrow(combined_GO_results) > 0) {
  
  # 1. Bar plot of top GO terms per trait
  top_GO_per_trait <- combined_GO_results %>%
    group_by(trait) %>%
    arrange(FDR) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(neg_log10_FDR = -log10(as.numeric(FDR)))
  
  if (nrow(top_GO_per_trait) > 0) {
    p1 <- ggplot(top_GO_per_trait, 
                 aes(x = reorder(Term, neg_log10_FDR), y = neg_log10_FDR, fill = trait)) +
      geom_col() +
      coord_flip() +
      facet_wrap(~trait, scales = "free_y") +
      theme_minimal() +
      labs(title = "Top GO Terms by Bioclimatic Variable",
           x = "GO Term",
           y = "-log10(FDR)",
           fill = "Trait") +
      theme(legend.position = "none",
            axis.text.y = element_text(size = 8))
    
    ggsave(file.path(OUTPUT_DIR, "top_GO_terms_by_trait.png"), 
           p1, width = 12, height = 10, dpi = 300)
    cat("  Saved: top_GO_terms_by_trait.png\n")
  }
  
  # 2. Dot plot of significant GO terms
  sig_GO_plot <- combined_GO_results[FDR < GO_FDR_THRESHOLD]
  
  if (nrow(sig_GO_plot) > 0) {
    sig_GO_plot$Significant <- as.numeric(sig_GO_plot$Significant)
    sig_GO_plot$neg_log10_FDR <- -log10(as.numeric(sig_GO_plot$FDR))
    
    p2 <- ggplot(sig_GO_plot, 
                 aes(x = trait, y = Term, size = Significant, color = neg_log10_FDR)) +
      geom_point() +
      scale_color_gradient(low = "blue", high = "red") +
      theme_minimal() +
      labs(title = "Significant GO Terms Across Traits",
           x = "Bioclimatic Variable",
           y = "GO Term",
           size = "# Genes",
           color = "-log10(FDR)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8))
    
    ggsave(file.path(OUTPUT_DIR, "significant_GO_dotplot.png"), 
           p2, width = 12, height = 8, dpi = 300)
    cat("  Saved: significant_GO_dotplot.png\n")
  }
}

# 3. Create summary statistics plot
if (nrow(sig_snps) > 0) {
  p3 <- ggplot(sig_summary, aes(x = trait, y = n_snps, fill = trait)) +
    geom_col() +
    theme_minimal() +
    labs(title = "Number of Significant SNPs per Bioclimatic Variable",
         x = "Bioclimatic Variable",
         y = "Number of Significant SNPs",
         caption = paste("P-value threshold:", PVAL_THRESHOLD)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(OUTPUT_DIR, "significant_SNPs_by_trait.png"), 
         p3, width = 10, height = 6, dpi = 300)
  cat("  Saved: significant_SNPs_by_trait.png\n")
}

################################################################################
# PART 9: GROUP ANALYSIS (Temperature vs Precipitation traits)
################################################################################

cat("\nStep 8: Running grouped trait analysis...\n")

# Group traits by type
temp_traits <- c("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_3",
                "wc2.1_30s_bio_4", "wc2.1_30s_bio_5", "wc2.1_30s_bio_6",
                "wc2.1_30s_bio_7", "wc2.1_30s_bio_8", "wc2.1_30s_bio_9",
                "wc2.1_30s_bio_10", "wc2.1_30s_bio_11")

precip_traits <- c("wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14",
                  "wc2.1_30s_bio_15", "wc2.1_30s_bio_16", "wc2.1_30s_bio_17",
                  "wc2.1_30s_bio_18", "wc2.1_30s_bio_19")

# Get genes for each group
temp_genes <- unique(snp_gene_map[trait %in% temp_traits, gene_id])
precip_genes <- unique(snp_gene_map[trait %in% precip_traits, gene_id])

cat("  Temperature-related traits: ", length(temp_genes), " genes\n")
cat("  Precipitation-related traits: ", length(precip_genes), " genes\n")

# Run enrichment for groups
if (length(gene2GO_list) > 0) {
  
  # Temperature group
  if (length(temp_genes) >= MIN_GENES) {
    temp_GO <- run_GO_enrichment(temp_genes, background_genes, gene2GO_list,
                                "Temperature_Group", ont = "BP")
    if (!is.null(temp_GO)) {
      fwrite(temp_GO, file.path(OUTPUT_DIR, "GO_enrichment_temperature_group.csv"))
    }
  }
  
  # Precipitation group
  if (length(precip_genes) >= MIN_GENES) {
    precip_GO <- run_GO_enrichment(precip_genes, background_genes, gene2GO_list,
                                  "Precipitation_Group", ont = "BP")
    if (!is.null(precip_GO)) {
      fwrite(precip_GO, file.path(OUTPUT_DIR, "GO_enrichment_precipitation_group.csv"))
    }
  }
}

################################################################################
# PART 10: FINAL SUMMARY
################################################################################

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=================================================================\n")
cat("\nFiles created in", OUTPUT_DIR, ":\n")
cat("  - significant_SNPs_all_traits.csv\n")
cat("  - significant_SNPs_summary.csv\n")
cat("  - SNP_to_gene_mapping.csv\n")
cat("  - GO_enrichment_results_all_traits.csv\n")
cat("  - GO_enrichment_significant.csv\n")
cat("  - Various visualization plots (.png)\n")
cat("\nNext steps:\n")
cat("  1. Review the significant GO terms in the output files\n")
cat("  2. Look for biological processes related to climate adaptation\n")
cat("  3. Identify candidate genes for follow-up studies\n")
cat("  4. Consider adjusting parameters and re-running:\n")
cat("     - Try different p-value thresholds (", PVAL_THRESHOLD, ")\n")
cat("     - Try different window sizes (", WINDOW_SIZE, " bp)\n")
cat("=================================================================\n")
