# GO Enrichment Analysis - User Guide
## For Intermediate Wheatgrass Climate Adaptation GWAS

### What This Script Does

This R script takes your GWAS results and:
1. **Identifies significant SNPs** based on p-value threshold
2. **Maps SNPs to nearby genes** using the GFF3 file
3. **Extracts GO terms** for those genes
4. **Tests for enrichment** - finds which biological processes are over-represented
5. **Creates visualizations** and summary reports

---

## Before You Run

### 1. Check Your File Structure

Make sure these files are in the right places:

```
C:\Users\Lizzy\OneDrive\Documents\UH_Manoa_Classes\IWG_Chapter_2\
├── Tintermedium_503_v2.1.gene.gff3
├── Tintermedium_503_v2.1.P14.annotation_info.txt  (or .txt.gz)
├── GO_enrichment_IWG_GWAS.R  (this script)
└── Analysis_Outputs_and_Visualizations\
    ├── gwas_results_bio1.csv
    ├── gwas_results_bio2.csv
    ├── ... (all 19 files)
    └── gwas_results_bio19.csv
```

### 2. Install Required R Packages

Open R or RStudio and run:

```r
# Install CRAN packages
install.packages(c("data.table", "dplyr", "ggplot2", "stringr"))

# Install Bioconductor package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("topGO")
```

This only needs to be done once!

---

## Running the Script

### Option 1: Run the Entire Script (Recommended First Time)

1. Open RStudio
2. File → Open File → Select `GO_enrichment_IWG_GWAS.R`
3. Click "Source" button (or Ctrl+Shift+S)
4. Watch the console for progress messages

### Option 2: Run Section by Section

1. Open the script in RStudio
2. Run each section (PART 1, PART 2, etc.) separately using Ctrl+Enter
3. This lets you check results at each step

---

## Key Parameters You Can Adjust

At the top of the script (PART 1), you'll see these settings:

### 1. **P-value Threshold**
```r
PVAL_THRESHOLD <- 0.001
```
- **Lower = stricter** (fewer genes, higher confidence)
- **Higher = more permissive** (more genes, may include false positives)
- Try: 0.001, 0.0001, or 0.00001

### 2. **Window Size**
```r
WINDOW_SIZE <- 25000
```
- How far from a SNP to look for genes (in base pairs)
- **Smaller window** = only very nearby genes (conservative)
- **Larger window** = genes further away (captures more regulatory effects)
- Try: 10000 (10 kb), 25000 (25 kb), 50000 (50 kb), 100000 (100 kb)

### 3. **Minimum Genes**
```r
MIN_GENES <- 5
```
- Minimum number of genes needed to run enrichment test
- Too few genes = unreliable statistics
- Usually keep at 5-10

### 4. **GO FDR Threshold**
```r
GO_FDR_THRESHOLD <- 0.05
```
- False Discovery Rate cutoff for significant GO terms
- Standard is 0.05 (5% FDR)
- More stringent: 0.01 (1% FDR)

---

## Understanding the Output

The script creates a folder called `GO_Analysis_Results` with these files:

### CSV Files:

1. **significant_SNPs_all_traits.csv**
   - All SNPs that pass your p-value threshold
   - Columns: marker, chr, pos, p_value, effect, trait

2. **significant_SNPs_summary.csv**
   - Summary statistics per trait
   - How many significant SNPs, minimum p-value, mean effect size

3. **SNP_to_gene_mapping.csv**
   - Which genes are near which SNPs
   - One row per SNP-gene pair

4. **GO_enrichment_results_all_traits.csv**
   - All GO enrichment results
   - Includes non-significant terms

5. **GO_enrichment_significant.csv**
   - **⭐ MAIN RESULTS ⭐**
   - Only GO terms with FDR < 0.05
   - These are the biological processes enriched in your data!

6. **GO_enrichment_temperature_group.csv**
   - GO enrichment for temperature-related traits combined

7. **GO_enrichment_precipitation_group.csv**
   - GO enrichment for precipitation-related traits combined

### Visualization Files (PNG):

1. **significant_SNPs_by_trait.png**
   - Bar chart showing how many significant SNPs per trait

2. **top_GO_terms_by_trait.png**
   - Top 10 GO terms for each bioclimatic variable

3. **significant_GO_dotplot.png**
   - Bubble plot showing enriched GO terms across traits
   - Size = number of genes
   - Color = significance level

---

## Interpreting GO Enrichment Results

### What to Look For:

The **GO_enrichment_significant.csv** file has these key columns:

- **GO.ID**: The GO term identifier (e.g., GO:0009414)
- **Term**: Human-readable description (e.g., "response to water deprivation")
- **Significant**: Number of your genes with this term
- **Expected**: Number expected by chance
- **Fisher**: P-value from Fisher's exact test
- **FDR**: False Discovery Rate (adjusted p-value)
- **trait**: Which bioclimatic variable
- **ontology**: BP (Biological Process), MF (Molecular Function), or CC (Cellular Component)

### Example Result:

```
GO.ID      Term                           Significant  Expected  FDR      trait
GO:0009414 response to water deprivation  15          2.3       0.001    bio12
```

**Interpretation:** 
- You found 15 genes involved in water deprivation response
- By chance, you'd expect only 2.3
- This is highly significant (FDR = 0.001)
- This was for bio12 (annual precipitation)

**Biological meaning:** 
Your GWAS identified genetic variants in genes involved in drought response, which makes sense for precipitation adaptation!

---

## Common Issues and Solutions

### Issue 1: "No GO annotations found"

**Cause:** The annotation file structure doesn't match what the script expects

**Solution:** 
1. Check the annotation file manually:
   ```r
   library(data.table)
   ann <- fread("Tintermedium_503_v2.1.P14.annotation_info.txt", nrows = 10)
   colnames(ann)  # See what columns exist
   ```

2. Look for columns with GO terms
3. Adjust lines 87-120 in the script to match your file structure

### Issue 2: "Not enough genes for enrichment"

**Cause:** Your p-value threshold is too strict, or window size is too small

**Solutions:**
- Increase PVAL_THRESHOLD (e.g., from 0.0001 to 0.001)
- Increase WINDOW_SIZE (e.g., from 25000 to 50000)
- Check if you have significant SNPs at all (look at significant_SNPs_summary.csv)

### Issue 3: "No significant GO terms"

**Possible reasons:**
1. Your genes don't have strong functional clustering
2. Sample size is too small for some traits
3. Need to adjust GO_FDR_THRESHOLD to be less strict

**Try:**
- Look at the full results file (GO_enrichment_results_all_traits.csv)
- Sort by Fisher p-value to see top terms even if not FDR-significant
- Combine related traits (the script does this automatically for temp vs precip)

### Issue 4: Error reading GFF3 or annotation files

**Solution:**
```r
# Test reading files manually
gff <- fread("Tintermedium_503_v2.1.gene.gff3", skip = "#", nrows = 100)
head(gff)

# If the .txt file doesn't work, try the .gz file
ann <- fread("Tintermedium_503_v2.1.P14.annotation_info.txt.gz", nrows = 100)
head(ann)
```

---

## Tips for Better Results

### 1. Start with Lenient Parameters
- Use p < 0.001 first
- Use 50 kb window
- See what you get, then tighten if needed

### 2. Compare Different Thresholds
Run the script multiple times with different settings:
```r
# Run 1: Lenient
PVAL_THRESHOLD <- 0.001
WINDOW_SIZE <- 50000

# Run 2: Moderate  
PVAL_THRESHOLD <- 0.0001
WINDOW_SIZE <- 25000

# Run 3: Strict
PVAL_THRESHOLD <- 0.00001
WINDOW_SIZE <- 10000
```

Rename the output folders each time:
```r
OUTPUT_DIR <- "GO_Analysis_Results_lenient"
OUTPUT_DIR <- "GO_Analysis_Results_moderate"
OUTPUT_DIR <- "GO_Analysis_Results_strict"
```

### 3. Focus on Grouped Analyses
Temperature and precipitation traits are correlated, so combining them gives more statistical power!

### 4. Look for Consistent Patterns
GO terms that appear across multiple related traits are more likely to be real

---

## What's Next After GO Analysis?

### 1. Biological Interpretation
- Do the enriched processes make biological sense?
- How do they relate to climate adaptation?
- Which processes are shared vs. unique across traits?

### 2. Candidate Gene Identification
- Look at the SNP_to_gene_mapping.csv
- Find genes in enriched GO categories
- These are your top candidates for follow-up!

### 3. Literature Search
- Search for enriched GO terms + "climate adaptation" + "grasses"
- See if your findings match known biology

### 4. Visualization for Your Thesis/Paper
- The script creates publication-ready plots
- You can customize them further in R

### 5. Functional Validation (Future Work)
- Expression studies (RNA-seq)
- Phenotyping experiments
- CRISPR/gene editing (if feasible)

---

## Quick Troubleshooting Checklist

Before asking for help, check:

- [ ] All required files are in the right directories
- [ ] File paths in the script match your actual paths
- [ ] R packages are installed (run the install commands)
- [ ] The annotation file has GO terms (check manually)
- [ ] You have significant SNPs in your GWAS results
- [ ] Your p-value threshold isn't too strict (try 0.001 first)

---

## Getting Help

If you're stuck:

1. **Check the console output** - The script prints detailed progress messages
2. **Look at intermediate files** - Open the CSV files to see what's being produced
3. **Try running section by section** - Helps identify exactly where it fails
4. **Examine the data structures** - Use `head()`, `str()`, `colnames()` in R

---

## Expected Runtime

- **File loading**: 30 seconds
- **SNP mapping**: 1-2 minutes
- **GO enrichment per trait**: 30-60 seconds each
- **Total**: 5-15 minutes depending on your computer

If it's taking much longer, check if you're trying to test too many GO terms or have a huge gene list.

---

## Success Criteria

You'll know it worked if:

✅ Script runs without errors  
✅ GO_Analysis_Results folder is created  
✅ CSV files have data (not empty)  
✅ At least some GO terms have FDR < 0.05  
✅ Plots are generated and look reasonable  

Even if you don't get significant results for all traits, that's scientifically informative! Some traits may not have strong GO enrichment, which tells you something about the genetic architecture.

---

Good luck with your analysis! 🌾
