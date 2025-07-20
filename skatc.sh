#!/bin/bash
#SBATCH --job-name=skat_c_analysis
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=skat_c_%j.out
#SBATCH --error=skat_c_%j.err

# Load R module (adjust version if needed)
module load R/4.4.2-gfbf-2024a

# Run R script
R --vanilla << 'EOF'

# Load required libraries
library(SKAT)

# Set working directory to the data folder
setwd("/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/05_SKAT/WRMT_Attack_Standard_t")

cat("=== SKAT-C ANALYSIS ===\n")
cat("Working directory:", getwd(), "\n")
cat("Starting analysis at:", Sys.time(), "\n\n")

# Define file paths
plink_prefix <- "plink_maf_filtered"
bed_file <- paste0(plink_prefix, ".bed")
bim_file <- paste0(plink_prefix, ".bim")
fam_file <- paste0(plink_prefix, ".fam")
cov_file <- "covariates.txt"
setid_file <- "regions_cleaned.txt"

# Define output file names
ssd_file <- paste0(plink_prefix, ".SSD")
ssd_info_file <- paste0(plink_prefix, ".SSD.info")

# Step 1: Check if all required files exist
cat("Step 1: Checking input files...\n")
required_files <- c(bed_file, bim_file, fam_file, cov_file, setid_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    cat("ERROR: Missing required files:\n")
    for (file in missing_files) {
        cat("  -", file, "\n")
    }
    stop("Please ensure all required files are present.")
}

cat("All required input files found:\n")
for (file in required_files) {
    size_mb <- round(file.info(file)$size / 1024 / 1024, 2)
    cat(sprintf("  - %s (%.2f MB)\n", file, size_mb))
}

# Step 2: Validate the SetID file format
cat("\nStep 2: Validating SetID file format...\n")

# Check if the file has a header and remove it if necessary
first_line <- readLines(setid_file, n = 1)
if (grepl("^Gene_name|^gene|^Gene", first_line, ignore.case = TRUE)) {
    cat("Header detected in SetID file. Removing header...\n")
    # Read all lines except the first
    all_lines <- readLines(setid_file)
    writeLines(all_lines[-1], setid_file)
}

# Read and validate the SetID file
setid_data <- read.table(setid_file, header = FALSE, stringsAsFactors = FALSE)
colnames(setid_data) <- c("Gene_name", "SNP_ID")

cat("SetID file summary:\n")
cat(sprintf("  - Total SNP assignments: %d\n", nrow(setid_data)))
cat(sprintf("  - Number of genes: %d\n", length(unique(setid_data$Gene_name))))
cat(sprintf("  - SNPs per gene (median): %.1f\n", median(table(setid_data$Gene_name))))

# Show top genes by SNP count
gene_counts <- table(setid_data$Gene_name)
top_genes <- head(sort(gene_counts, decreasing = TRUE), 5)
cat("Top 5 genes by SNP count:\n")
for (i in 1:length(top_genes)) {
    cat(sprintf("  - %s: %d SNPs\n", names(top_genes)[i], top_genes[i]))
}

# Step 3: Generate SSD files if they don't exist
cat("\nStep 3: Generating SSD files...\n")

if (!file.exists(ssd_file) || !file.exists(ssd_info_file)) {
    cat("Creating SSD files from PLINK data...\n")
    
    Generate_SSD_SetID(
        File.Bed = bed_file,
        File.Bim = bim_file, 
        File.Fam = fam_file,
        File.SetID = setid_file,
        File.SSD = ssd_file,
        File.Info = ssd_info_file
    )
    
    cat("SSD files generated successfully.\n")
} else {
    cat("SSD files already exist, skipping generation.\n")
}

# Step 4: Read phenotype and covariate data
cat("\nStep 4: Reading phenotype and covariate data...\n")

# Read the covariate file to check format
cov_sample <- readLines(cov_file, n = 3)
cat("First few lines of covariate file:\n")
for (i in 1:length(cov_sample)) {
    cat(sprintf("  Line %d: %s\n", i, substr(cov_sample[i], 1, 100)))
}

# Auto-detect if covariate file has header
first_line_fields <- strsplit(cov_sample[1], "\\s+")[[1]]
# Check if any of the first fields (after FID/IID) are non-numeric
has_header <- FALSE
if (length(first_line_fields) > 2) {
    test_fields <- first_line_fields[3:min(5, length(first_line_fields))]
    has_header <- any(is.na(suppressWarnings(as.numeric(test_fields))))
}

cat(sprintf("Covariate file header detected: %s\n", has_header))

# Read FAM and covariate data
FAM_Cov <- Read_Plink_FAM_Cov(
    Filename = fam_file,
    File_Cov = cov_file,
    Is.binary = FALSE,  # Continuous phenotype
    cov_header = has_header
)

cat(sprintf("Data loaded: %d samples, %d variables\n", nrow(FAM_Cov), ncol(FAM_Cov)))
cat("Column names:", paste(colnames(FAM_Cov), collapse = ", "), "\n")

# Step 5: Set up the null model
cat("\nStep 5: Setting up null model...\n")

# Standard PLINK FAM file has phenotype in column 6 (named "Phenotype" by Read_Plink_FAM_Cov)
phenotype_col <- "Phenotype"
covariate_cols <- colnames(FAM_Cov)[7:ncol(FAM_Cov)]  # Covariates start from column 7

cat(sprintf("Using phenotype column: %s\n", phenotype_col))
cat(sprintf("Using %d covariate columns: %s\n", 
           length(covariate_cols), paste(covariate_cols, collapse = ", ")))

# Check phenotype distribution
y <- FAM_Cov[, phenotype_col]
cat("Phenotype summary:\n")
cat(sprintf("  - N samples: %d\n", sum(!is.na(y))))
cat(sprintf("  - Mean: %.3f\n", mean(y, na.rm = TRUE)))
cat(sprintf("  - SD: %.3f\n", sd(y, na.rm = TRUE)))
cat(sprintf("  - Range: %.3f to %.3f\n", min(y, na.rm = TRUE), max(y, na.rm = TRUE)))

# Create formula for null model
if (length(covariate_cols) > 0) {
    formula_str <- paste(phenotype_col, "~", paste(covariate_cols, collapse = " + "))
} else {
    formula_str <- paste(phenotype_col, "~ 1")  # Intercept only model
}

cat(sprintf("Null model formula: %s\n", formula_str))

# Fit null model for continuous phenotype
cat("Fitting null model...\n")
obj <- SKAT_Null_Model(
    formula = as.formula(formula_str), 
    data = FAM_Cov,
    out_type = "C"  # Continuous phenotype
)

cat("Null model fitted successfully.\n")

# Step 6: Open SSD file and perform SKAT-C analysis
cat("\nStep 6: Performing SKAT-C analysis...\n")

# Open SSD file
SSD.INFO <- Open_SSD(File.SSD = ssd_file, File.Info = ssd_info_file)

cat(sprintf("SSD file opened: %d gene sets available\n", SSD.INFO$nSets))

# Run SKAT-C analysis (Combined sum test for common and rare variants)
cat("Running SKAT-C analysis across all gene sets...\n")

results <- SKAT_CommonRare.SSD.All(
    SSD.INFO = SSD.INFO,
    obj = obj,
    method = "C"  # Combined sum test (SKAT-C)
)

cat("SKAT-C analysis completed.\n")

# Step 7: Process and save results
cat("\nStep 7: Processing results...\n")

# Display results summary
cat("=== SKAT-C ANALYSIS RESULTS ===\n")
cat(sprintf("Total gene sets analyzed: %d\n", nrow(results$results)))

# Count significant results at different thresholds
sig_001 <- sum(results$results$P.value < 0.001, na.rm = TRUE)
sig_01 <- sum(results$results$P.value < 0.01, na.rm = TRUE)
sig_05 <- sum(results$results$P.value < 0.05, na.rm = TRUE)

cat(sprintf("Significant results (p < 0.001): %d\n", sig_001))
cat(sprintf("Significant results (p < 0.01): %d\n", sig_01))
cat(sprintf("Significant results (p < 0.05): %d\n", sig_05))

if (nrow(results$results) > 0) {
    min_p <- min(results$results$P.value, na.rm = TRUE)
    most_sig_idx <- which.min(results$results$P.value)
    most_sig_gene <- results$results$SetID[most_sig_idx]
    
    cat(sprintf("Minimum p-value: %.2e\n", min_p))
    cat(sprintf("Most significant gene: %s (p = %.2e)\n", most_sig_gene, min_p))
}

# Sort results by p-value
results_sorted <- results$results[order(results$results$P.value), ]

# Display top 10 results
cat("\nTop 10 most significant results:\n")
top10 <- head(results_sorted, 10)
for (i in 1:nrow(top10)) {
    cat(sprintf("%2d. %s: p = %.2e (N.Rare = %d, N.Common = %d)\n", 
               i, top10$SetID[i], top10$P.value[i], 
               top10$N.Marker.Rare[i], top10$N.Marker.Common[i]))
}

# Save detailed results
output_file <- "skat_c_results.txt"
write.table(results_sorted, output_file, 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

cat(sprintf("\nDetailed results saved to: %s\n", output_file))

# Save summary results for significant genes
if (sig_05 > 0) {
    sig_results <- results_sorted[results_sorted$P.value < 0.05, ]
    sig_output_file <- "skat_c_significant_results.txt"
    write.table(sig_results, sig_output_file, 
                row.names = FALSE, col.names = TRUE, 
                quote = FALSE, sep = "\t")
    cat(sprintf("Significant results (p < 0.05) saved to: %s\n", sig_output_file))
}

# Create a summary statistics file
summary_stats <- data.frame(
    Analysis = "SKAT-C",
    Total_Genes = nrow(results$results),
    Significant_p001 = sig_001,
    Significant_p01 = sig_01,
    Significant_p05 = sig_05,
    Min_Pvalue = min_p,
    Most_Significant_Gene = most_sig_gene,
    Analysis_Date = Sys.time()
)

write.table(summary_stats, "skat_c_summary.txt", 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

cat("Summary statistics saved to: skat_c_summary.txt\n")

# Close SSD file
Close_SSD()

cat("\n=== ANALYSIS COMPLETED SUCCESSFULLY ===\n")
cat("Analysis finished at:", as.character(Sys.time()), "\n")
cat("Check the following output files:\n")
cat("  - skat_c_results.txt (all results)\n")
if (sig_05 > 0) {
    cat("  - skat_c_significant_results.txt (significant results)\n")
}
cat("  - skat_c_summary.txt (summary statistics)\n")

EOF
