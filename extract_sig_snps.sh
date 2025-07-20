#!/bin/bash
#SBATCH --job-name=extract_sig_snps
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=extract_sig_snps_%j.out
#SBATCH --error=extract_sig_snps_%j.err

# Load R module
module load R/4.4.2-gfbf-2024a

# Set working directory
cd /mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/05_SKAT/CTPP_PA_Elision_Z_t

echo "=== EXTRACTING ACTUALLY TESTED SNPs FOR SIGNIFICANT GENES ==="
echo "Working directory: $(pwd)"
echo "Starting analysis at: $(date)"
echo

# Run R script to extract actually tested SNPs
R --vanilla << 'EOF'

# Load required libraries
library(SKAT)
library(data.table)

cat("=== EXTRACTING ACTUALLY TESTED SNPs FOR SIGNIFICANT GENES ===\n")
cat("Analysis started at:", as.character(Sys.time()), "\n\n")

# Define input files
results_file <- "skat_c_results.txt"
ssd_file <- "plink_maf_filtered.SSD"
ssd_info_file <- "plink_maf_filtered.SSD.info"
bim_file <- "plink_maf_filtered.bim"

# Define significance threshold
sig_threshold <- 5e-5  # p < 0.05 * 10^-4
maf_cutoff <- 0.05     # MAF cutoff for rare vs common

cat("Step 1: Checking input files...\n")

# Check if all required files exist
required_files <- c(results_file, ssd_file, ssd_info_file, bim_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    cat("ERROR: Missing required files:\n")
    for (file in missing_files) {
        cat("  -", file, "\n")
    }
    stop("Please ensure all required files are present.")
}

# Show file sizes
cat("Input files found:\n")
for (file in required_files) {
    size_mb <- round(file.info(file)$size / 1024 / 1024, 2)
    cat(sprintf("  - %s (%.2f MB)\n", file, size_mb))
}

cat("\nStep 2: Reading SKAT-C results...\n")

# Read SKAT-C results
skat_results <- fread(results_file)
cat(sprintf("Total genes tested: %d\n", nrow(skat_results)))
cat(sprintf("Significance threshold: p < %.0e\n", sig_threshold))

# Filter for significant genes
sig_genes <- skat_results[P.value < sig_threshold]
cat(sprintf("Significant genes found: %d\n", nrow(sig_genes)))

if (nrow(sig_genes) == 0) {
    cat("No significant genes found at the specified threshold.\n")
    cat("Most significant p-value:", min(skat_results$P.value, na.rm = TRUE), "\n")
    stop("Analysis stopped - no significant results to extract SNPs for.")
}

# Sort by p-value (most significant first)
sig_genes <- sig_genes[order(P.value)]

cat("\nSignificant genes:\n")
for (i in 1:nrow(sig_genes)) {
    cat(sprintf("%2d. %s: p = %.2e (tested %d/%d SNPs, %d rare, %d common)\n", 
               i, sig_genes$SetID[i], sig_genes$P.value[i], 
               sig_genes$N.Marker.Test[i], sig_genes$N.Marker.All[i],
               sig_genes$N.Marker.Rare[i], sig_genes$N.Marker.Common[i]))
}

cat("\nStep 3: Opening SSD file to extract actually tested SNPs...\n")

# Open SSD file
SSD.INFO <- Open_SSD(File.SSD = ssd_file, File.Info = ssd_info_file)
cat(sprintf("SSD file opened successfully with %d gene sets\n", SSD.INFO$nSets))

# Get information about sets in SSD file
set_info <- SSD.INFO$SetInfo
cat(sprintf("Available sets in SSD: %d\n", nrow(set_info)))

cat("\nStep 4: Reading BIM file for SNP details...\n")

# Read BIM file
bim_data <- fread(bim_file, header = FALSE)
colnames(bim_data) <- c("CHR", "SNP_ID", "cM", "POS", "A1", "A2")
cat(sprintf("SNPs in BIM file: %d\n", nrow(bim_data)))

# Create lookup table for SNP details
setkey(bim_data, SNP_ID)

cat("\nStep 5: Reading frequency file for MAF calculation...\n")

# Check if frequency file exists
freq_file <- "freq_data.frq"
if (file.exists(freq_file)) {
    freq_data <- fread(freq_file)
    setnames(freq_data, c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"))
    setkey(freq_data, SNP)
    cat(sprintf("Frequency data loaded: %d SNPs\n", nrow(freq_data)))
    use_real_maf <- TRUE
} else {
    cat("Frequency file not found. Using placeholder MAF values.\n")
    cat("To get real MAF: plink --bfile plink_maf_filtered --freq --out freq_data\n")
    use_real_maf <- FALSE
}

cat("\nStep 6: Extracting actually tested SNPs for each significant gene...\n")

# Initialize list to store results
all_snp_data <- list()

for (i in 1:nrow(sig_genes)) {
    gene_name <- sig_genes$SetID[i]
    gene_pvalue <- sig_genes$P.value[i]
    
    cat(sprintf("Processing %d/%d: %s...\n", i, nrow(sig_genes), gene_name))
    
    # Find the set index for this gene
    set_index <- which(set_info$SetID == gene_name)
    
    if (length(set_index) == 0) {
        cat(sprintf("  Warning: Gene %s not found in SSD file\n", gene_name))
        next
    }
    
    if (length(set_index) > 1) {
        cat(sprintf("  Warning: Multiple entries for %s, using first\n", gene_name))
        set_index <- set_index[1]
    }
    
    # Extract genotype matrix for this gene set (this gives us the actual tested SNPs)
    tryCatch({
        genotype_matrix <- Get_Genotypes_SSD(SSD.INFO, set_index, is_ID = TRUE)
        
        if (is.null(genotype_matrix) || ncol(genotype_matrix) == 0) {
            cat(sprintf("  Warning: No genotype data for %s\n", gene_name))
            next
        }
        
        # Get the SNP IDs that were actually tested
        tested_snp_ids <- colnames(genotype_matrix)
        n_tested <- length(tested_snp_ids)
        
        cat(sprintf("  Found %d actually tested SNPs\n", n_tested))
        
        # Verify this matches the SKAT results
        expected_tested <- sig_genes$N.Marker.Test[i]
        if (n_tested != expected_tested) {
            cat(sprintf("  Note: Expected %d SNPs, found %d SNPs\n", expected_tested, n_tested))
        }
        
        # Create data frame for this gene
        gene_snp_data <- data.table(
            Gene = gene_name,
            Gene_Pvalue = gene_pvalue,
            SNP_ID = tested_snp_ids
        )
        
        # Add to results list
        all_snp_data[[i]] <- gene_snp_data
        
    }, error = function(e) {
        cat(sprintf("  Error extracting SNPs for %s: %s\n", gene_name, e$message))
    })
}

# Close SSD file
Close_SSD()

cat("\nStep 7: Combining results and adding SNP details...\n")

# Combine all gene data
if (length(all_snp_data) == 0) {
    stop("No SNP data was successfully extracted!")
}

final_snp_data <- rbindlist(all_snp_data)
cat(sprintf("Total actually tested SNPs across significant genes: %d\n", nrow(final_snp_data)))

# Merge with BIM file to get SNP details
final_data <- merge(final_snp_data, bim_data, by = "SNP_ID", all.x = TRUE)

# Check for missing SNP details
missing_snp_details <- sum(is.na(final_data$CHR))
if (missing_snp_details > 0) {
    cat(sprintf("Warning: %d SNPs not found in BIM file\n", missing_snp_details))
}

# Remove SNPs not found in BIM file
final_data <- final_data[!is.na(CHR)]
cat(sprintf("SNPs with complete details: %d\n", nrow(final_data)))

cat("\nStep 8: Adding MAF information and rare/common classification...\n")

if (use_real_maf) {
    # Merge with frequency data
    final_data <- merge(final_data, freq_data[, .(SNP, MAF)], 
                       by.x = "SNP_ID", by.y = "SNP", all.x = TRUE)
    
    # Check for missing MAF
    missing_maf <- sum(is.na(final_data$MAF))
    if (missing_maf > 0) {
        cat(sprintf("Warning: %d SNPs missing MAF data, using 0.1 as placeholder\n", missing_maf))
        final_data[is.na(MAF), MAF := 0.1]
    }
    
    cat("Using real MAF values from frequency file\n")
} else {
    # Use placeholder MAF
    final_data[, MAF := 0.1]
    cat("Using placeholder MAF values (0.1 for all SNPs)\n")
}

# Classify rare vs common
final_data[, Rare_Common := ifelse(MAF < maf_cutoff, "Rare", "Common")]

# Reorder columns
col_order <- c("Gene", "Gene_Pvalue", "SNP_ID", "CHR", "POS", "A1", "A2", "MAF", "Rare_Common")
final_data <- final_data[, ..col_order]

# Sort by gene p-value, then by chromosome and position
setorder(final_data, Gene_Pvalue, CHR, POS)

cat(sprintf("Final dataset: %d actually tested SNPs across %d significant genes\n", 
           nrow(final_data), length(unique(final_data$Gene))))

cat("\nStep 9: Creating summary statistics...\n")

# Create gene-level summary with actual tested SNP counts
gene_summary <- final_data[, .(
    Gene_Pvalue = first(Gene_Pvalue),
    Actually_Tested_SNPs = .N,
    Rare_SNPs = sum(Rare_Common == "Rare"),
    Common_SNPs = sum(Rare_Common == "Common"),
    Chromosomes = paste(unique(CHR), collapse = ","),
    Position_Range = paste(min(POS), max(POS), sep = "-")
), by = Gene]

# Add original SKAT counts for comparison
gene_summary <- merge(gene_summary, 
                     sig_genes[, .(SetID, N.Marker.All, N.Marker.Test, N.Marker.Rare, N.Marker.Common)],
                     by.x = "Gene", by.y = "SetID", all.x = TRUE)

setorder(gene_summary, Gene_Pvalue)

cat("Gene-level summary (actual tested vs expected):\n")
for (i in 1:nrow(gene_summary)) {
    row <- gene_summary[i]
    cat(sprintf("  %s: p=%.2e\n", row$Gene, row$Gene_Pvalue))
    cat(sprintf("    Actually tested: %d SNPs (%d rare, %d common)\n", 
               row$Actually_Tested_SNPs, row$Rare_SNPs, row$Common_SNPs))
    cat(sprintf("    SKAT reported: %d/%d tested (%d rare, %d common)\n",
               row$N.Marker.Test, row$N.Marker.All, row$N.Marker.Rare, row$N.Marker.Common))
    
    if (row$Actually_Tested_SNPs != row$N.Marker.Test) {
        cat(sprintf("    *** MISMATCH: Found %d, expected %d ***\n", 
                   row$Actually_Tested_SNPs, row$N.Marker.Test))
    }
}

cat("\nStep 10: Saving results...\n")

# Save detailed results
output_file <- "significant_genes_actually_tested_snps.txt"
fwrite(final_data, output_file, sep = "\t", quote = FALSE)
cat(sprintf("Detailed results saved: %s\n", output_file))

# Save summary results
summary_file <- "significant_genes_tested_summary.txt"
fwrite(gene_summary, summary_file, sep = "\t", quote = FALSE)
cat(sprintf("Gene summary saved: %s\n", summary_file))

# Create overall summary
overall_summary <- data.frame(
    Total_Significant_Genes = nrow(gene_summary),
    Total_Actually_Tested_SNPs = nrow(final_data),
    Total_Rare_SNPs = sum(final_data$Rare_Common == "Rare"),
    Total_Common_SNPs = sum(final_data$Rare_Common == "Common"),
    Most_Significant_Pvalue = min(final_data$Gene_Pvalue),
    Most_Significant_Gene = final_data$Gene[which.min(final_data$Gene_Pvalue)],
    Analysis_Date = Sys.time(),
    Significance_Threshold = sig_threshold,
    MAF_Cutoff = maf_cutoff,
    Used_Real_MAF = use_real_maf
)

overall_summary_file <- "actually_tested_extraction_summary.txt"
fwrite(overall_summary, overall_summary_file, sep = "\t", quote = FALSE)
cat(sprintf("Overall summary saved: %s\n", overall_summary_file))

cat("\n=== SUMMARY ===\n")
cat(sprintf("Significant genes (p < %.0e): %d\n", sig_threshold, nrow(gene_summary)))
cat(sprintf("Total actually tested SNPs: %d\n", nrow(final_data)))
cat(sprintf("Rare SNPs (MAF < %.2f): %d\n", maf_cutoff, sum(final_data$Rare_Common == "Rare")))
cat(sprintf("Common SNPs (MAF >= %.2f): %d\n", maf_cutoff, sum(final_data$Rare_Common == "Common")))

# Check for discrepancies
total_expected <- sum(sig_genes$N.Marker.Test)
total_found <- nrow(final_data)
if (total_expected != total_found) {
    cat(sprintf("\n*** IMPORTANT DISCREPANCY ***\n"))
    cat(sprintf("Expected total tested SNPs: %d\n", total_expected))
    cat(sprintf("Actually found tested SNPs: %d\n", total_found))
    cat(sprintf("Difference: %d\n", total_found - total_expected))
}

cat("\nOutput files created:\n")
cat(sprintf("  - %s (detailed actually tested SNP information)\n", output_file))
cat(sprintf("  - %s (gene-level summary with comparisons)\n", summary_file))
cat(sprintf("  - %s (overall analysis summary)\n", overall_summary_file))

if (!use_real_maf) {
    cat("\n=== TO GET REAL MAF VALUES ===\n")
    cat("Run: plink --bfile plink_maf_filtered --freq --out freq_data\n")
    cat("Then re-run this script for accurate MAF classification.\n")
}

cat("\n=== ANALYSIS COMPLETED ===\n")
cat("Analysis finished at:", as.character(Sys.time()), "\n")

EOF

# Check if R script completed successfully
if [ $? -eq 0 ]; then
    echo
    echo "=== EXTRACTION COMPLETED SUCCESSFULLY ==="
    echo "Analysis finished at: $(date)"
    echo
    echo "Output files:"
    echo "  - significant_genes_actually_tested_snps.txt"
    echo "  - significant_genes_tested_summary.txt" 
    echo "  - actually_tested_extraction_summary.txt"
    echo
    echo "This script extracts the EXACT SNPs that SKAT actually tested,"
    echo "not just the SNPs that were supposed to be tested."
else
    echo "ERROR: R script failed!"
    exit 1
fi
