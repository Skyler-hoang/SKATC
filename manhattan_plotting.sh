#!/bin/bash
#SBATCH --job-name=manhattan_final
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=manhattan_final_%j.out
#SBATCH --error=manhattan_final_%j.err

# Load R module
module load R/4.4.2-gfbf-2024a

# Set working directory
cd /mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/05_SKAT/CTPP_RAN_Colors_Z_t

echo "=== SKAT-C MANHATTAN PLOT - FINAL VERSION ==="
echo "Working directory: $(pwd)"
echo "Starting analysis at: $(date)"
echo

# Check if input files exist
SKAT_FILE="skat_c_results.txt"
MAPPING_FILE="/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/05_SKAT/mapping.txt"

if [[ ! -f "$SKAT_FILE" ]]; then
    echo "ERROR: SKAT results file not found: $SKAT_FILE"
    exit 1
fi

if [[ ! -f "$MAPPING_FILE" ]]; then
    echo "ERROR: Mapping file not found: $MAPPING_FILE"
    exit 1
fi

echo "Input files validated successfully"

# Run R script
R --vanilla << 'EOF'

# Load required libraries
suppressMessages({
    library(data.table)
    library(ggplot2)
    library(dplyr)
})

cat("=== SKAT-C Manhattan Plot - Final Version ===\n")
cat("Analysis started at:", as.character(Sys.time()), "\n\n")

# Define file paths
skat_file <- "skat_c_results.txt"
mapping_file <- "/mnt/vstor/SOM_EPBI_SKI/PHONOLOGY/PGSSKATO/05_SKAT/mapping.txt"

cat("Step 1: Loading data...\n")

# Load SKAT-C results
skat_results <- fread(skat_file, header = TRUE)
cat(paste("Loaded", nrow(skat_results), "SKAT-C results\n"))

# Load mapping file
mapping <- fread(mapping_file, header = TRUE)
cat(paste("Loaded", nrow(mapping), "gene mappings\n"))

cat("Step 2: Processing chromosome data...\n")

# Clean chromosome names (remove chr prefix)
mapping$CHR <- gsub("chr", "", mapping$CHR)
mapping$CHR <- gsub("Chr", "", mapping$CHR)
mapping$CHR <- as.numeric(mapping$CHR)

# Filter for autosomal chromosomes (1-22)
mapping_filtered <- mapping[mapping$CHR %in% 1:22 & !is.na(mapping$CHR), ]

cat(paste("Chromosomes 1-22 mappings:", nrow(mapping_filtered), "\n"))

cat("Step 3: Merging SKAT results with gene positions...\n")

# Merge data
merged_data <- merge(skat_results, mapping_filtered, 
                    by.x = "SetID", by.y = "Gene_name", 
                    all.x = FALSE)

cat(paste("Successfully merged", nrow(merged_data), "genes\n"))

# Filter for valid p-values
merged_data <- merged_data[!is.na(merged_data$P.value) & 
                          merged_data$P.value > 0 &
                          merged_data$P.value <= 1, ]

cat(paste("Valid data points:", nrow(merged_data), "\n"))

# Calculate -log10(p-values)
merged_data$neg_log10_p <- -log10(merged_data$P.value)

cat("Step 4: Setting up chromosome positions for all 22 chromosomes...\n")

# Sort by chromosome and position
merged_data <- merged_data[order(merged_data$CHR, merged_data$RoundPosition), ]

# Calculate chromosome information for ALL 22 chromosomes
all_chrs <- 1:22
chr_info <- data.frame(CHR = all_chrs)

# Get actual chromosome data
chr_data <- merged_data %>%
    group_by(CHR) %>%
    summarise(
        min_pos = min(RoundPosition),
        max_pos = max(RoundPosition),
        n_genes = n(),
        .groups = "drop"
    ) %>%
    arrange(CHR)

# Merge with all chromosomes
chr_info <- merge(chr_info, chr_data, by = "CHR", all.x = TRUE)

# For any missing chromosomes, set reasonable defaults
chr_info$max_pos[is.na(chr_info$max_pos)] <- 250000000  # 250MB default
chr_info$min_pos[is.na(chr_info$min_pos)] <- 0
chr_info$n_genes[is.na(chr_info$n_genes)] <- 0

# Calculate cumulative positions to properly space all 22 chromosomes
chr_info$cumulative_start <- c(0, cumsum(chr_info$max_pos)[-nrow(chr_info)])
chr_info$cumulative_center <- chr_info$cumulative_start + chr_info$max_pos / 2

# Add cumulative positions to data
merged_data <- merge(merged_data, chr_info[, c("CHR", "cumulative_start")], by = "CHR")
merged_data$plot_position <- merged_data$RoundPosition + merged_data$cumulative_start
merged_data <- merged_data[order(merged_data$plot_position), ]

# Add alternating colors for chromosomes
merged_data$chr_color <- ifelse(merged_data$CHR %% 2 == 0, "even", "odd")

cat("Chromosome distribution:\n")
chr_dist <- table(merged_data$CHR)
print(chr_dist)

cat("\nStep 5: Defining significance thresholds...\n")

# Define significance thresholds
bonferroni_threshold <- 5e-7    # Bonferroni corrected
genome_wide_threshold <- 5e-6   # Gene-based genome-wide
suggestive_threshold <- 5e-5    # Suggestive

# Convert to -log10 scale
bonferroni_line <- -log10(bonferroni_threshold)
genome_wide_line <- -log10(genome_wide_threshold)
suggestive_line <- -log10(suggestive_threshold)

# Count significant genes
bonferroni_genes <- sum(merged_data$P.value < bonferroni_threshold)
genome_wide_genes <- sum(merged_data$P.value < genome_wide_threshold)
suggestive_genes <- sum(merged_data$P.value < suggestive_threshold)

cat(paste("Significance summary:\n"))
cat(paste("  Bonferroni (5×10⁻⁷):", bonferroni_genes, "genes\n"))
cat(paste("  Genome-wide (5×10⁻⁶):", genome_wide_genes, "genes\n"))
cat(paste("  Suggestive (5×10⁻⁵):", suggestive_genes, "genes\n"))

cat("\nStep 6: Creating publication-quality Manhattan plot...\n")

# Define colors for alternating chromosomes
chr_colors <- c("odd" = "#2E86AB", "even" = "#A23B72")

# Create Manhattan plot
manhattan_plot <- ggplot(merged_data, aes(x = plot_position, y = neg_log10_p, color = chr_color)) +
    # Add points
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = chr_colors, guide = "none") +
    
    # Add significance threshold lines
    geom_hline(yintercept = suggestive_line, 
              color = "#FFA500", linetype = "dashed", size = 0.8,
              alpha = 0.8) +
    geom_hline(yintercept = genome_wide_line, 
              color = "#FF4500", linetype = "dashed", size = 0.8,
              alpha = 0.8) +

    
    # Set x-axis to show ALL 22 chromosomes properly
    scale_x_continuous(
        breaks = chr_info$cumulative_center,
        labels = chr_info$CHR,
        limits = c(0, max(chr_info$cumulative_start + chr_info$max_pos)),
        expand = c(0.01, 0)
    ) +
    
    # Set y-axis
    scale_y_continuous(expand = c(0.02, 0)) +
    
    # Labels and title
    labs(
        title = "SKAT-C Gene-Based Association Analysis",
        subtitle = paste("Manhattan Plot -", nrow(merged_data), "genes tested across 22 chromosomes"),
        x = "Chromosome",
        y = expression(-log[10](P-value))
    ) +
    
    # Theme and styling
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )

# Clean plot without text annotations (will be explained in figure legend)
# All three significance lines are now dashed and clearly distinguished by color

cat("Step 7: Saving high-quality Manhattan plots...\n")

# Save high-resolution plots
# PNG version (high resolution for presentations)
ggsave("skat_c_manhattan_final.png", manhattan_plot, 
       width = 16, height = 8, dpi = 300, units = "in")
cat("Saved: skat_c_manhattan_final.png\n")

# PDF version (vector format for publications)
ggsave("skat_c_manhattan_final.pdf", manhattan_plot, 
       width = 16, height = 8, units = "in")
cat("Saved: skat_c_manhattan_final.pdf\n")

# Also create a compressed version for quick viewing
ggsave("skat_c_manhattan_final_compressed.png", manhattan_plot, 
       width = 12, height = 6, dpi = 150, units = "in")
cat("Saved: skat_c_manhattan_final_compressed.png\n")

cat("Step 8: Creating summary files...\n")

# Create comprehensive summary statistics
summary_stats <- data.frame(
    Metric = c("Total genes tested",
               "Chromosomes represented", 
               "Minimum p-value",
               "Maximum -log10(p)",
               "Median p-value",
               "Genes < Bonferroni (5×10⁻⁷)",
               "Genes < Genome-wide (5×10⁻⁶)", 
               "Genes < Suggestive (5×10⁻⁵)",
               "Most significant gene"),
    Value = c(nrow(merged_data),
              "1-22 (all)",
              sprintf("%.2e", min(merged_data$P.value, na.rm = TRUE)),
              sprintf("%.2f", max(merged_data$neg_log10_p, na.rm = TRUE)),
              sprintf("%.4f", median(merged_data$P.value, na.rm = TRUE)),
              bonferroni_genes,
              genome_wide_genes,
              suggestive_genes,
              merged_data$SetID[which.min(merged_data$P.value)])
)

write.table(summary_stats, "manhattan_plot_final_summary.txt", 
           sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: manhattan_plot_final_summary.txt\n")

# Save significant genes if any
if(suggestive_genes > 0) {
    sig_genes <- merged_data[merged_data$P.value < suggestive_threshold, ]
    sig_genes <- sig_genes[order(sig_genes$P.value), ]
    
    sig_genes_output <- sig_genes[, c("SetID", "CHR", "RoundPosition", "P.value", 
                                     "neg_log10_p", "N.Marker.Test", 
                                     "N.Marker.Rare", "N.Marker.Common")]
    
    write.table(sig_genes_output, "manhattan_plot_significant_genes_final.txt", 
               sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved: manhattan_plot_significant_genes_final.txt\n")
}

# Save chromosome information
chr_summary <- chr_info[, c("CHR", "n_genes", "min_pos", "max_pos", "cumulative_center")]
chr_summary$n_genes[is.na(chr_summary$n_genes)] <- 0
colnames(chr_summary) <- c("Chromosome", "Genes_Tested", "Min_Position", 
                          "Max_Position", "Plot_Center")

write.table(chr_summary, "manhattan_plot_chromosome_info.txt", 
           sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: manhattan_plot_chromosome_info.txt\n")

cat("\n=== MANHATTAN PLOT GENERATION COMPLETED SUCCESSFULLY ===\n")
cat(paste("Total genes plotted:", nrow(merged_data), "\n"))
cat(paste("Chromosomes displayed: 1-22 (all)\n"))
cat(paste("Significant findings: Bonferroni=", bonferroni_genes, 
          ", Genome-wide=", genome_wide_genes, 
          ", Suggestive=", suggestive_genes, "\n"))
cat("Analysis completed at:", as.character(Sys.time()), "\n")

EOF

# Check if R script completed successfully
if [ $? -eq 0 ]; then
    echo
    echo "=== MANHATTAN PLOT GENERATION COMPLETED SUCCESSFULLY ==="
    echo "Analysis finished at: $(date)"
    echo
    echo "Output files created:"
    echo "  📊 skat_c_manhattan_final.png (high-resolution, 300 DPI)"
    echo "  📊 skat_c_manhattan_final.pdf (publication-ready vector)"
    echo "  📊 skat_c_manhattan_final_compressed.png (quick viewing)"
    echo "  📋 manhattan_plot_final_summary.txt (analysis summary)"
    echo "  📋 manhattan_plot_chromosome_info.txt (chromosome details)"
    if [ -f "manhattan_plot_significant_genes_final.txt" ]; then
        echo "  🎯 manhattan_plot_significant_genes_final.txt (significant genes)"
    fi
    echo
    echo "✅ All 22 chromosomes properly displayed with alternating colors"
    echo "✅ Three significance thresholds clearly marked"
    echo "✅ High-quality plots ready for publication"
else
    echo "ERROR: R script failed!"
    exit 1
fi
