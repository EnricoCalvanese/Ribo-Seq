#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=predict_ref_uorfs
#SBATCH --output=predict_ref_uorfs_%j.out
#SBATCH --error=predict_ref_uorfs_%j.err
module load r
# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Create the R script for predictions
cat << 'EOF' > predict_reference_uorfs.R
# Load required libraries for uORF prediction
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)

print("Starting uORF prediction on representative UTRs...")

# Function to validate ORFs with length constraint
validateORFs <- function(orfs, max_length = 1100) {
    if (length(orfs) == 0) return(orfs)
    
    # Get widths of all ORFs
    widths <- width(ranges(orfs))
    
    # Keep only ORFs within length constraints
    valid_orfs <- orfs[widths >= 6 & widths <= max_length]
    
    return(valid_orfs)
}

# Load filtered UTR sequences
print("Loading representative UTR sequences...")
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")
print(paste("Processing", length(utr_sequences), "UTR sequences"))

# Extract gene IDs for reference
gene_ids <- sub("transcript:([^.]+).*", "\\1", names(utr_sequences))
print(paste("Number of unique genes:", length(unique(gene_ids))))

# Predict ORFs using predORF
print("Predicting uORFs...")
uorf_predictions <- predORF(utr_sequences, 
                          n="all",           # Find all possible ORFs
                          mode="orf",        # Complete ORFs only
                          longest_disjoint=TRUE,  # No overlapping ORFs
                          strand="sense")    # Sense strand only

# Apply length validation
print("Validating predictions...")
validated_predictions <- lapply(uorf_predictions, validateORFs)

# Create comprehensive results data frame
results_df <- data.frame(
    transcript_id = rep(names(validated_predictions), 
                       sapply(validated_predictions, length)),
    start = unlist(lapply(validated_predictions, function(x) {
        if(length(x) == 0) numeric(0) else start(ranges(x))
    })),
    end = unlist(lapply(validated_predictions, function(x) {
        if(length(x) == 0) numeric(0) else end(ranges(x))
    })),
    width = unlist(lapply(validated_predictions, function(x) {
        if(length(x) == 0) numeric(0) else width(ranges(x))
    }))
)

# Add gene IDs and sort by position
if (nrow(results_df) > 0) {
    results_df$gene_id <- sub("transcript:([^.]+).*", "\\1", results_df$transcript_id)
    results_df <- results_df[order(results_df$gene_id, results_df$start), ]
}

# Calculate detailed statistics
genes_with_uorfs <- unique(results_df$gene_id)
num_genes_with_uorfs <- length(genes_with_uorfs)
total_uorfs <- nrow(results_df)
genes_per_transcript <- table(results_df$gene_id)

# Create detailed summary report
summary_report <- paste0(
    "uORF Prediction Results for Reference Genes\n",
    "========================================\n\n",
    "Input Data:\n",
    "- Total UTR sequences analyzed: ", length(utr_sequences), "\n",
    "- All sequences are representative UTRs (one per gene)\n\n",
    "Prediction Parameters:\n",
    "- Maximum uORF length: 1100 nucleotides\n",
    "- Minimum length: 6 nucleotides\n",
    "- No overlapping ORFs allowed\n",
    "- Sense strand only\n\n",
    "Results:\n",
    "- Number of genes with predicted uORFs: ", num_genes_with_uorfs, "\n",
    "- Total number of predicted uORFs: ", total_uorfs, "\n",
    "- Average uORFs per gene: ", round(total_uorfs/num_genes_with_uorfs, 2), "\n\n",
    "Gene Distribution:\n",
    "- Genes with 1 uORF: ", sum(genes_per_transcript == 1), "\n",
    "- Genes with 2 uORFs: ", sum(genes_per_transcript == 2), "\n",
    "- Genes with 3+ uORFs: ", sum(genes_per_transcript >= 3), "\n\n",
    "Length Distribution of Predicted uORFs:\n",
    "- Minimum length: ", min(results_df$width), " nt\n",
    "- Maximum length: ", max(results_df$width), " nt\n",
    "- Mean length: ", round(mean(results_df$width), 2), " nt\n",
    "- Median length: ", round(median(results_df$width), 2), " nt\n\n",
    "Comparison with Reference Study:\n",
    "- Expected genes with uORFs: 6033\n",
    "- Our prediction: ", num_genes_with_uorfs, "\n",
    "- Difference: ", num_genes_with_uorfs - 6033, " genes\n",
    "- Match percentage: ", 
        round(100 * min(num_genes_with_uorfs, 6033)/max(num_genes_with_uorfs, 6033), 2), "%\n"
)

# Save results
print("Saving results...")
saveRDS(validated_predictions, 
        "results/reference_analysis/uorf_predictions.rds")
write.table(results_df,
            "results/reference_analysis/uorf_results.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
writeLines(summary_report, 
          "results/reference_analysis/prediction_summary.txt")
write.table(genes_with_uorfs,
            "results/reference_analysis/predicted_genes.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

print("Analysis complete. Check results/reference_analysis/ for output files.")
EOF

# Make R script executable
chmod +x predict_reference_uorfs.R

# Run R script
Rscript predict_reference_uorfs.R
