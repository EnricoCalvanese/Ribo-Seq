#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=uorf_predict
#SBATCH --output=uorf_predict_%j.out
#SBATCH --error=uorf_predict_%j.err

module load r 

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Create the R script
cat << 'EOF' > predict_uorfs.R
# Load required libraries
library(GenomicFeatures)
library(Biostrings)

print("Starting uORF prediction analysis...")

# Load our previously filtered and processed UTR sequences
print("Loading filtered 5' UTR sequences...")
utr_sequences <- readDNAStringSet("results/sequences/utr_sequences.fa")
print(paste("Loaded", length(utr_sequences), "UTR sequences"))

# Predict ORFs using predORF
print("Predicting ORFs...")
uorf_predictions <- predORF(utr_sequences, 
                          n="all",           # Find all possible ORFs
                          mode="orf",        # Look for complete ORFs
                          longest_disjoint=TRUE,  # No overlapping ORFs
                          strand="sense")    # Only look on sense strand

# Process results to get gene-level information
print("Processing results...")

# Function to extract gene ID from transcript ID
get_gene_id <- function(transcript_id) {
    sub("transcript:([^.]+).*", "\\1", transcript_id)
}

# Convert results to data frame for easier processing
results_df <- data.frame(
    transcript_id = rep(names(uorf_predictions), 
                       sapply(uorf_predictions, length)),
    start = unlist(lapply(uorf_predictions, function(x) start(ranges(x)))),
    end = unlist(lapply(uorf_predictions, function(x) end(ranges(x)))),
    width = unlist(lapply(uorf_predictions, function(x) width(ranges(x))))
)

# Add gene IDs
results_df$gene_id <- get_gene_id(results_df$transcript_id)

# Filter ORFs
valid_orfs <- results_df$width >= 6  # Minimum length requirement
results_df <- results_df[valid_orfs,]

# Count unique genes and their uORFs
genes_with_uorfs <- unique(results_df$gene_id)
num_genes_with_uorfs <- length(genes_with_uorfs)
total_uorfs <- nrow(results_df)

# Calculate statistics
uorf_lengths <- results_df$width
genes_per_transcript <- table(results_df$gene_id)

# Create summary report
summary_report <- paste0(
    "uORF Analysis Summary\n",
    "===================\n\n",
    "Input Data:\n",
    "- Using filtered 5' UTR sequences from previous steps\n",
    "- Total UTR sequences analyzed: ", length(utr_sequences), "\n\n",
    "Prediction Parameters:\n",
    "- Start codons: ATG\n",
    "- Stop codons: TAA, TAG, TGA\n",
    "- Minimum length: 6 nucleotides\n",
    "- No overlapping ORFs allowed\n",
    "- Sense strand only\n\n",
    "uORF Prediction Results:\n",
    "- Number of genes with predicted uORFs: ", num_genes_with_uorfs, "\n",
    "- Total number of predicted uORFs: ", total_uorfs, "\n",
    "- Average uORFs per gene: ", round(total_uorfs/num_genes_with_uorfs, 2), "\n\n",
    "Gene Distribution:\n",
    "- Genes with 1 uORF: ", sum(genes_per_transcript == 1), "\n",
    "- Genes with 2 uORFs: ", sum(genes_per_transcript == 2), "\n",
    "- Genes with 3+ uORFs: ", sum(genes_per_transcript >= 3), "\n\n",
    "Length Distribution of Predicted uORFs:\n",
    "- Minimum length: ", min(uorf_lengths), " nt\n",
    "- Maximum length: ", max(uorf_lengths), " nt\n",
    "- Mean length: ", round(mean(uorf_lengths), 2), " nt\n",
    "- Median length: ", round(median(uorf_lengths), 2), " nt\n\n",
    "Validation Notes:\n",
    "- Using predORF function on previously filtered UTRs\n",
    "- All ORFs are complete (start to stop)\n",
    "- No overlapping ORFs included\n",
    "- All ORFs >= 6 nucleotides\n"
)

# Save results
saveRDS(results_df, "results/predictions/uorf_predictions.rds")
writeLines(summary_report, "results/predictions/prediction_summary.txt")
write.table(genes_with_uorfs, "results/predictions/genes_with_uorfs.txt", 
            row.names=FALSE, col.names=FALSE, quote=FALSE)

print("Analysis complete. Check results/predictions/ for output files.")
EOF

# Make R script executable
chmod +x predict_uorfs.R

# Run R script
Rscript predict_uorfs.R
