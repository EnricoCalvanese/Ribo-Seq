#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio
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

# Explicitly create directories
mkdir -p results/predictions

# Create the R script
cat << 'EOF' > predict_uorfs.R
# Load required libraries
library(GenomicFeatures)
library(Biostrings)
library(systemPipeR)

print("Starting uORF prediction analysis...")

# Function to trim sequences to 500 nt
trimSequences <- function(sequences, max_length = 500) {
    trimmed <- sequences
    long_seqs <- width(sequences) > max_length
    if (any(long_seqs)) {
        # For sequences longer than 500nt, keep only the last 500nt 
        # (closest to the CDS)
        trimmed[long_seqs] <- subseq(sequences[long_seqs], 
                                   start = pmax(1, width(sequences[long_seqs]) - max_length + 1),
                                   end = width(sequences[long_seqs]))
        print(paste("Trimmed", sum(long_seqs), "sequences to 500 nt"))
    }
    return(trimmed)
}

# Function to check spacing between uORFs
checkSpacing <- function(orfs, min_distance = 6) {
    if (length(orfs) <= 1) return(TRUE)
    
    # Get all stop positions and start positions
    stops <- end(ranges(orfs))
    starts <- start(ranges(orfs))[-1]  # All starts except the first one
    
    # Check distances between consecutive ORFs
    distances <- starts - stops[-length(stops)]
    return(all(distances >= min_distance))
}

# Load and trim sequences
print("Loading and trimming 5' UTR sequences...")
utr_sequences <- readDNAStringSet("results/sequences/utr_sequences.fa")
print(paste("Loaded", length(utr_sequences), "UTR sequences"))

# Trim sequences to 500 nt
trimmed_sequences <- trimSequences(utr_sequences)

# Predict ORFs
print("Predicting ORFs...")
uorf_predictions <- predORF(trimmed_sequences, 
                          n="all",           
                          mode="orf",        
                          longest_disjoint=TRUE,  
                          strand="sense")    

print(paste("Found predictions for", length(uorf_predictions), "sequences"))

# Process results
print("Processing results...")
results_df <- data.frame(
    transcript_id = rep(names(uorf_predictions), 
                       sapply(uorf_predictions, length)),
    start = unlist(lapply(uorf_predictions, function(x) start(ranges(x)))),
    end = unlist(lapply(uorf_predictions, function(x) end(ranges(x)))),
    width = unlist(lapply(uorf_predictions, function(x) width(ranges(x))))
)

# Verify no overlaps exist
print("Checking for overlaps...")
overlap_count <- 0
for (transcript in unique(results_df$transcript_id)) {
    transcript_orfs <- results_df[results_df$transcript_id == transcript,]
    if (nrow(transcript_orfs) > 1) {
        for (i in 1:(nrow(transcript_orfs)-1)) {
            if (transcript_orfs$end[i] >= transcript_orfs$start[i+1]) {
                overlap_count <- overlap_count + 1
            }
        }
    }
}
print(paste("Found", overlap_count, "overlapping uORF pairs"))

# Check spacing between uORFs
print("Checking uORF spacing...")
spacing_violations <- 0
for (transcript in names(uorf_predictions)) {
    if (!checkSpacing(uorf_predictions[[transcript]])) {
        spacing_violations <- spacing_violations + 1
    }
}
print(paste("Found", spacing_violations, "spacing violations"))

# Add gene IDs and filter
results_df$gene_id <- sub("transcript:([^.]+).*", "\\1", results_df$transcript_id)
valid_orfs <- results_df$width >= 6
results_df <- results_df[valid_orfs,]

# Generate statistics
genes_with_uorfs <- unique(results_df$gene_id)
num_genes_with_uorfs <- length(genes_with_uorfs)
total_uorfs <- nrow(results_df)
uorf_lengths <- results_df$width
genes_per_transcript <- table(results_df$gene_id)

# Create summary report
summary_report <- paste0(
    "uORF Analysis Summary\n",
    "===================\n\n",
    "Input Data:\n",
    "- Using filtered 5' UTR sequences (trimmed to 500 nt)\n",
    "- Total UTR sequences analyzed: ", length(trimmed_sequences), "\n",
    "- Sequences trimmed to 500 nt: ", sum(width(utr_sequences) > 500), "\n\n",
    "Prediction Parameters:\n",
    "- Maximum UTR length considered: 500 nt\n",
    "- Start codons: ATG\n",
    "- Stop codons: TAA, TAG, TGA\n",
    "- Minimum length: 6 nucleotides\n",
    "- No overlapping ORFs allowed\n",
    "- Minimum spacing between uORFs: 6 nt\n",
    "- Sense strand only\n\n",
    "Quality Control:\n",
    "- Overlapping uORF pairs found: ", overlap_count, "\n",
    "- Spacing violations found: ", spacing_violations, "\n\n",
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
    "- Median length: ", round(median(uorf_lengths), 2), " nt\n"
)

# Save results
print("Saving results...")
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
