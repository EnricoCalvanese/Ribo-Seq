#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
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
library(systemPipeR)  # Required for workflow management
library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)

print("Starting uORF prediction analysis...")

# Function to trim sequences to 500 nt
trimSequences <- function(sequences, max_length = 500) {
    trimmed <- sequences
    long_seqs <- width(sequences) > max_length
    if (any(long_seqs)) {
        # For sequences longer than 500nt, keep only the last 500nt (closest to CDS)
        trimmed[long_seqs] <- subseq(sequences[long_seqs], 
                                   start = pmax(1, width(sequences[long_seqs]) - max_length + 1),
                                   end = width(sequences[long_seqs]))
    }
    return(trimmed)
}

# Function to validate and filter ORFs within a transcript
validateORFs <- function(orfs, min_spacing = 6) {
    if (length(orfs) <= 1) return(orfs)
    
    # Convert to GRanges for easier manipulation
    gr <- GRanges(seqnames = rep("seq", length(orfs)),
                  ranges = IRanges(start = start(ranges(orfs)),
                                 end = end(ranges(orfs))))
    
    # Check for overlaps
    overlaps <- findOverlaps(gr, drop.self=TRUE, drop.redundant=TRUE)
    if (length(overlaps) > 0) {
        # Remove all overlapping ORFs
        to_remove <- unique(c(queryHits(overlaps), subjectHits(overlaps)))
        gr <- gr[-to_remove]
    }
    
    # Check spacing between consecutive ORFs
    if (length(gr) > 1) {
        # Sort by position
        gr <- sort(gr)
        # Calculate distances between consecutive ORFs
        distances <- start(gr)[-1] - end(gr)[-length(gr)]
        # Keep only ORFs with sufficient spacing
        valid_indices <- c(TRUE, distances >= min_spacing)
        gr <- gr[valid_indices]
    }
    
    # Convert back to original format
    return(orfs[start(ranges(orfs)) %in% start(gr)])
}

# Function to get the most downstream uORF
getMostDownstreamORF <- function(orfs) {
    if (length(orfs) <= 1) return(orfs)
    # Get the last (most downstream) ORF
    orfs[length(orfs)]
}

# Load and trim sequences
print("Loading and trimming 5' UTR sequences...")
utr_sequences <- readDNAStringSet("results/sequences/utr_sequences.fa")
trimmed_sequences <- trimSequences(utr_sequences)

# Predict ORFs
print("Predicting ORFs...")
uorf_predictions <- predORF(trimmed_sequences, 
                          n="all",           
                          mode="orf",        
                          longest_disjoint=TRUE,  
                          strand="sense")

# Apply validation and filtering
print("Applying validation and filtering...")
validated_predictions <- lapply(uorf_predictions, validateORFs)

# Create downstream-only version
downstream_predictions <- lapply(validated_predictions, getMostDownstreamORF)

# Process results for both versions
processResults <- function(predictions, label) {
    results_df <- data.frame(
        transcript_id = rep(names(predictions), 
                          sapply(predictions, length)),
        start = unlist(lapply(predictions, function(x) start(ranges(x)))),
        end = unlist(lapply(predictions, function(x) end(ranges(x)))),
        width = unlist(lapply(predictions, function(x) width(ranges(x))))
    )
    
    if (nrow(results_df) == 0) {
        print(paste("No predictions found for", label))
        return(NULL)
    }
    
    results_df$gene_id <- sub("transcript:([^.]+).*", "\\1", results_df$transcript_id)
    valid_orfs <- results_df$width >= 6
    results_df <- results_df[valid_orfs,]
    
    genes_with_uorfs <- unique(results_df$gene_id)
    num_genes_with_uorfs <- length(genes_with_uorfs)
    total_uorfs <- nrow(results_df)
    genes_per_transcript <- table(results_df$gene_id)
    
    stats <- list(
        num_genes = num_genes_with_uorfs,
        total_uorfs = total_uorfs,
        genes_1_uorf = sum(genes_per_transcript == 1),
        genes_2_uorf = sum(genes_per_transcript == 2),
        genes_3plus_uorf = sum(genes_per_transcript >= 3),
        lengths = results_df$width
    )
    
    return(list(results = results_df, stats = stats))
}

# Process both versions
validated_results <- processResults(validated_predictions, "validated")
downstream_results <- processResults(downstream_predictions, "downstream-only")

# Create summary report
summary_report <- paste0(
    "uORF Analysis Summary with Multiple Filtering Scenarios\n",
    "================================================\n\n",
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
    "Results with All Valid uORFs:\n",
    "----------------------------\n",
    "- Number of genes with predicted uORFs: ", validated_results$stats$num_genes, "\n",
    "- Total number of predicted uORFs: ", validated_results$stats$total_uorfs, "\n",
    "- Average uORFs per gene: ", 
        round(validated_results$stats$total_uorfs/validated_results$stats$num_genes, 2), "\n",
    "Gene Distribution:\n",
    "- Genes with 1 uORF: ", validated_results$stats$genes_1_uorf, "\n",
    "- Genes with 2 uORFs: ", validated_results$stats$genes_2_uorf, "\n",
    "- Genes with 3+ uORFs: ", validated_results$stats$genes_3plus_uorf, "\n",
    "Length Distribution:\n",
    "- Minimum length: ", min(validated_results$stats$lengths), " nt\n",
    "- Maximum length: ", max(validated_results$stats$lengths), " nt\n",
    "- Mean length: ", round(mean(validated_results$stats$lengths), 2), " nt\n",
    "- Median length: ", round(median(validated_results$stats$lengths), 2), " nt\n\n",
    "Results with Only Most Downstream uORF:\n",
    "------------------------------------\n",
    "- Number of genes with predicted uORFs: ", downstream_results$stats$num_genes, "\n",
    "- Total number of predicted uORFs: ", downstream_results$stats$total_uorfs, "\n",
    "Length Distribution:\n",
    "- Minimum length: ", min(downstream_results$stats$lengths), " nt\n",
    "- Maximum length: ", max(downstream_results$stats$lengths), " nt\n",
    "- Mean length: ", round(mean(downstream_results$stats$lengths), 2), " nt\n",
    "- Median length: ", round(median(downstream_results$stats$lengths), 2), " nt\n"
)

# Save results
print("Saving results...")
saveRDS(validated_results$results, "results/predictions/uorf_predictions.rds")
saveRDS(downstream_results$results, "results/predictions/downstream_uorf_predictions.rds")
writeLines(summary_report, "results/predictions/prediction_summary.txt")

print("Analysis complete. Check results/predictions/ for output files.")
EOF

# Make R script executable
chmod +x predict_uorfs.R

# Run R script
Rscript predict_uorfs.R
