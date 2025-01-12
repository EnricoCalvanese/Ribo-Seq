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
# Load required libraries
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)
library(GenomicRanges)

print("Starting uORF prediction with reference length matching...")

# Function to read and parse reference uORF data
read_reference_uorfs <- function(file_path) {
    # Read the reference file
    ref_data <- read.table(file_path, header=TRUE, stringsAsFactors=FALSE)
    
    # Split uORF_ID into gene and uORF number
    ref_data$gene_id <- sub("_.*$", "", ref_data$uORF_ID)
    
    # Create a list where each element contains the expected lengths for a gene
    ref_list <- split(ref_data$Length_uORF, ref_data$gene_id)
    
    # Sort lengths within each gene for easier matching
    ref_list <- lapply(ref_list, sort)
    
    return(ref_list)
}

# Function to find best matching uORFs for a gene
match_uorf_lengths <- function(predicted_orfs, expected_lengths, tolerance=0) {
    if (length(predicted_orfs) == 0) return(NULL)
    
    # Get widths of predicted ORFs
    pred_lengths <- width(ranges(predicted_orfs))
    
    # Initialize result vector
    matches <- integer(0)
    used_pred <- logical(length(pred_lengths))
    
    # For each expected length, find the closest match
    for (exp_len in expected_lengths) {
        # Calculate differences for unused predictions
        diffs <- abs(pred_lengths[!used_pred] - exp_len)
        
        # If we found matches within tolerance
        if (any(diffs <= tolerance)) {
            # Get index of best match among unused predictions
            best_idx <- which(!used_pred)[which.min(diffs)]
            matches <- c(matches, best_idx)
            used_pred[best_idx] <- TRUE
        }
    }
    
    if (length(matches) > 0) {
        return(predicted_orfs[matches])
    } else {
        return(NULL)
    }
}

# Load reference uORF data
print("Loading reference uORF data...")
ref_uorfs <- read_reference_uorfs("results/xu2017uORFs.txt")
print(paste("Loaded", length(ref_uorfs), "reference genes"))

# Load UTR sequences
print("Loading UTR sequences...")
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")
print(paste("Processing", length(utr_sequences), "UTR sequences"))

# Predict initial ORFs
print("Predicting uORFs...")
uorf_predictions <- predORF(utr_sequences, 
                          n="all",
                          mode="orf",
                          longest_disjoint=TRUE,
                          strand="sense")

# Match predictions to reference lengths
print("Matching predictions to reference lengths...")
matched_predictions <- list()
matching_stats <- list(
    total_matches = 0,
    perfect_genes = 0,
    partial_genes = 0,
    no_match_genes = 0,
    missing_genes = character(0)
)

for (gene_id in names(ref_uorfs)) {
    # Get transcript ID for this gene (if available)
    transcript_pattern <- paste0(">transcript:", gene_id, "\\.")
    transcript_idx <- grep(transcript_pattern, names(utr_sequences))
    
    if (length(transcript_idx) > 0) {
        # Get predictions for this transcript
        transcript_id <- sub(">", "", names(utr_sequences)[transcript_idx[1]])
        predictions <- uorf_predictions[[transcript_id]]
        
        # Match to reference lengths
        matched_orfs <- match_uorf_lengths(predictions, ref_uorfs[[gene_id]], tolerance=0)
        
        if (!is.null(matched_orfs)) {
            matched_predictions[[transcript_id]] <- matched_orfs
            if (length(matched_orfs) == length(ref_uorfs[[gene_id]])) {
                matching_stats$perfect_genes <- matching_stats$perfect_genes + 1
            } else {
                matching_stats$partial_genes <- matching_stats$partial_genes + 1
            }
            matching_stats$total_matches <- matching_stats$total_matches + length(matched_orfs)
        } else {
            matching_stats$no_match_genes <- matching_stats$no_match_genes + 1
        }
    } else {
        matching_stats$missing_genes <- c(matching_stats$missing_genes, gene_id)
    }
}

# Create results data frame
results_df <- data.frame(
    transcript_id = rep(names(matched_predictions), 
                       sapply(matched_predictions, length)),
    start = unlist(lapply(matched_predictions, function(x) start(ranges(x)))),
    end = unlist(lapply(matched_predictions, function(x) end(ranges(x)))),
    width = unlist(lapply(matched_predictions, function(x) width(ranges(x))))
)

if (nrow(results_df) > 0) {
    results_df$gene_id <- sub("transcript:([^.]+).*", "\\1", results_df$transcript_id)
}

# Create detailed summary report
summary_report <- paste0(
    "uORF Prediction Results with Length Matching\n",
    "=========================================\n\n",
    "Reference Data:\n",
    "- Total reference genes: ", length(ref_uorfs), "\n",
    "- Total reference uORFs: ", sum(sapply(ref_uorfs, length)), "\n\n",
    "Matching Results:\n",
    "- Total matched uORFs: ", matching_stats$total_matches, "\n",
    "- Genes with all uORFs matched: ", matching_stats$perfect_genes, "\n",
    "- Genes with partial matches: ", matching_stats$partial_genes, "\n",
    "- Genes with no matches: ", matching_stats$no_match_genes, "\n",
    "- Missing genes: ", length(matching_stats$missing_genes), "\n\n",
    "Length Distribution of Matched uORFs:\n",
    "- Minimum length: ", min(results_df$width), " nt\n",
    "- Maximum length: ", max(results_df$width), " nt\n",
    "- Mean length: ", round(mean(results_df$width), 2), " nt\n",
    "- Median length: ", round(median(results_df$width), 2), " nt\n\n",
    "Missing Genes:\n",
    paste(matching_stats$missing_genes, collapse="\n")
)

# Save results
print("Saving results...")
saveRDS(matched_predictions, 
        "results/reference_analysis/matched_uorf_predictions.rds")
write.table(results_df,
            "results/reference_analysis/matched_uorf_results.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
writeLines(summary_report, 
          "results/reference_analysis/matching_summary.txt")

# Save list of problematic genes for investigation
writeLines(matching_stats$missing_genes,
           "results/reference_analysis/missing_genes.txt")

print("Analysis complete. Check results/reference_analysis/ for output files.")
EOF

# Make R script executable
chmod +x predict_reference_uorfs.R

# Run R script
Rscript predict_reference_uorfs.R
