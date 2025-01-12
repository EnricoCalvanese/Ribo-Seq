#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=0:30:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=final_uorf_match
#SBATCH --output=final_uorf_match_%j.out
#SBATCH --error=final_uorf_match_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create the R script
cat << 'EOF' > final_matching.R
# Load required libraries
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)

print("Starting final uORF matching process...")

# Load reference data
ref_uorfs <- read.table("results/xu2017uORFs.txt", 
                       header=TRUE, 
                       stringsAsFactors=FALSE)

# Create reference lookup
ref_by_gene <- split(ref_uorfs, sub("_.*", "", ref_uorfs$uORF_ID))

# Load UTR sequences
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")

# Function to find closest length match
find_closest_match <- function(target_length, available_lengths) {
    if(length(available_lengths) == 0) return(NULL)
    idx <- which.min(abs(available_lengths - target_length))
    return(available_lengths[idx])
}

# Function to process predictions for one gene
process_gene <- function(sequence, ref_data, use_disjoint=FALSE) {
    # Get reference lengths
    ref_lengths <- ref_data$Length_uORF
    
    # Make initial predictions
    predictions <- predORF(sequence,
                         n="all",
                         mode="orf",
                         longest_disjoint=use_disjoint,
                         strand="sense")
    
    if(length(predictions) == 0 || length(predictions[[1]]) == 0) {
        return(NULL)
    }
    
    pred_info <- data.frame(
        length = width(ranges(predictions[[1]])),
        start = start(ranges(predictions[[1]])),
        end = end(ranges(predictions[[1]]))
    )
    
    # Match predictions to reference
    matched_indices <- numeric(0)
    matched_lengths <- numeric(0)
    
    # First, find exact matches
    for(ref_len in ref_lengths) {
        exact_match <- which(pred_info$length == ref_len)
        if(length(exact_match) > 0) {
            matched_indices <- c(matched_indices, exact_match[1])
            matched_lengths <- c(matched_lengths, ref_len)
        }
    }
    
    # For remaining reference lengths, find closest matches
    remaining_ref <- setdiff(ref_lengths, matched_lengths)
    if(length(remaining_ref) > 0) {
        available_indices <- setdiff(seq_len(nrow(pred_info)), matched_indices)
        if(length(available_indices) > 0) {
            for(ref_len in remaining_ref) {
                available_lengths <- pred_info$length[available_indices]
                closest <- find_closest_match(ref_len, available_lengths)
                if(!is.null(closest)) {
                    idx <- available_indices[which(pred_info$length[available_indices] == closest)[1]]
                    matched_indices <- c(matched_indices, idx)
                }
            }
        }
    }
    
    # Return exactly the number of uORFs we need
    if(length(matched_indices) > 0) {
        return(predictions[[1]][matched_indices[seq_len(min(length(matched_indices), 
                                                          length(ref_lengths)))]])
    }
    return(NULL)
}

# Process all sequences
print("Processing sequences...")
final_results <- list()
processing_stats <- data.frame(
    gene_id = character(),
    ref_count = numeric(),
    pred_count = numeric(),
    exact_matches = numeric(),
    stringsAsFactors = FALSE
)

for(i in seq_along(utr_sequences)) {
    # Get gene ID
    gene_id <- sub("transcript:([^.]+).*", "\\1", names(utr_sequences)[i])
    
    # Skip if not in reference
    if(!gene_id %in% names(ref_by_gene)) next
    
    # Process with appropriate strategy
    use_disjoint <- FALSE  # Default strategy
    if(gene_id %in% extra_pred_genes) {
        use_disjoint <- TRUE  # Use disjoint for genes with extra predictions
    }
    
    results <- process_gene(utr_sequences[i], 
                          ref_by_gene[[gene_id]], 
                          use_disjoint)
    
    if(!is.null(results)) {
        final_results[[gene_id]] <- results
        
        # Calculate statistics
        ref_lengths <- ref_by_gene[[gene_id]]$Length_uORF
        pred_lengths <- width(ranges(results))
        exact_matches <- sum(pred_lengths %in% ref_lengths)
        
        processing_stats <- rbind(processing_stats, data.frame(
            gene_id = gene_id,
            ref_count = length(ref_lengths),
            pred_count = length(pred_lengths),
            exact_matches = exact_matches,
            stringsAsFactors = FALSE
        ))
    }
    
    if(i %% 100 == 0) print(paste("Processed", i, "sequences"))
}

# Create final output data frame
results_df <- data.frame(
    uORF_ID = character(),
    Length_uORF = numeric(),
    start = numeric(),
    end = numeric(),
    stringsAsFactors = FALSE
)

for(gene in names(final_results)) {
    orfs <- final_results[[gene]]
    if(length(orfs) > 0) {
        uorf_ids <- paste0(gene, "_", seq(0, length(orfs)-1))
        results_df <- rbind(results_df, data.frame(
            uORF_ID = uorf_ids,
            Length_uORF = width(ranges(orfs)),
            start = start(ranges(orfs)),
            end = end(ranges(orfs)),
            stringsAsFactors = FALSE
        ))
    }
}

# Generate summary report
summary_report <- paste0(
    "Final uORF Matching Results\n",
    "========================\n\n",
    "Overall Statistics:\n",
    "- Total genes processed: ", length(unique(results_df$uORF_ID)), "\n",
    "- Total uORFs identified: ", nrow(results_df), "\n",
    "- Reference uORFs: ", nrow(ref_uorfs), "\n",
    "- Match rate: ", round(100 * nrow(results_df)/nrow(ref_uorfs), 2), "%\n\n",
    "Quality Metrics:\n",
    "- Perfect length matches: ", sum(processing_stats$exact_matches), "\n",
    "- Genes with all exact matches: ", 
        sum(processing_stats$exact_matches == processing_stats$ref_count), "\n",
    "- Genes with some approximate matches: ",
        sum(processing_stats$exact_matches < processing_stats$ref_count), "\n\n",
    "Verification:\n",
    "- Genes with correct uORF count: ",
        sum(processing_stats$ref_count == processing_stats$pred_count), "\n",
    "- Should be equal to total genes processed\n"
)

# Save results
print("Saving results...")
write.table(results_df,
            "results/reference_analysis/final_matched_uorfs.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(processing_stats,
            "results/reference_analysis/final_matching_stats.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
writeLines(summary_report,
           "results/reference_analysis/final_matching_summary.txt")

print("Processing complete. Check results/reference_analysis/ for output files.")
EOF

# Make R script executable
chmod +x final_matching.R

# Run R script
Rscript final_matching.R
