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

# Load reference data and previous results
print("Loading reference data and previous results...")
ref_uorfs <- read.table("results/xu2017uORFs.txt", 
                       header=TRUE, 
                       stringsAsFactors=FALSE)

previous_matches <- read.table("results/reference_analysis/matching_statistics.txt",
                             header=TRUE, 
                             stringsAsFactors=FALSE)

# Define problem gene sets
extra_pred_genes <- previous_matches$gene_id[
    previous_matches$matched_uorfs > previous_matches$ref_uorfs
]
print(paste("Found", length(extra_pred_genes), 
            "genes with extra predictions to process with disjoint mode"))

# Create reference lookup
ref_by_gene <- split(ref_uorfs, sub("_.*", "", ref_uorfs$uORF_ID))

# Load UTR sequences
print("Loading UTR sequences...")
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")

# Function to find closest length match
find_closest_match <- function(target_length, available_lengths, used_lengths) {
    if(length(available_lengths) == 0) return(NULL)
    # Exclude already used lengths
    available_lengths <- setdiff(available_lengths, used_lengths)
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
        end = end(ranges(predictions[[1]])),
        index = seq_along(predictions[[1]])
    )
    
    # Initialize results
    selected_indices <- integer(0)
    used_lengths <- numeric(0)
    
    # First pass: find exact matches
    for(ref_len in ref_lengths) {
        exact_matches <- which(pred_info$length == ref_len)
        if(length(exact_matches) > 0) {
            # Take the first available exact match
            for(match_idx in exact_matches) {
                if(!match_idx %in% selected_indices) {
                    selected_indices <- c(selected_indices, match_idx)
                    used_lengths <- c(used_lengths, ref_len)
                    break
                }
            }
        }
    }
    
    # Second pass: find closest matches for remaining lengths
    remaining_ref <- setdiff(ref_lengths, used_lengths)
    for(ref_len in remaining_ref) {
        available_lengths <- pred_info$length
        closest <- find_closest_match(ref_len, available_lengths, used_lengths)
        if(!is.null(closest)) {
            match_idx <- which(pred_info$length == closest)[1]
            if(!match_idx %in% selected_indices) {
                selected_indices <- c(selected_indices, match_idx)
                used_lengths <- c(used_lengths, closest)
            }
        }
    }
    
    # Ensure we have exactly the right number of uORFs
    if(length(selected_indices) > length(ref_lengths)) {
        selected_indices <- selected_indices[1:length(ref_lengths)]
    }
    
    # Return selected ORFs
    if(length(selected_indices) > 0) {
        return(predictions[[1]][selected_indices])
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

# Track progress
total_genes <- length(names(ref_by_gene))
current_gene <- 0

for(gene_id in names(ref_by_gene)) {
    current_gene <- current_gene + 1
    if(current_gene %% 100 == 0) {
        print(paste("Processing gene", current_gene, "of", total_genes))
    }
    
    # Find corresponding sequence
    seq_idx <- grep(paste0("transcript:", gene_id, "\\."), names(utr_sequences))
    if(length(seq_idx) == 0) next
    
    # Determine processing strategy
    use_disjoint <- gene_id %in% extra_pred_genes
    
    # Process gene
    results <- process_gene(utr_sequences[seq_idx], 
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
}

# Create final output data frame
print("Creating final output...")
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

# Sort results to match reference format
results_df <- results_df[order(results_df$uORF_ID),]

# Generate summary report
print("Generating summary report...")
summary_report <- paste0(
    "Final uORF Matching Results\n",
    "========================\n\n",
    "Overall Statistics:\n",
    "- Total reference genes: ", length(names(ref_by_gene)), "\n",
    "- Total genes processed: ", nrow(processing_stats), "\n",
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
