#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=match_uorfs
#SBATCH --output=match_uorfs_%j.out
#SBATCH --error=match_uorfs_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create the R script
cat << 'EOF' > match_predictions.R
# Load required libraries
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)

print("Starting uORF prediction and matching process...")

# Load reference uORF data
print("Loading reference uORF data...")
ref_uorfs <- read.table("results/xu2017uORFs.txt", 
                       header=TRUE, 
                       stringsAsFactors=FALSE)

# Process reference data
ref_genes <- unique(sub("_.*", "", ref_uorfs$uORF_ID))
print(paste("Reference contains", length(ref_genes), "unique genes"))

# Create lookup table for reference lengths
ref_lookup <- split(ref_uorfs$Length_uORF, sub("_.*", "", ref_uorfs$uORF_ID))

# Load UTR sequences
print("Loading UTR sequences...")
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")
print(paste("Loaded", length(utr_sequences), "UTR sequences"))

# Function to match predictions with reference lengths
match_predictions <- function(predictions, ref_lengths) {
    if(length(predictions) == 0) return(NULL)
    
    pred_lengths <- width(ranges(predictions))
    matches <- pred_lengths %in% ref_lengths
    
    return(predictions[matches])
}

# Initialize results storage
print("Processing sequences and matching predictions...")
all_results <- list()
match_stats <- data.frame(
    gene_id = character(),
    ref_uorfs = numeric(),
    pred_uorfs = numeric(),
    matched_uorfs = numeric(),
    stringsAsFactors = FALSE
)

# Process each sequence
for(i in seq_along(utr_sequences)) {
    # Extract gene ID
    gene_id <- sub("transcript:([^.]+).*", "\\1", names(utr_sequences)[i])
    
    # Skip if gene not in reference
    if(!gene_id %in% ref_genes) next
    
    # Get reference lengths for this gene
    ref_lengths <- ref_lookup[[gene_id]]
    
    # Predict uORFs
    predictions <- predORF(utr_sequences[i],
                         n="all",
                         mode="orf",
                         longest_disjoint=FALSE,
                         strand="sense")
    
    # Match with reference lengths
    if(length(predictions) > 0) {
        matched_orfs <- match_predictions(predictions[[1]], ref_lengths)
        
        # Store results
        if(!is.null(matched_orfs)) {
            all_results[[gene_id]] <- matched_orfs
        }
        
        # Update statistics
        match_stats <- rbind(match_stats, data.frame(
            gene_id = gene_id,
            ref_uorfs = length(ref_lengths),
            pred_uorfs = length(predictions[[1]]),
            matched_uorfs = if(is.null(matched_orfs)) 0 else length(matched_orfs),
            stringsAsFactors = FALSE
        ))
    }
    
    # Progress update
    if(i %% 100 == 0) {
        print(paste("Processed", i, "sequences"))
    }
}

# Calculate summary statistics
total_matched_uorfs <- sum(match_stats$matched_uorfs)
genes_with_matches <- sum(match_stats$matched_uorfs > 0)
perfect_matches <- sum(match_stats$ref_uorfs == match_stats$matched_uorfs)

# Create output data frame with matched uORFs
results_df <- data.frame(
    uORF_ID = character(),
    Length_uORF = numeric(),
    start = numeric(),
    end = numeric(),
    stringsAsFactors = FALSE
)

for(gene in names(all_results)) {
    orfs <- all_results[[gene]]
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

# Create detailed summary report
summary_report <- paste0(
    "uORF Prediction and Reference Matching Summary\n",
    "===========================================\n\n",
    "Reference Data Summary:\n",
    "- Total unique genes in reference: ", length(ref_genes), "\n",
    "- Total uORFs in reference: ", nrow(ref_uorfs), "\n",
    "- Average uORFs per gene: ", round(nrow(ref_uorfs)/length(ref_genes), 2), "\n\n",
    "Matching Results:\n",
    "- Total matched uORFs found: ", total_matched_uorfs, "\n",
    "- Reference recovery rate: ", 
        round(100 * total_matched_uorfs/nrow(ref_uorfs), 2), "%\n",
    "- Genes with at least one match: ", genes_with_matches, "\n",
    "- Genes with perfect matches: ", perfect_matches, "\n\n",
    "Match Distribution:\n",
    "- Genes with no matches: ", 
        sum(match_stats$matched_uorfs == 0), "\n",
    "- Genes with partial matches: ", 
        sum(match_stats$matched_uorfs > 0 & match_stats$matched_uorfs < match_stats$ref_uorfs), "\n",
    "- Genes with complete matches: ", perfect_matches, "\n\n",
    "Detailed Statistics:\n",
    "- Mean matched uORFs per gene: ", 
        round(mean(match_stats$matched_uorfs), 2), "\n",
    "- Median matched uORFs per gene: ", 
        median(match_stats$matched_uorfs), "\n",
    "- Maximum matched uORFs for a gene: ", 
        max(match_stats$matched_uorfs), "\n"
)

# Save results
print("Saving results...")
saveRDS(all_results, "results/reference_analysis/matched_predictions.rds")
write.table(results_df,
            "results/reference_analysis/matched_uorfs.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(match_stats,
            "results/reference_analysis/matching_statistics.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
writeLines(summary_report,
           "results/reference_analysis/matching_summary.txt")

print("Analysis complete. Check results/reference_analysis/ for output files.")
EOF

# Make R script executable
chmod +x match_predictions.R

# Run R script
Rscript match_predictions.R
