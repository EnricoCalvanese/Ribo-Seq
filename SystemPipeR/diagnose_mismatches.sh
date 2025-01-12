#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=0:30:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=diagnose_uorfs
#SBATCH --output=diagnose_uorfs_%j.out
#SBATCH --error=diagnose_uorfs_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create the R script
cat << 'EOF' > diagnose_mismatches.R
# Load required libraries
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)

print("Starting mismatch diagnosis...")

# Load all necessary data
print("Loading data...")
ref_uorfs <- read.table("results/xu2017uORFs.txt", 
                       header=TRUE, 
                       stringsAsFactors=FALSE)
match_stats <- read.table("results/reference_analysis/matching_statistics.txt",
                         header=TRUE, 
                         stringsAsFactors=FALSE)
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")
matched_results <- read.table("results/reference_analysis/matched_uorfs.txt",
                            header=TRUE, 
                            stringsAsFactors=FALSE)

# Create lookup tables
ref_by_gene <- split(ref_uorfs, sub("_.*", "", ref_uorfs$uORF_ID))
matched_by_gene <- split(matched_results, sub("_.*", "", matched_results$uORF_ID))

# Function to analyze mismatches for a gene
analyze_gene_mismatch <- function(gene_id, ref_data, matched_data, utr_seq) {
    ref_lengths <- ref_data$Length_uORF
    matched_lengths <- if(gene_id %in% names(matched_data)) {
        matched_data[[gene_id]]$Length_uORF
    } else {
        numeric(0)
    }
    
    # Get UTR sequence if available
    seq_idx <- grep(paste0("transcript:", gene_id, "\\."), names(utr_seq))
    utr_info <- if(length(seq_idx) > 0) {
        paste("UTR length:", width(utr_seq[seq_idx]))
    } else {
        "UTR sequence not found"
    }
    
    return(list(
        gene = gene_id,
        ref_uorfs = length(ref_lengths),
        ref_lengths = paste(ref_lengths, collapse=", "),
        matched_uorfs = length(matched_lengths),
        matched_lengths = paste(matched_lengths, collapse=", "),
        utr_info = utr_info
    ))
}

# 1. Analyze genes with no matches
print("\nAnalyzing genes with no matches...")
no_match_genes <- match_stats$gene_id[match_stats$matched_uorfs == 0]
no_match_analysis <- lapply(no_match_genes, function(gene) {
    analyze_gene_mismatch(gene, ref_by_gene[[gene]], 
                         matched_by_gene, utr_sequences)
})

# 2. Analyze genes with partial matches
print("Analyzing genes with partial matches...")
partial_match_genes <- match_stats$gene_id[
    match_stats$matched_uorfs > 0 & 
    match_stats$matched_uorfs < match_stats$ref_uorfs
]
partial_match_analysis <- lapply(partial_match_genes, function(gene) {
    analyze_gene_mismatch(gene, ref_by_gene[[gene]], 
                         matched_by_gene, utr_sequences)
})

# 3. Find genes where we predicted more uORFs than reference
print("Analyzing genes with extra predictions...")
extra_pred_genes <- match_stats$gene_id[
    match_stats$matched_uorfs > match_stats$ref_uorfs
]
extra_pred_analysis <- lapply(extra_pred_genes, function(gene) {
    analyze_gene_mismatch(gene, ref_by_gene[[gene]], 
                         matched_by_gene, utr_sequences)
})

# Convert analysis results to data frames
no_match_df <- do.call(rbind, lapply(no_match_analysis, as.data.frame))
partial_match_df <- do.call(rbind, lapply(partial_match_analysis, as.data.frame))
extra_pred_df <- do.call(rbind, lapply(extra_pred_analysis, as.data.frame))

# Generate detailed report
mismatch_report <- paste0(
    "uORF Prediction Mismatch Analysis\n",
    "==============================\n\n",
    "1. Genes with No Matches (", length(no_match_genes), " genes)\n",
    "----------------------------------------\n",
    capture.output(head(no_match_df, 10)), "\n",
    "\n2. Genes with Partial Matches (", length(partial_match_genes), " genes)\n",
    "----------------------------------------\n",
    capture.output(head(partial_match_df, 10)), "\n",
    "\n3. Genes with Extra Predictions (", length(extra_pred_genes), " genes)\n",
    "----------------------------------------\n",
    capture.output(head(extra_pred_df, 10)), "\n",
    "\nCommon Issues Found:\n",
    "1. UTR Length Issues:\n",
    "   - Number of genes with missing UTR: ",
    sum(grepl("not found", c(no_match_df$utr_info, 
                            partial_match_df$utr_info, 
                            extra_pred_df$utr_info))), "\n",
    "\n2. Prediction Pattern Issues:\n",
    "   - Genes with completely missing predictions: ", nrow(no_match_df), "\n",
    "   - Genes with partial predictions: ", nrow(partial_match_df), "\n",
    "   - Genes with extra predictions: ", nrow(extra_pred_df), "\n"
)

# Save detailed results
print("Saving detailed analysis...")
write.table(no_match_df,
            "results/reference_analysis/no_match_analysis.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(partial_match_df,
            "results/reference_analysis/partial_match_analysis.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
write.table(extra_pred_df,
            "results/reference_analysis/extra_predictions_analysis.txt",
            row.names=FALSE, quote=FALSE, sep="\t")
writeLines(mismatch_report,
           "results/reference_analysis/mismatch_analysis_report.txt")

print("Analysis complete. Check results/reference_analysis/ for detailed reports.")
EOF

# Make R script executable
chmod +x diagnose_mismatches.R

# Run R script
Rscript diagnose_mismatches.R
