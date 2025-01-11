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
library(GenomicFeatures)
library(Biostrings)

# Function to ensure directory exists and is writable
ensure_directory <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
        print(paste("Created directory:", dir_path))
    }
    if (!file.access(dir_path, mode = 2) == 0) {
        stop(paste("Directory", dir_path, "is not writable"))
    }
    print(paste("Directory", dir_path, "is ready for writing"))
}

print("Starting uORF prediction analysis...")

# Ensure output directory exists and is writable
ensure_directory("results/predictions")

# Load our previously filtered and processed UTR sequences
print("Loading filtered 5' UTR sequences...")
utr_sequences <- readDNAStringSet("results/sequences/utr_sequences.fa")
print(paste("Loaded", length(utr_sequences), "UTR sequences"))

# Predict ORFs using predORF
print("Predicting ORFs...")
uorf_predictions <- predORF(utr_sequences, 
                          n="all",           
                          mode="orf",        
                          longest_disjoint=TRUE,  
                          strand="sense")    

print(paste("Found predictions for", length(uorf_predictions), "sequences"))

# Process results to get gene-level information
print("Processing results...")

# Function to extract gene ID from transcript ID
get_gene_id <- function(transcript_id) {
    sub("transcript:([^.]+).*", "\\1", transcript_id)
}

# Safety check before creating data frame
if (length(uorf_predictions) == 0) {
    stop("No ORF predictions found!")
}

# Convert results to data frame with error checking
print("Converting predictions to data frame...")
results_df <- tryCatch({
    data.frame(
        transcript_id = rep(names(uorf_predictions), 
                          sapply(uorf_predictions, length)),
        start = unlist(lapply(uorf_predictions, function(x) start(ranges(x)))),
        end = unlist(lapply(uorf_predictions, function(x) end(ranges(x)))),
        width = unlist(lapply(uorf_predictions, function(x) width(ranges(x))))
    )
}, error = function(e) {
    print(paste("Error creating data frame:", e$message))
    print("Examining uorf_predictions structure:")
    print(str(uorf_predictions))
    stop("Failed to create results data frame")
})

print(paste("Created data frame with", nrow(results_df), "rows"))

# Add gene IDs
results_df$gene_id <- get_gene_id(results_df$transcript_id)

# Filter ORFs
valid_orfs <- results_df$width >= 6  
results_df <- results_df[valid_orfs,]

print(paste("After filtering, data frame has", nrow(results_df), "rows"))

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
    "- Median length: ", round(median(uorf_lengths), 2), " nt\n"
)

print("Saving results...")

# Try saving results with error catching
tryCatch({
    saveRDS(results_df, "results/predictions/uorf_predictions.rds")
    print("Saved uorf_predictions.rds")
    
    writeLines(summary_report, "results/predictions/prediction_summary.txt")
    print("Saved prediction_summary.txt")
    
    write.table(genes_with_uorfs, "results/predictions/genes_with_uorfs.txt", 
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    print("Saved genes_with_uorfs.txt")
}, error = function(e) {
    print(paste("Error saving files:", e$message))
    stop("Failed to save output files")
})

print("Analysis complete. Check results/predictions/ for output files.")
EOF

# Make R script executable
chmod +x predict_uorfs.R

# Run R script
Rscript predict_uorfs.R
