#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=test_uorf
#SBATCH --output=test_uorf_%j.out
#SBATCH --error=test_uorf_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create the R script
cat << 'EOF' > test_prediction.R
# Load required libraries
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)

# Test gene
test_gene <- "AT1G01030"

# Load reference uORF data
ref_uorfs <- read.table("results/xu2017uORFs.txt", 
                       header=TRUE, 
                       stringsAsFactors=FALSE)

# Filter reference data for our test gene
test_ref_uorfs <- ref_uorfs[grep(test_gene, ref_uorfs$uORF_ID),]
print("Reference uORFs:")
print(test_ref_uorfs)

# Load UTR sequences
print("\nLoading UTR sequences...")
utr_sequences <- readDNAStringSet("results/reference_analysis/filtered_utrs.fa")

# Find the sequence for our test gene
test_seq_idx <- grep(paste0("transcript:", test_gene, "\\."), names(utr_sequences))
if(length(test_seq_idx) == 0) {
    stop("Test gene sequence not found!")
}
test_sequence <- utr_sequences[test_seq_idx]

print(paste("\nAnalyzing UTR sequence for", test_gene))
print(paste("UTR length:", width(test_sequence)))

# Predict uORFs
print("\nPredicting uORFs...")
uorf_predictions <- predORF(test_sequence,
                          n="all",
                          mode="orf",
                          longest_disjoint=TRUE,
                          strand="sense")

# Process predictions
if(length(uorf_predictions) > 0 && length(uorf_predictions[[1]]) > 0) {
    # Extract lengths
    pred_lengths <- width(ranges(uorf_predictions[[1]]))
    
    # Create output data frame
    results <- data.frame(
        uORF_ID = paste0(test_gene, "_", seq(0, length(pred_lengths)-1)),
        Length_uORF = pred_lengths
    )
    
    print("\nPredicted uORFs:")
    print(results)
    
    # Save results
    write.table(results, 
                "results/test_predictions.txt", 
                row.names=FALSE, 
                quote=FALSE, 
                sep="\t")
    
    # Compare with reference
    print("\nComparison with reference:")
    print("Reference lengths:")
    print(sort(test_ref_uorfs$Length_uORF))
    print("Predicted lengths:")
    print(sort(results$Length_uORF))
} else {
    print("No uORFs predicted!")
}
EOF

# Make R script executable
chmod +x test_prediction.R

# Run R script
Rscript test_prediction.R
