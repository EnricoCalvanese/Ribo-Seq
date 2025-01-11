#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
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
# SystemPipeR script for uORF prediction
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)

# Function to count entries in FASTA file
countFastaEntries <- function(fasta_file) {
    fasta <- readDNAStringSet(fasta_file)
    return(length(fasta))
}

# Function to predict uORFs
predictUORFs <- function(seq) {
    start_codons <- c("ATG")  
    stop_codons <- c("TAA", "TAG", "TGA")
    min_length <- 6
    
    start_pos <- sapply(start_codons, function(codon) {
        start(matchPattern(codon, seq))
    })
    start_pos <- sort(unique(unlist(start_pos)))
    
    stop_pos <- sapply(stop_codons, function(codon) {
        end(matchPattern(codon, seq))
    })
    stop_pos <- sort(unique(unlist(stop_pos)))
    
    orfs <- list()
    for(start in start_pos) {
        valid_stops <- stop_pos[stop_pos > start & (stop_pos - start + 3) %% 3 == 0]
        if(length(valid_stops) > 0) {
            first_stop <- min(valid_stops)
            length <- first_stop - start + 3
            if(length >= min_length) {
                orfs[[length(orfs) + 1]] <- c(start, first_stop)
            }
        }
    }
    return(orfs)
}

# Main execution
main <- function() {
    print("Starting analysis...")
    
    # Create results directory
    dir.create("results/predictions", recursive = TRUE, showWarnings = FALSE)
    
    # Count UTRs in our input file
    print("Counting UTR sequences...")
    utr_count <- countFastaEntries("results/sequences/utr_sequences.fa")
    print(paste("Number of UTR sequences:", utr_count))
    
    # Load 5' UTR sequences
    print("Loading UTR sequences...")
    utr_sequences <- readDNAStringSet("results/sequences/utr_sequences.fa")
    
    # Process sequences and predict uORFs
    print("Predicting uORFs...")
    uorf_predictions <- list()
    for(i in seq_along(utr_sequences)) {
        seq_name <- names(utr_sequences)[i]
        seq <- utr_sequences[[i]]
        orfs <- predictUORFs(seq)
        if(length(orfs) > 0) {
            uorf_predictions[[seq_name]] <- orfs
        }
        if(i %% 1000 == 0) {
            print(paste0("Processed ", i, " sequences..."))
        }
    }
    
    # Extract mORF information
    print("Processing mORF information...")
    txdb <- makeTxDbFromGFF("/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf")
    cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
    
    # Validation: Count mORFs
    morf_count <- length(cds_by_tx)
    print(paste("Number of mORFs found:", morf_count))
    
    # Generate statistics
    num_uorf_genes <- length(unique(sub("transcript:([^.]+).*", "\\1", names(uorf_predictions))))
    total_uorfs <- sum(sapply(uorf_predictions, length))
    
    # Save predictions
    saveRDS(uorf_predictions, "results/predictions/uorf_predictions.rds")
    saveRDS(cds_by_tx, "results/predictions/morf_annotations.rds")
    
    # Calculate length distribution
    uorf_lengths <- unlist(lapply(uorf_predictions, function(orfs) {
        sapply(orfs, function(orf) orf[2] - orf[1] + 1)
    }))
    
    # Create comprehensive summary report
    summary_report <- paste0(
        "uORF and mORF Analysis Summary\n",
        "===========================\n\n",
        "Input Data Validation:\n",
        "- Total 5' UTR sequences analyzed: ", utr_count, "\n",
        "- Total mORFs found: ", morf_count, "\n",
        "- UTR to mORF ratio: ", round(utr_count/morf_count, 2), "\n\n",
        "Prediction Parameters:\n",
        "- Minimum ORF length: 6 nucleotides\n",
        "- Start codons considered: ATG\n",
        "- Stop codons considered: TAA, TAG, TGA\n\n",
        "uORF Prediction Results:\n",
        "- Number of genes with predicted uORFs: ", num_uorf_genes, "\n",
        "- Total number of predicted uORFs: ", total_uorfs, "\n",
        "- Average uORFs per gene: ", round(total_uorfs/num_uorf_genes, 2), "\n\n",
        "Length Distribution of Predicted uORFs:\n",
        "- Minimum length: ", min(uorf_lengths), " nt\n",
        "- Maximum length: ", max(uorf_lengths), " nt\n",
        "- Mean length: ", round(mean(uorf_lengths), 2), " nt\n",
        "- Median length: ", round(median(uorf_lengths), 2), " nt\n\n",
        "Validation Notes:\n",
        "- Each transcript should have one mORF\n",
        "- The number of mORFs should approximately match the number of UTRs\n",
        "- Significant deviation might indicate annotation issues\n"
    )
    
    writeLines(summary_report, "results/predictions/prediction_summary.txt")
    print("Analysis complete. Check results/predictions/ for output files.")
}

# Execute main function
main()
EOF

# Make R script executable
chmod +x predict_uorfs.R

# Run R script
Rscript predict_uorfs.R
