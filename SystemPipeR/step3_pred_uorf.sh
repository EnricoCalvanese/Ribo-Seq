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
# SystemPipeR script for uORF prediction with validation
library(systemPipeR)
library(GenomicFeatures)
library(Biostrings)

# Function to predict uORFs with strict validation
predictUORFs <- function(seq, seq_length) {
    # Define codons
    start_codons <- c("ATG")
    stop_codons <- c("TAA", "TAG", "TGA")
    min_length <- 6  # Strict enforcement
    
    # Find all potential start positions
    start_pos <- sort(unique(unlist(lapply(start_codons, function(codon) {
        start(matchPattern(codon, seq))
    }))))
    
    # Find all stop positions
    stop_pos <- sort(unique(unlist(lapply(stop_codons, function(codon) {
        end(matchPattern(codon, seq))
    }))))
    
    # Initialize list for valid ORFs
    valid_orfs <- list()
    
    # Track used positions to avoid overlapping ORFs in the same frame
    used_positions <- matrix(FALSE, nrow=3, ncol=seq_length)
    
    for(start in start_pos) {
        frame <- (start - 1) %% 3 + 1
        
        # Skip if this position is already used in this frame
        if(used_positions[frame, start]) next
        
        # Find valid stop codons
        valid_stops <- stop_pos[stop_pos > start & (stop_pos - start + 3) %% 3 == 0]
        
        if(length(valid_stops) > 0) {
            first_stop <- min(valid_stops)
            orf_length <- first_stop - start + 3
            
            # Validate ORF
            if(validateORF(start, first_stop, orf_length, seq_length, min_length)) {
                valid_orfs[[length(valid_orfs) + 1]] <- list(
                    start = start,
                    stop = first_stop,
                    length = orf_length,
                    frame = frame
                )
                
                # Mark positions as used in this frame
                used_positions[frame, start:first_stop] <- TRUE
            }
        }
    }
    
    return(valid_orfs)
}

# Function to validate individual ORFs
validateORF <- function(start, stop, length, seq_length, min_length) {
    # Check minimum length
    if(length < min_length) return(FALSE)
    
    # Check if ORF is entirely within UTR
    if(start < 1 || stop > seq_length) return(FALSE)
    
    # Add any additional validation criteria here
    
    return(TRUE)
}

# Main execution
main <- function() {
    print("Starting analysis with improved validation...")
    
    # Create results directory
    dir.create("results/predictions", recursive = TRUE, showWarnings = FALSE)
    
    # Load and process 5' UTR sequences
    print("Loading UTR sequences...")
    utr_sequences <- readDNAStringSet("results/sequences/utr_sequences.fa")
    
    # Initialize counters for detailed statistics
    stats <- list(
        total_utrs = length(utr_sequences),
        utrs_with_uorfs = 0,
        total_uorfs = 0,
        uorfs_by_frame = c(0,0,0),
        lengths = numeric()
    )
    
    # Process sequences and predict uORFs
    print("Predicting uORFs with validation...")
    uorf_predictions <- list()
    
    for(i in seq_along(utr_sequences)) {
        seq_name <- names(utr_sequences)[i]
        seq <- utr_sequences[[i]]
        seq_length <- length(seq)
        
        orfs <- predictUORFs(seq, seq_length)
        
        if(length(orfs) > 0) {
            uorf_predictions[[seq_name]] <- orfs
            stats$utrs_with_uorfs <- stats$utrs_with_uorfs + 1
            stats$total_uorfs <- stats$total_uorfs + length(orfs)
            
            # Collect frame statistics
            frames <- sapply(orfs, function(x) x$frame)
            stats$uorfs_by_frame <- stats$uorfs_by_frame + tabulate(frames, 3)
            
            # Collect length statistics
            stats$lengths <- c(stats$lengths, sapply(orfs, function(x) x$length))
        }
        
        if(i %% 1000 == 0) {
            print(paste0("Processed ", i, " sequences..."))
        }
    }
    
    # Process mORF information
    print("Processing mORF information...")
    txdb <- makeTxDbFromGFF("/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf")
    cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
    morf_count <- length(cds_by_tx)
    
    # Save predictions
    saveRDS(uorf_predictions, "results/predictions/uorf_predictions.rds")
    saveRDS(cds_by_tx, "results/predictions/morf_annotations.rds")
    
    # Create detailed summary report
    summary_report <- paste0(
        "uORF and mORF Analysis Summary (with validation)\n",
        "==========================================\n\n",
        "Input Data Validation:\n",
        "- Total 5' UTR sequences analyzed: ", stats$total_utrs, "\n",
        "- Total mORFs found: ", morf_count, "\n",
        "- UTR to mORF ratio: ", round(stats$total_utrs/morf_count, 2), "\n\n",
        "Prediction Parameters:\n",
        "- Minimum ORF length: 6 nucleotides (strictly enforced)\n",
        "- Start codons considered: ATG\n",
        "- Stop codons considered: TAA, TAG, TGA\n",
        "- Frame overlap prevention: enabled\n",
        "- Position validation: enforced\n\n",
        "uORF Prediction Results:\n",
        "- Number of UTRs with predicted uORFs: ", stats$utrs_with_uorfs, "\n",
        "- Total number of predicted uORFs: ", stats$total_uorfs, "\n",
        "- Average uORFs per UTR (when present): ", 
            round(stats$total_uorfs/stats$utrs_with_uorfs, 2), "\n\n",
        "Frame Distribution:\n",
        "- Frame 1: ", stats$uorfs_by_frame[1], " uORFs\n",
        "- Frame 2: ", stats$uorfs_by_frame[2], " uORFs\n",
        "- Frame 3: ", stats$uorfs_by_frame[3], " uORFs\n\n",
        "Length Distribution of Predicted uORFs:\n",
        "- Minimum length: ", min(stats$lengths), " nt\n",
        "- Maximum length: ", max(stats$lengths), " nt\n",
        "- Mean length: ", round(mean(stats$lengths), 2), " nt\n",
        "- Median length: ", round(median(stats$lengths), 2), " nt\n\n",
        "Validation Notes:\n",
        "- All uORFs are confirmed to be entirely within UTR boundaries\n",
        "- Frame overlaps have been prevented\n",
        "- Minimum length of 6 nt is strictly enforced\n"
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
