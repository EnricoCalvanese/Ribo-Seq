#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=analyze_counts
#SBATCH --output=analyze_counts_%j.out
#SBATCH --error=analyze_counts_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script for analysis
cat << 'EOF' > analyze_counts.R
# This script analyzes the count data from HTSeq output files
# It checks for non-zero counts and calculates basic statistics

# Function to analyze a single count file
analyze_counts <- function(filename) {
    # Read the count data
    counts <- read.table(filename, stringsAsFactors=FALSE)
    
    # Extract counts (second column)
    count_values <- counts[,2]
    
    # Calculate statistics
    total_features <- length(count_values)
    nonzero_counts <- sum(count_values > 0)
    max_count <- max(count_values)
    mean_count <- mean(count_values)
    
    # Return results
    list(
        filename = basename(filename),
        total_features = total_features,
        nonzero_features = nonzero_counts,
        percent_nonzero = round(100 * nonzero_counts/total_features, 2),
        max_count = max_count,
        mean_count = round(mean_count, 2)
    )
}

# Get list of uORF count files
count_files <- list.files("results/counts", 
                         pattern="*_uorf_counts.txt", 
                         full.names=TRUE)

# Analyze each file
results <- lapply(count_files, analyze_counts)

# Create summary report
report <- paste0(
    "uORF Count Data Analysis\n",
    "=====================\n\n"
)

for(res in results) {
    report <- paste0(
        report,
        "File: ", res$filename, "\n",
        "------------------------\n",
        "Total features: ", res$total_features, "\n",
        "Features with non-zero counts: ", res$nonzero_features, "\n",
        "Percent features with counts: ", res$percent_nonzero, "%\n",
        "Maximum count: ", res$max_count, "\n",
        "Mean count: ", res$mean_count, "\n\n"
    )
}

# Write report to file
writeLines(report, "results/counts/count_analysis_report.txt")

# Also display on console
cat(report)
EOF

# Run R script
Rscript analyze_counts.R
