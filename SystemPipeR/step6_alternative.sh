#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=compare_counting
#SBATCH --output=compare_counting_%j.out
#SBATCH --error=compare_counting_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required modules
module load r
module load python

# Create directories for different approaches
mkdir -p results/counts/approach1
mkdir -p results/counts/approach2

# Create R script to generate modified GTF
cat << 'EOF' > modify_gtf.R
# This script creates two versions of our GTF file:
# 1. Original version with uORF features
# 2. Modified version with CDS features

print("Reading original GTF file...")
gtf_data <- read.table("results/counts/predicted_uorfs.gtf", 
                      sep="\t", 
                      quote="", 
                      stringsAsFactors=FALSE)

# Save original version (approach 1)
print("Saving GTF for approach 1 (uORF features)...")
write.table(gtf_data, 
            "results/counts/approach1/uorf_features.gtf",
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE)

# Modify feature type to CDS and save (approach 2)
print("Creating and saving GTF for approach 2 (CDS features)...")
gtf_data[,3] <- "CDS"  # Change feature type from uORF to CDS
write.table(gtf_data, 
            "results/counts/approach2/cds_features.gtf",
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE)

print("GTF modifications complete.")
EOF

# Run R script to create modified GTF files
echo "Creating modified GTF files..."
Rscript modify_gtf.R

# Function to run htseq-count with detailed logging
run_htseq_count() {
    local bam=$1
    local gtf=$2
    local output=$3
    local feature_type=$4
    local log_file="${output%.txt}_log.txt"
    
    echo "Starting HTSeq counting for ${bam} with feature type ${feature_type}"
    echo "Parameters used:" > "${log_file}"
    echo "- BAM file: ${bam}" >> "${log_file}"
    echo "- GTF file: ${gtf}" >> "${log_file}"
    echo "- Feature type: ${feature_type}" >> "${log_file}"
    echo "- Output file: ${output}" >> "${log_file}"
    echo "-----------------" >> "${log_file}"
    
    # Run HTSeq count with verbose output
    htseq-count \
        -f bam \
        -r name \
        -s yes \
        -t ${feature_type} \
        -i gene_id \
        --additional-attr=transcript_id \
        -m union \
        --nonunique all \
        --verbose \
        ${bam} ${gtf} 2>> "${log_file}" > "${output}"
    
    # Add basic statistics to log file
    echo "-----------------" >> "${log_file}"
    echo "Count Statistics:" >> "${log_file}"
    total_features=$(grep -v "__" "${output}" | wc -l)
    nonzero_features=$(grep -v "__" "${output}" | awk '$2 > 0' | wc -l)
    max_count=$(grep -v "__" "${output}" | awk 'BEGIN{max=0} {if($2>max) max=$2} END{print max}')
    echo "Total features: ${total_features}" >> "${log_file}"
    echo "Features with counts > 0: ${nonzero_features}" >> "${log_file}"
    echo "Maximum count: ${max_count}" >> "${log_file}"
}

# Process each BAM file with both approaches
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
for sample in LZT103-1 LZT103-2 LZT104-1 LZT104-2; do
    echo "Processing ${sample}..."
    
    # Approach 1: Using uORF feature type
    echo "Running Approach 1 (uORF features)..."
    run_htseq_count \
        "${BAM_DIR}/${sample}_uniq_sort.bam" \
        "results/counts/approach1/uorf_features.gtf" \
        "results/counts/approach1/${sample}_counts.txt" \
        "uORF"
    
    # Approach 2: Using CDS feature type
    echo "Running Approach 2 (CDS features)..."
    run_htseq_count \
        "${BAM_DIR}/${sample}_uniq_sort.bam" \
        "results/counts/approach2/cds_features.gtf" \
        "results/counts/approach2/${sample}_counts.txt" \
        "CDS"
done

# Create summary report comparing both approaches
echo "Creating comparison report..."
cat << 'EOF' > compare_results.R
# Load required libraries
library(data.table)

# Function to analyze count file
analyze_counts <- function(file) {
    counts <- fread(file)
    total <- nrow(counts[!grepl("^__", V1)])
    nonzero <- nrow(counts[!grepl("^__", V1) & V2 > 0])
    max_count <- max(counts[!grepl("^__", V1)]$V2)
    mean_count <- mean(counts[!grepl("^__", V1)]$V2)
    
    list(
        total = total,
        nonzero = nonzero,
        percent = round(100 * nonzero/total, 2),
        max = max_count,
        mean = round(mean_count, 2)
    )
}

# Analyze both approaches
samples <- c("LZT103-1", "LZT103-2", "LZT104-1", "LZT104-2")
results <- data.frame(
    sample = character(),
    approach = character(),
    total_features = numeric(),
    nonzero_features = numeric(),
    percent_nonzero = numeric(),
    max_count = numeric(),
    mean_count = numeric(),
    stringsAsFactors = FALSE
)

for(sample in samples) {
    # Analyze Approach 1
    a1 <- analyze_counts(paste0("results/counts/approach1/", sample, "_counts.txt"))
    results <- rbind(results, data.frame(
        sample = sample,
        approach = "uORF",
        total_features = a1$total,
        nonzero_features = a1$nonzero,
        percent_nonzero = a1$percent,
        max_count = a1$max,
        mean_count = a1$mean
    ))
    
    # Analyze Approach 2
    a2 <- analyze_counts(paste0("results/counts/approach2/", sample, "_counts.txt"))
    results <- rbind(results, data.frame(
        sample = sample,
        approach = "CDS",
        total_features = a2$total,
        nonzero_features = a2$nonzero,
        percent_nonzero = a2$percent,
        max_count = a2$max,
        mean_count = a2$mean
    ))
}

# Create detailed report
report <- paste0(
    "Comparison of HTSeq Counting Approaches\n",
    "===================================\n\n"
)

for(sample in samples) {
    report <- paste0(
        report,
        "Sample: ", sample, "\n",
        "----------------\n\n"
    )
    
    for(approach in c("uORF", "CDS")) {
        r <- results[results$sample == sample & results$approach == approach,]
        report <- paste0(
            report,
            approach, " approach:\n",
            "- Total features: ", r$total_features, "\n",
            "- Features with counts: ", r$nonzero_features, " (", r$percent_nonzero, "%)\n",
            "- Maximum count: ", r$max_count, "\n",
            "- Mean count: ", r$mean_count, "\n\n"
        )
    }
    report <- paste0(report, "\n")
}

# Add overall insights
report <- paste0(
    report,
    "\nKey Insights:\n",
    "------------\n",
    "1. Feature type impact: ",
    ifelse(all(results$nonzero_features[results$approach == "uORF"] == 
              results$nonzero_features[results$approach == "CDS"]),
           "No difference between approaches\n",
           "Approaches show different results\n"),
    "2. Consistency across samples: ",
    ifelse(length(unique(results$nonzero_features)) == 1,
           "All samples show identical feature counts\n",
           "Samples show varying feature counts\n"),
    "\nRecommendation:\n",
    "--------------\n",
    ifelse(mean(results$nonzero_features[results$approach == "CDS"]) >
           mean(results$nonzero_features[results$approach == "uORF"]),
           "The CDS approach appears more effective\n",
           "The uORF approach appears more effective\n")
)

writeLines(report, "results/counts/approach_comparison.txt")
EOF

# Run comparison script
Rscript compare_results.R

echo "Analysis complete. Check results/counts/approach_comparison.txt for detailed comparison."
