#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --job-name=validate_filter
#SBATCH --output=validate_filter_%j.out
#SBATCH --error=validate_filter_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script
cat << 'EOF' > validate_filtering.R
# Create validation directory
dir.create("results/counts/validation", showWarnings=FALSE, recursive=TRUE)

# Function to manually calculate RPKM
manual_rpkm <- function(counts, length, total_reads) {
    (counts * 10^9) / (length * total_reads)
}

# Load complete analysis files
samples <- c("LZT103-1", "LZT103-2", "LZT104-1", "LZT104-2")
all_data <- list()
filtered_genes <- list()

for(sample in samples) {
    # Load complete analysis data
    all_data[[sample]] <- read.table(
        paste0("results/counts/diagnostics/", sample, "_complete_analysis.txt"),
        header=TRUE, stringsAsFactors=FALSE
    )
    
    # Load filtered genes
    filtered_genes[[sample]] <- readLines(
        paste0("results/counts/", sample, "_filtered_genes.txt")
    )
}

# Open validation report file
sink("results/counts/validation/validation_report.txt")

cat("Validation Report for Filtering Results\n")
cat("=====================================\n\n")

# 1. Manual RPKM calculation check
cat("1. RPKM Calculation Verification\n")
cat("--------------------------------\n")
for(sample in samples) {
    cat("\nSample:", sample, "\n")
    
    # Take top 5 genes by RPKM for verification
    top_genes <- head(all_data[[sample]][order(-all_data[[sample]]$rpkm), ], 5)
    
    for(i in 1:nrow(top_genes)) {
        manual <- manual_rpkm(
            top_genes$counts[i],
            top_genes$length[i],
            sum(all_data[[sample]]$counts)
        )
        
        cat(sprintf(
            "Gene: %s\n  Counts: %d\n  Length: %d\n  Calculated RPKM: %.2f\n  Reported RPKM: %.2f\n  Difference: %.2f%%\n",
            top_genes$gene_id[i],
            top_genes$counts[i],
            top_genes$length[i],
            manual,
            top_genes$rpkm[i],
            abs(manual - top_genes$rpkm[i])/manual * 100
        ))
    }
}

# 2. Verify filtering criteria
cat("\n2. Filtering Criteria Verification\n")
cat("--------------------------------\n")
for(sample in samples) {
    cat("\nSample:", sample, "\n")
    
    # Get data for filtered genes
    filtered_data <- all_data[[sample]][all_data[[sample]]$gene_id %in% filtered_genes[[sample]], ]
    
    # Check RPKM criterion
    rpkm_fail <- filtered_data$rpkm < 1
    if(any(rpkm_fail)) {
        cat("WARNING: Found", sum(rpkm_fail), "genes with RPKM < 1\n")
        print(filtered_data[rpkm_fail, ])
    }
    
    # Check uORF reads criterion
    uorf_fail <- filtered_data$uorf_counts == 0 | is.na(filtered_data$uorf_counts)
    if(any(uorf_fail)) {
        cat("WARNING: Found", sum(uorf_fail), "genes without uORF reads\n")
        print(filtered_data[uorf_fail, ])
    }
    
    # Print summary stats for filtered genes
    cat(sprintf(
        "Summary for filtered genes:\n  RPKM range: %.2f - %.2f\n  uORF reads range: %d - %d\n",
        min(filtered_data$rpkm),
        max(filtered_data$rpkm),
        min(filtered_data$uorf_counts),
        max(filtered_data$uorf_counts)
    ))
}

# 3. Check replicate overlap
cat("\n3. Replicate Overlap Analysis\n")
cat("--------------------------------\n")

# WT replicates
wt_overlap <- length(intersect(
    filtered_genes[["LZT103-1"]],
    filtered_genes[["LZT103-2"]]
))
cat(sprintf(
    "\nWT replicates overlap:\n  Rep1: %d genes\n  Rep2: %d genes\n  Overlap: %d genes (%.1f%%)\n",
    length(filtered_genes[["LZT103-1"]]),
    length(filtered_genes[["LZT103-2"]]),
    wt_overlap,
    100 * wt_overlap/min(length(filtered_genes[["LZT103-1"]]), 
                        length(filtered_genes[["LZT103-2"]]))
))

# imb2 replicates
imb2_overlap <- length(intersect(
    filtered_genes[["LZT104-1"]],
    filtered_genes[["LZT104-2"]]
))
cat(sprintf(
    "\nimb2 replicates overlap:\n  Rep1: %d genes\n  Rep2: %d genes\n  Overlap: %d genes (%.1f%%)\n",
    length(filtered_genes[["LZT104-1"]]),
    length(filtered_genes[["LZT104-2"]]),
    imb2_overlap,
    100 * imb2_overlap/min(length(filtered_genes[["LZT104-1"]]), 
                          length(filtered_genes[["LZT104-2"]]))
))

# Create 4-way Venn diagram data
cat("\n4. Four-way comparison:\n")
all_genes <- unique(unlist(filtered_genes))
presence_matrix <- matrix(0, nrow=length(all_genes), ncol=4)
colnames(presence_matrix) <- samples
rownames(presence_matrix) <- all_genes

for(i in seq_along(samples)) {
    presence_matrix[,i] <- all_genes %in% filtered_genes[[samples[i]]]
}

# Count genes in each combination
combinations <- apply(presence_matrix, 1, paste, collapse="")
combo_counts <- table(combinations)
cat("\nGene presence patterns (1=present, 0=absent):\n")
print(combo_counts)

sink()

print("Validation complete. Check results/counts/validation/validation_report.txt")
EOF

# Run R script
Rscript validate_filtering.R
