#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=filter_translated_genes
#SBATCH --output=rpkm_filter_%j.out
#SBATCH --error=rpkm_filter_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script
cat << 'EOF' > filter_by_rpkm.R
# Load required libraries
library(rtracklayer)

print("Starting RPKM-based filtering...")

# Function to safely read a file
safe_read <- function(filepath, ...) {
    if (!file.exists(filepath)) {
        stop(paste("File not found:", filepath))
    }
    read.table(filepath, ...)
}

# Function to calculate RPKM
calculate_rpkm <- function(counts, lengths, total_reads) {
    # RPKM = (read_count * 10^9) / (gene_length * total_reads)
    rpkm <- (counts * 10^9) / (lengths * total_reads)
    return(rpkm)
}

# Read GTF files to get feature lengths
print("Loading GTF files...")
morf_gtf <- safe_read("results/counts/morfs.gtf", 
                      header=FALSE, sep="\t", quote="")
uorf_gtf <- safe_read("results/counts/predicted_uorfs.gtf",
                      header=FALSE, sep="\t", quote="")

# Calculate lengths
print("Calculating feature lengths...")
# For mORFs (one per gene)
morf_lengths <- data.frame(
    gene_id = sub('gene_id "([^"]+)".*', "\\1", morf_gtf$V9),
    length = morf_gtf$V5 - morf_gtf$V4 + 1
)

# For uORFs (multiple possible per gene)
uorf_lengths <- data.frame(
    gene_id = sub('gene_id "([^"]+)".*', "\\1", uorf_gtf$V9),
    length = uorf_gtf$V5 - uorf_gtf$V4 + 1
)

# Process each sample
samples <- c("LZT103-1", "LZT103-2", "LZT104-1", "LZT104-2")
results_summary <- data.frame()

for(sample in samples) {
    print(paste("\nProcessing sample:", sample))
    
    # Read count files
    morf_counts <- safe_read(paste0("results/counts/", sample, "_morf_counts.txt"))
    uorf_counts <- safe_read(paste0("results/counts/", sample, "_uorf_counts.txt"))
    
    # Calculate total reads for RPKM
    total_reads <- sum(morf_counts$V2) + sum(uorf_counts$V2)
    
    # Calculate RPKM for mORFs
    print("Calculating mORF RPKMs...")
    morf_data <- data.frame(
        gene_id = morf_counts$V1,
        counts = morf_counts$V2,
        length = morf_lengths$length[match(morf_counts$V1, morf_lengths$gene_id)]
    )
    morf_data$rpkm <- calculate_rpkm(morf_data$counts, morf_data$length, total_reads)
    
    # Filter for genes with RPKM >= 1 in mORF
    genes_pass_rpkm <- morf_data$gene_id[morf_data$rpkm >= 1]
    print(paste("Genes with mORF RPKM >= 1:", length(genes_pass_rpkm)))
    
    # Among these, find genes with uORF reads
    uorf_data <- data.frame(
        gene_id = sub("_.*", "", uorf_counts$V1),
        counts = uorf_counts$V2
    )
    # Sum counts for all uORFs of the same gene
    uorf_by_gene <- aggregate(counts ~ gene_id, data=uorf_data, sum)
    
    # Filter for genes that pass both criteria
    genes_pass_both <- intersect(
        genes_pass_rpkm,
        uorf_by_gene$gene_id[uorf_by_gene$counts > 0]
    )
    
    print(paste("Genes passing both filters:", length(genes_pass_both)))
    
    # Save filtered gene list
    write.table(genes_pass_both,
                file=paste0("results/counts/", sample, "_filtered_genes.txt"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    # Update results summary
    results_summary <- rbind(results_summary, data.frame(
        sample = sample,
        total_genes = nrow(morf_data),
        genes_rpkm = length(genes_pass_rpkm),
        genes_with_uorf_reads = sum(uorf_by_gene$counts > 0),
        genes_pass_both = length(genes_pass_both)
    ))
}

# Create comprehensive summary report
print("\nCreating summary report...")
summary_report <- paste0(
    "RPKM and Read Coverage Filtering Summary\n",
    "=====================================\n\n"
)

for(i in 1:nrow(results_summary)) {
    summary_report <- paste0(
        summary_report,
        "Sample: ", results_summary$sample[i], "\n",
        "- Total genes analyzed: ", results_summary$total_genes[i], "\n",
        "- Genes with mORF RPKM >= 1: ", results_summary$genes_rpkm[i], "\n",
        "- Genes with uORF reads: ", results_summary$genes_with_uorf_reads[i], "\n",
        "- Genes passing both filters: ", results_summary$genes_pass_both[i], "\n\n"
    )
}

# Save summary report
writeLines(summary_report, "results/counts/filtering_summary.txt")
print("Analysis complete. Check results/counts/filtering_summary.txt for results.")

# Print summary to console
cat(summary_report)
EOF

# Run R script
Rscript filter_by_rpkm.R
