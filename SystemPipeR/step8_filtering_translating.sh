#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=0:30:00
#SBATCH --job-name=rpkm_filter_diag
#SBATCH --output=rpkm_filter_diag_%j.out
#SBATCH --error=rpkm_filter_diag_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script
cat << 'EOF' > filter_by_rpkm_diagnostic.R
# Load required libraries
library(rtracklayer)

# Create directory for diagnostics
dir.create("results/counts/diagnostics", showWarnings=FALSE, recursive=TRUE)

# Function to write diagnostic info
write_diagnostic <- function(data, filename, sample="", note="") {
    sink(file.path("results/counts/diagnostics", filename), append=TRUE)
    cat("\n=================================\n")
    cat(sample, "\n")
    cat(note, "\n")
    cat("=================================\n")
    print(data)
    sink()
}

# Function to calculate RPKM
calculate_rpkm <- function(counts, lengths, total_reads) {
    rpkm <- (counts * 10^9) / (lengths * total_reads)
    return(rpkm)
}

print("Loading GTF files and calculating lengths...")

# Read GTF files
morf_gtf <- read.table("results/counts/morfs.gtf", header=FALSE, sep="\t", quote="")
write_diagnostic(head(morf_gtf), "gtf_content.txt", note="First few lines of mORF GTF:")

# Calculate mORF lengths
morf_lengths <- data.frame(
    gene_id = sub('gene_id "([^"]+)".*', "\\1", morf_gtf$V9),
    length = morf_gtf$V5 - morf_gtf$V4 + 1
)
morf_lengths$gene_id <- sub("gene:", "", morf_lengths$gene_id)

write_diagnostic(summary(morf_lengths$length), "length_stats.txt", 
                note="mORF length statistics:")

# Process each sample
samples <- c("LZT103-1", "LZT103-2", "LZT104-1", "LZT104-2")

for(sample in samples) {
    print(paste("\nProcessing sample:", sample))
    
    # Read count files
    morf_counts <- read.table(paste0("results/counts/", sample, "_morf_counts.txt"))
    uorf_counts <- read.table(paste0("results/counts/", sample, "_uorf_counts.txt"))
    
    write_diagnostic(head(morf_counts), "raw_counts.txt", sample,
                    "First few mORF counts before processing:")
    
    # Clean gene IDs
    morf_counts$V1 <- sub("gene:", "", morf_counts$V1)
    
    # Diagnostic: Check for HTSeq categories
    htseq_cats <- morf_counts[grep("^__", morf_counts$V1), ]
    write_diagnostic(htseq_cats, "htseq_stats.txt", sample,
                    "HTSeq statistics categories:")
    
    # Calculate total mapped reads
    total_reads <- sum(morf_counts$V2[!grepl("^__", morf_counts$V1)]) + 
                  sum(uorf_counts$V2[!grepl("^__", uorf_counts$V1)])
    write_diagnostic(total_reads, "total_reads.txt", sample,
                    "Total mapped reads:")
    
    # Create mORF data frame with all information
    morf_data <- data.frame(
        gene_id = morf_counts$V1[!grepl("^__", morf_counts$V1)],
        counts = morf_counts$V2[!grepl("^__", morf_counts$V1)]
    )
    morf_data$length <- morf_lengths$length[match(morf_data$gene_id, morf_lengths$gene_id)]
    
    # Check for any missing lengths
    missing_lengths <- is.na(morf_data$length)
    write_diagnostic(sum(missing_lengths), "missing_lengths.txt", sample,
                    "Number of genes with missing length information:")
    if(sum(missing_lengths) > 0) {
        write_diagnostic(head(morf_data[missing_lengths,]), "missing_lengths.txt", sample,
                        "Examples of genes with missing lengths:")
    }
    
    # Calculate RPKM
    morf_data$rpkm <- calculate_rpkm(morf_data$counts, morf_data$length, total_reads)
    
    # Detailed RPKM diagnostics
    write_diagnostic(summary(morf_data$rpkm), "rpkm_stats.txt", sample,
                    "RPKM summary statistics:")
    write_diagnostic(sum(morf_data$rpkm >= 1), "rpkm_filter.txt", sample,
                    "Number of genes with RPKM >= 1:")
    write_diagnostic(head(morf_data[order(-morf_data$rpkm),]), "top_rpkm.txt", sample,
                    "Top RPKM values:")
    write_diagnostic(head(morf_data[morf_data$rpkm < 1,]), "low_rpkm.txt", sample,
                    "Examples of genes with RPKM < 1:")
    
    # Filter for RPKM >= 1
    genes_pass_rpkm <- morf_data$gene_id[morf_data$rpkm >= 1]
    
    # Process uORF counts
    uorf_data <- data.frame(
        gene_id = sub("_.*", "", uorf_counts$V1[!grepl("^__", uorf_counts$V1)]),
        counts = uorf_counts$V2[!grepl("^__", uorf_counts$V1)]
    )
    
    # Diagnostic for uORF counts
    write_diagnostic(table(uorf_data$counts > 0), "uorf_read_dist.txt", sample,
                    "Distribution of genes with/without uORF reads:")
    
    # Sum counts for all uORFs of the same gene
    uorf_by_gene <- aggregate(counts ~ gene_id, data=uorf_data, sum)
    
    # Diagnostic for gene ID matching
    write_diagnostic(head(genes_pass_rpkm), "gene_ids.txt", sample,
                    "Sample of gene IDs passing RPKM filter:")
    write_diagnostic(head(uorf_by_gene$gene_id[uorf_by_gene$counts > 0]), "gene_ids.txt", sample,
                    "Sample of gene IDs with uORF reads:")
    
    # Check gene ID overlap
    genes_in_both <- intersect(genes_pass_rpkm, uorf_by_gene$gene_id)
    write_diagnostic(length(genes_in_both), "gene_overlap.txt", sample,
                    "Number of genes present in both sets:")
    
    # Final filtering
    genes_pass_both <- intersect(
        genes_pass_rpkm,
        uorf_by_gene$gene_id[uorf_by_gene$counts > 0]
    )
    
    # Save detailed info for genes passing RPKM but not uORF filter
    rpkm_only <- setdiff(genes_pass_rpkm, uorf_by_gene$gene_id[uorf_by_gene$counts > 0])
    if(length(rpkm_only) > 0) {
        write_diagnostic(head(rpkm_only), "rpkm_only_genes.txt", sample,
                        "Genes passing RPKM but no uORF reads:")
    }
    
    # Save all analysis results
    write.table(data.frame(
        gene_id = morf_data$gene_id,
        morf_counts = morf_data$counts,
        morf_length = morf_data$length,
        morf_rpkm = morf_data$rpkm,
        uorf_counts = uorf_by_gene$counts[match(morf_data$gene_id, uorf_by_gene$gene_id)]
    ), file=paste0("results/counts/diagnostics/", sample, "_complete_analysis.txt"),
    row.names=FALSE, quote=FALSE, sep="\t")
    
    # Update results summary
    write_diagnostic(c(
        paste("Total genes analyzed:", nrow(morf_data)),
        paste("Genes with RPKM >= 1:", length(genes_pass_rpkm)),
        paste("Genes with uORF reads:", sum(uorf_by_gene$counts > 0)),
        paste("Genes passing both filters:", length(genes_pass_both))
    ), "step_by_step.txt", sample, "Step-by-step results:")
}

print("Analysis complete. Check results/counts/diagnostics/ for detailed information.")
EOF

# Run R script
Rscript filter_by_rpkm_diagnostic.R
