#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --job-name=rpkm_filter
#SBATCH --output=rpkm_filter_%j.out
#SBATCH --error=rpkm_filter_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script
cat << 'EOF' > fixed_rpkm.R
# Create diagnostic directory
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

# Function to clean gene IDs from GTF attributes
clean_gene_id <- function(attr) {
    # Extract text between gene_id "gene: and first "
    id <- sub('.*gene_id "gene:([^"]+)".*', "\\1", attr)
    return(id)
}

print("Loading and processing GTF file...")
# Read GTF files
morf_gtf <- read.table("results/counts/morfs.gtf", header=FALSE, sep="\t", quote="")

# Calculate mORF lengths with proper gene ID extraction
morf_lengths <- data.frame(
    gene_id = sapply(morf_gtf$V9, clean_gene_id),
    length = morf_gtf$V5 - morf_gtf$V4 + 1
)

# Verify gene ID cleaning worked
write_diagnostic(head(morf_lengths), "length_extraction.txt", 
                note="First few mORF lengths with cleaned gene IDs:")

# Function to calculate RPKM
calculate_rpkm <- function(counts, lengths, total_reads) {
    rpkm <- (counts * 10^9) / (lengths * total_reads)
    return(rpkm)
}

# Process each sample
samples <- c("LZT103-1", "LZT103-2", "LZT104-1", "LZT104-2")
results_summary <- data.frame()

for(sample in samples) {
    print(paste("\nProcessing sample:", sample))
    
    # Read count files
    morf_counts <- read.table(paste0("results/counts/", sample, "_morf_counts.txt"))
    uorf_counts <- read.table(paste0("results/counts/", sample, "_uorf_counts.txt"))
    
    # Clean gene IDs in count data
    morf_counts$V1 <- sub("gene:", "", morf_counts$V1)
    
    # Filter out HTSeq summary statistics
    morf_data <- morf_counts[!grepl("^__", morf_counts$V1),]
    
    # Calculate total mapped reads (excluding HTSeq stats)
    total_reads <- sum(morf_data$V2) + 
                  sum(uorf_counts$V2[!grepl("^__", uorf_counts$V1)])
    
    write_diagnostic(total_reads, "total_reads.txt", sample, 
                    "Total mapped reads (excluding HTSeq stats):")
    
    # Match lengths to count data
    morf_data <- data.frame(
        gene_id = morf_data$V1,
        counts = morf_data$V2
    )
    
    # Add lengths with diagnostic output
    morf_data$length <- morf_lengths$length[match(morf_data$gene_id, morf_lengths$gene_id)]
    
    # Check for missing length matches
    missing_lengths <- is.na(morf_data$length)
    if(sum(missing_lengths) > 0) {
        write_diagnostic(
            morf_data[missing_lengths, ],
            "missing_length_details.txt",
            sample,
            paste("Genes with missing lengths:", sum(missing_lengths))
        )
    }
    
    # Calculate RPKM
    morf_data$rpkm <- calculate_rpkm(morf_data$counts, morf_data$length, total_reads)
    
    # RPKM distribution check
    write_diagnostic(summary(morf_data$rpkm), "rpkm_distribution.txt", sample,
                    "RPKM distribution:")
    write_diagnostic(head(morf_data[order(-morf_data$rpkm),]), "top_rpkm.txt", sample,
                    "Top RPKM values:")
    
    # Filter for RPKM >= 1
    genes_pass_rpkm <- morf_data$gene_id[morf_data$rpkm >= 1]
    write_diagnostic(length(genes_pass_rpkm), "rpkm_filter.txt", sample,
                    "Number of genes passing RPKM filter:")
    
    # Process uORF data
    uorf_data <- data.frame(
        gene_id = sub("_.*", "", uorf_counts$V1[!grepl("^__", uorf_counts$V1)]),
        counts = uorf_counts$V2[!grepl("^__", uorf_counts$V1)]
    )
    
    # Sum counts by gene
    uorf_by_gene <- aggregate(counts ~ gene_id, data=uorf_data, sum)
    
    # Find genes passing both filters
    genes_pass_both <- intersect(
        genes_pass_rpkm,
        uorf_by_gene$gene_id[uorf_by_gene$counts > 0]
    )
    
    # Save filtered gene list
    write.table(genes_pass_both,
                paste0("results/counts/", sample, "_filtered_genes.txt"),
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    # Save complete analysis
    analysis_data <- merge(
        morf_data,
        uorf_by_gene,
        by="gene_id",
        all.x=TRUE
    )
    names(analysis_data)[5] <- "uorf_counts"
    write.table(analysis_data,
                paste0("results/counts/diagnostics/", sample, "_complete_analysis.txt"),
                quote=FALSE, row.names=FALSE, sep="\t")
    
    # Update results summary
    results_summary <- rbind(results_summary, data.frame(
        sample = sample,
        total_genes = nrow(morf_data),
        genes_rpkm = length(genes_pass_rpkm),
        genes_with_uorf_reads = sum(uorf_by_gene$counts > 0),
        genes_pass_both = length(genes_pass_both)
    ))
}

# Create final summary
summary_report <- "RPKM and Read Coverage Filtering Summary\n"
summary_report <- paste0(summary_report, "=====================================\n\n")

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

writeLines(summary_report, "results/counts/filtering_summary.txt")
print(summary_report)
EOF

# Run R script
Rscript fixed_rpkm.R
