#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=create_gtf
#SBATCH --output=create_gtf_%j.out
#SBATCH --error=create_gtf_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Load required module
module load r

# Create R script
cat << 'EOF' > create_filtered_gtf.R
# Load required libraries
library(rtracklayer)
library(GenomicRanges)

print("Starting GTF file creation with filtered mORFs...")

# Load uORF predictions
print("Loading uORF predictions...")
uorf_data <- read.table("results/reference_analysis/final_matched_uorfs.txt",
                       header=TRUE, stringsAsFactors=FALSE)

# Get list of genes with uORFs
genes_with_uorfs <- unique(sub("_.*", "", uorf_data$uORF_ID))
print(paste("Found", length(genes_with_uorfs), "genes with uORFs"))

# Load 5' UTR coordinates
print("Loading UTR coordinates...")
utr_bed <- read.table("/global/scratch/users/enricocalvane/riboseq/imb2/ribotish/reference/tair10_5utr.sorted.bed",
                      col.names=c("chr", "start", "end", "transcript_id", "score", "strand"))
utr_bed$transcript_id <- sub("transcript:", "", utr_bed$transcript_id)
utr_bed$gene_id <- sub("\\..*", "", utr_bed$transcript_id)

# Load GTF
print("Loading reference GTF...")
gtf <- read.table("/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf",
                 header=FALSE, sep="\t", quote="")

# Function to convert UTR-relative to genomic coordinates
convert_to_genomic <- function(uorf_row, utr_info) {
    gene_id <- sub("_.*", "", uorf_row$uORF_ID)
    gene_utrs <- utr_info[utr_info$gene_id == gene_id,]
    
    if(nrow(gene_utrs) == 0) return(NULL)
    
    # Use the first UTR entry for the gene
    utr <- gene_utrs[1,]
    
    if(utr$strand == "+") {
        start <- utr$start + uorf_row$start
        end <- utr$start + uorf_row$end
    } else {
        start <- utr$end - uorf_row$end
        end <- utr$end - uorf_row$start
    }
    
    list(
        chr = utr$chr,
        start = start,
        end = end,
        strand = utr$strand
    )
}

# Create uORF GTF entries
print("Creating uORF GTF entries...")
uorf_gtf <- data.frame()
skipped_uorfs <- 0

for(i in seq_len(nrow(uorf_data))) {
    coords <- convert_to_genomic(uorf_data[i,], utr_bed)
    
    if(!is.null(coords)) {
        gene_id <- sub("_.*", "", uorf_data$uORF_ID[i])
        uorf_gtf <- rbind(uorf_gtf, data.frame(
            seqname = coords$chr,
            source = "systemPipeR",
            feature = "uORF",
            start = coords$start,
            end = coords$end,
            score = ".",
            strand = coords$strand,
            frame = ".",
            attribute = sprintf('gene_id "%s"; transcript_id "%s";', 
                              gene_id, uorf_data$uORF_ID[i])
        ))
    } else {
        skipped_uorfs <- skipped_uorfs + 1
    }
}

# Filter and process mORF entries
print("Processing mORF entries...")
# Extract gene IDs from GTF
gtf$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)
gtf$gene_id <- sub("gene:", "", gtf$gene_id)
gtf$transcript_id <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf$V9)
gtf$transcript_id <- sub("transcript:", "", gtf$transcript_id)

# Filter for CDS entries of genes with uORFs
morf_gtf <- gtf[gtf$V3 == "CDS" & gtf$gene_id %in% genes_with_uorfs,]

# Keep only primary transcript for each gene
primary_transcripts <- sapply(genes_with_uorfs, function(gene) {
    gene_rows <- which(morf_gtf$gene_id == gene)
    if(length(gene_rows) == 0) return(NULL)
    
    # Try to find .1 transcript first
    primary_pattern <- paste0(gene, "\\.1$")
    primary_rows <- gene_rows[grep(primary_pattern, morf_gtf$transcript_id[gene_rows])]
    
    # If no .1 transcript, take the first available transcript
    if(length(primary_rows) == 0) {
        primary_rows <- gene_rows[1]
    }
    
    return(primary_rows)
})

# Flatten the list of row indices
primary_rows <- unlist(primary_transcripts)
morf_gtf <- morf_gtf[primary_rows,]

# Format mORF GTF
colnames(morf_gtf) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute", "gene_id", "transcript_id")
morf_gtf <- morf_gtf[,1:9]  # Keep only GTF columns

print(paste("Original number of mORF entries:", nrow(gtf[gtf$V3 == "CDS",])))
print(paste("Filtered mORF entries:", nrow(morf_gtf)))

# Write GTF files
print("Writing GTF files...")
write.table(uorf_gtf, "results/counts/predicted_uorfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(morf_gtf, "results/counts/morfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

print("GTF file creation complete.")
print(paste("Created uORF GTF with", nrow(uorf_gtf), "entries"))
print(paste("Created mORF GTF with", nrow(morf_gtf), "entries"))

# Print example entries for verification
print("\nFirst few entries of uORF GTF:")
print(head(uorf_gtf))
print("\nFirst few entries of mORF GTF:")
print(head(morf_gtf))
EOF

# Run R script
Rscript create_filtered_gtf.R
