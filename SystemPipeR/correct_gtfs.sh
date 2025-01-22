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
cat << 'EOF' > create_correct_gtf.R
# Load required libraries
library(rtracklayer)
library(GenomicRanges)

print("Starting GTF file creation with correct genomic coordinates...")

# Load uORF predictions
print("Loading uORF predictions...")
uorf_data <- read.table("results/reference_analysis/final_matched_uorfs.txt",
                       header=TRUE, stringsAsFactors=FALSE)

# Load 5' UTR coordinates
print("Loading UTR coordinates...")
utr_bed <- read.table("/global/scratch/users/enricocalvane/riboseq/imb2/ribotish/reference/tair10_5utr.sorted.bed",
                      col.names=c("chr", "start", "end", "transcript_id", "score", "strand"))
# Clean transcript IDs
utr_bed$transcript_id <- sub("transcript:", "", utr_bed$transcript_id)
utr_bed$gene_id <- sub("\\..*", "", utr_bed$transcript_id)

# Load GTF for additional gene information
print("Loading reference GTF...")
gtf <- read.table("/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf",
                 header=FALSE, sep="\t", quote="")

# Function to convert UTR-relative to genomic coordinates
convert_to_genomic <- function(uorf_row, utr_info) {
    gene_id <- sub("_.*", "", uorf_row$uORF_ID)
    gene_utrs <- utr_info[utr_info$gene_id == gene_id,]
    
    if(nrow(gene_utrs) == 0) {
        return(NULL)
    }
    
    # Use the first UTR entry for the gene (should be the primary transcript)
    utr <- gene_utrs[1,]
    
    if(utr$strand == "+") {
        # For positive strand, add UTR start to uORF coordinates
        start <- utr$start + uorf_row$start
        end <- utr$start + uorf_row$end
    } else {
        # For negative strand, subtract from UTR end
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

# Create uORF GTF entries with correct coordinates
print("Converting coordinates and creating GTF entries...")
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
    
    if(i %% 1000 == 0) {
        print(paste("Processed", i, "uORFs"))
    }
}

# Sort the GTF entries
uorf_gtf <- uorf_gtf[order(uorf_gtf$seqname, uorf_gtf$start),]

print(paste("Skipped", skipped_uorfs, "uORFs due to missing UTR information"))

# Create mORF GTF
print("Processing mORF GTF entries...")
morf_gtf <- gtf[gtf$V3 == "CDS",]
colnames(morf_gtf) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute")

# Write GTF files
print("Writing GTF files...")
write.table(uorf_gtf, "results/counts/predicted_uorfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(morf_gtf, "results/counts/morfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

print("GTF file creation complete.")
print(paste("Created uORF GTF with", nrow(uorf_gtf), "entries"))
print(paste("Created mORF GTF with", nrow(morf_gtf), "entries"))

# Verify coordinate ranges
print("\nCoordinate range verification:")
print("uORF coordinate ranges:")
print(summary(uorf_gtf$start))
print(summary(uorf_gtf$end))
print("\nmORF coordinate ranges:")
print(summary(morf_gtf$start))
print(summary(morf_gtf$end))
EOF

# Run R script
Rscript create_correct_gtf.R
