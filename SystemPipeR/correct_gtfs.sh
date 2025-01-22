#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=0:10:00
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
cat << 'EOF' > create_morf_gtf.R
# Load required libraries
library(rtracklayer)
library(GenomicRanges)

# Function to safely read a file with error checking
safe_read_table <- function(filepath, ...) {
    if (!file.exists(filepath)) {
        stop(paste("File not found:", filepath))
    }
    tryCatch({
        read.table(filepath, ...)
    }, error = function(e) {
        stop(paste("Error reading file:", filepath, "\nError:", e$message))
    })
}

print("Starting GTF file creation with proper mORF definitions...")

# Step 1: Load and validate reference data
print("Loading input files...")

# Load uORF predictions
uorf_data <- safe_read_table("results/reference_analysis/final_matched_uorfs.txt", 
                            header=TRUE, stringsAsFactors=FALSE)
print(paste("Loaded", nrow(uorf_data), "uORF predictions"))

# Get unique genes with uORFs
genes_with_uorfs <- unique(sub("_.*", "", uorf_data$uORF_ID))
print(paste("Found", length(genes_with_uorfs), "genes with uORFs"))

# Load 5' UTR coordinates
utr_bed <- safe_read_table(
    "/global/scratch/users/enricocalvane/riboseq/imb2/ribotish/reference/tair10_5utr.sorted.bed",
    col.names=c("chr", "start", "end", "transcript_id", "score", "strand"),
    stringsAsFactors=FALSE
)
print(paste("Loaded", nrow(utr_bed), "UTR entries"))

# Clean UTR transcript IDs
utr_bed$transcript_id <- sub("transcript:", "", utr_bed$transcript_id)
utr_bed$gene_id <- sub("\\..*", "", utr_bed$transcript_id)

# Load GTF file
print("Loading reference GTF...")
gtf_file <- "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"
gtf <- safe_read_table(gtf_file, sep="\t", quote="", comment.char="#", stringsAsFactors=FALSE)
names(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
print(paste("Loaded GTF file with", nrow(gtf), "entries"))

# Step 2: Process GTF entries
print("Processing GTF entries...")

# Extract gene and transcript IDs from GTF attributes
gtf$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$attribute)
gtf$gene_id <- sub("gene:", "", gtf$gene_id)
gtf$transcript_id <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf$attribute)
gtf$transcript_id <- sub("transcript:", "", gtf$transcript_id)

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
    
    if(i %% 1000 == 0) {
        print(paste("Processed", i, "uORFs"))
    }
}

print(paste("Created", nrow(uorf_gtf), "uORF entries"))

# Step 3: Create complete mORFs
# Function to create complete mORF entry
create_complete_morf <- function(cds_entries) {
    if(nrow(cds_entries) == 0) return(NULL)
    
    # Sort CDS entries by position
    cds_entries <- cds_entries[order(cds_entries$start),]
    
    # Create single mORF entry
    data.frame(
        seqname = cds_entries$seqname[1],
        source = "systemPipeR",
        feature = "mORF",
        start = min(cds_entries$start),
        end = max(cds_entries$end),
        score = ".",
        strand = cds_entries$strand[1],
        frame = ".",
        attribute = cds_entries$attribute[1]
    )
}

# Process mORFs
print("Processing mORF entries...")
cds_entries <- gtf[gtf$feature == "CDS" & gtf$gene_id %in% genes_with_uorfs,]
print(paste("Found", nrow(cds_entries), "CDS entries for", length(genes_with_uorfs), "genes"))

# Create complete mORFs
morf_gtf <- data.frame()
processed_genes <- 0

for(gene in genes_with_uorfs) {
    # Get all CDS entries for this gene
    gene_cds <- cds_entries[cds_entries$gene_id == gene,]
    if(nrow(gene_cds) == 0) next
    
    # Find primary transcript (prefer .1)
    transcripts <- unique(gene_cds$transcript_id)
    primary_transcript <- transcripts[grep("\\.1$", transcripts)]
    if(length(primary_transcript) == 0) {
        primary_transcript <- transcripts[1]
    }
    
    # Get CDS entries for primary transcript
    transcript_cds <- gene_cds[gene_cds$transcript_id == primary_transcript,]
    
    # Create complete mORF entry
    morf_entry <- create_complete_morf(transcript_cds)
    if(!is.null(morf_entry)) {
        morf_gtf <- rbind(morf_gtf, morf_entry)
    }
    
    processed_genes <- processed_genes + 1
    if(processed_genes %% 1000 == 0) {
        print(paste("Processed", processed_genes, "genes"))
    }
}

print("\nSummary Statistics:")
print(paste("Total genes with uORFs:", length(genes_with_uorfs)))
print(paste("Genes with mORF entries:", nrow(morf_gtf)))
print(paste("Original CDS entries:", nrow(cds_entries)))
print(paste("Final mORF entries:", nrow(morf_gtf)))

# Write GTF files
print("\nWriting GTF files...")
write.table(uorf_gtf, "results/counts/predicted_uorfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(morf_gtf, "results/counts/morfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

print("GTF file creation complete.")
print(paste("Created uORF GTF with", nrow(uorf_gtf), "entries"))
print(paste("Created mORF GTF with", nrow(morf_gtf), "entries"))

# Print example entries and coordinate ranges
print("\nFirst few entries of uORF GTF:")
print(head(uorf_gtf))
print("\nFirst few entries of mORF GTF:")
print(head(morf_gtf))

print("\nCoordinate ranges for verification:")
print("uORF coordinate ranges:")
print(summary(uorf_gtf$start))
print(summary(uorf_gtf$end))
print("\nmORF coordinate ranges:")
print(summary(morf_gtf$start))
print(summary(morf_gtf$end))
EOF

# Run R script
Rscript create_morf_gtf.R
