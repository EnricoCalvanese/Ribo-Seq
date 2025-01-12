#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=count_reads
#SBATCH --output=count_reads_%j.out
#SBATCH --error=count_reads_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Create directory for counts
mkdir -p results/counts

# Create R script to generate GTF files
cat << 'EOF' > create_gtf.R
# Load required libraries
library(rtracklayer)

# Function to create GTF entry
create_gtf_entry <- function(chr, start, end, strand, feature_type, gene_id, transcript_id) {
    data.frame(
        seqname = chr,
        source = "systemPipeR",
        feature = feature_type,
        start = as.integer(start),
        end = as.integer(end),
        score = ".",
        strand = strand,
        frame = ".",
        attribute = sprintf('gene_id "%s"; transcript_id "%s";', gene_id, transcript_id)
    )
}

# Load uORF predictions
uorf_data <- read.table("results/reference_analysis/final_matched_uorfs.txt",
                       header=TRUE, stringsAsFactors=FALSE)

# Load original GTF to get chromosome and strand information
gtf <- read.table("/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf",
                 header=FALSE, sep="\t", quote="")

# Extract chromosome and strand information for each gene
gene_info <- unique(gtf[gtf$V3 == "transcript",
                       c("V1", "V7", "V9")])  # chr, strand, attributes
gene_info$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gene_info$V9)
gene_info$gene_id <- sub("gene:", "", gene_info$gene_id)

# Create uORF GTF entries
uorf_gtf <- data.frame()
for(i in seq_len(nrow(uorf_data))) {
    gene_id <- sub("_.*", "", uorf_data$uORF_ID[i])
    gene_row <- gene_info[gene_info$gene_id == gene_id,]
    
    if(nrow(gene_row) > 0) {
        uorf_gtf <- rbind(uorf_gtf,
            create_gtf_entry(
                chr = gene_row$V1[1],
                start = uorf_data$start[i],
                end = uorf_data$end[i],
                strand = gene_row$V7[1],
                feature_type = "uORF",
                gene_id = gene_id,
                transcript_id = uorf_data$uORF_ID[i]
            )
        )
    }
}

# Create mORF GTF entries
morf_gtf <- gtf[gtf$V3 == "CDS",]
colnames(morf_gtf) <- c("seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute")

# Write GTF files
write.table(uorf_gtf, "results/counts/predicted_uorfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(morf_gtf, "results/counts/morfs.gtf",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
EOF

# Run R script to create GTF files
Rscript create_gtf.R

# Function to run htseq-count
run_htseq_count() {
    local bam=$1
    local gtf=$2
    local output=$3
    local feature_type=$4
    
    htseq-count \
        -f bam \
        -r name \
        -s yes \
        -t ${feature_type} \
        -i gene_id \
        ${bam} ${gtf} > ${output}
}

# Process each BAM file for both uORFs and mORFs
BAM_DIR="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
for sample in LZT103-1 LZT103-2 LZT104-1 LZT104-2; do
    echo "Processing ${sample}..."
    
    # Count reads for uORFs
    echo "Counting uORF reads..."
    run_htseq_count \
        ${BAM_DIR}/${sample}_uniq_sort.bam \
        results/counts/predicted_uorfs.gtf \
        results/counts/${sample}_uorf_counts.txt \
        uORF
    
    # Count reads for mORFs
    echo "Counting mORF reads..."
    run_htseq_count \
        ${BAM_DIR}/${sample}_uniq_sort.bam \
        results/counts/morfs.gtf \
        results/counts/${sample}_morf_counts.txt \
        CDS
done

echo "Counting complete. Results are in results/counts/"
