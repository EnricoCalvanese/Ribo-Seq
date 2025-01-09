#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=uorf_analysis
#SBATCH --output=uorf_analysis_%j.out
#SBATCH --error=uorf_analysis_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create output directories
mkdir -p results/AUG_uORFs
mkdir -p results/nonAUG_uORFs

# Define input files
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"
WT_REPS="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam"
IMB2_REPS="/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam"

# Analysis for AUG-initiated uORFs
echo "Starting analysis of AUG-initiated uORFs"
ribo-TISH -t aTIS \
    --ctrl ${WT_REPS} \
    --treat ${IMB2_REPS} \
    -g ${GENOME} \
    -s ${GTF} \
    --anno-start AUG \
    --output results/AUG_uORFs

# Analysis for non-canonical start codons
echo "Starting analysis of non-canonical start codons"
ribo-TISH -t aTIS \
    --ctrl ${WT_REPS} \
    --treat ${IMB2_REPS} \
    -g ${GENOME} \
    -s ${GTF} \
    --anno-start CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT \
    --output results/nonAUG_uORFs

echo "uORF analysis complete"

# Create a summary of the results
echo "Creating results summary"
Rscript - <<'EOF'
# Read results
aug_results <- read.table("results/AUG_uORFs/result_aTIS.txt", header=TRUE, sep="\t")
nonaug_results <- read.table("results/nonAUG_uORFs/result_aTIS.txt", header=TRUE, sep="\t")

# Filter for 5' UTR regions and significant changes
aug_utr <- aug_results[aug_results$region == "5UTR" & aug_results$pvalue < 0.05,]
nonaug_utr <- nonaug_results[nonaug_results$region == "5UTR" & nonaug_results$pvalue < 0.05,]

# Create summary
sink("results/uORF_analysis_summary.txt")
cat("Summary of differential uORF translation analysis\n\n")
cat("AUG-initiated uORFs:\n")
cat("Total significant differences:", nrow(aug_utr), "\n")
cat("Upregulated in imb2:", sum(aug_utr$log2FC > 0), "\n")
cat("Downregulated in imb2:", sum(aug_utr$log2FC < 0), "\n\n")

cat("Non-AUG-initiated uORFs:\n")
cat("Total significant differences:", nrow(nonaug_utr), "\n")
cat("Upregulated in imb2:", sum(nonaug_utr$log2FC > 0), "\n")
cat("Downregulated in imb2:", sum(nonaug_utr$log2FC < 0), "\n")
sink()
EOF
