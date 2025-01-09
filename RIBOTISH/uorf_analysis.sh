#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=uorf_analysis
#SBATCH --output=uorf_analysis_%j.out
#SBATCH --error=uorf_analysis_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create output directories
mkdir -p uorf_results/{WT,imb2}/{rep1,rep2}

# Define paths
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"
UTR_BED="reference/tair10_5utr.sorted.bed"

# Step 1: Predict uORFs for each sample
echo "Predicting uORFs in WT replicate 1..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam \
    -g ${GENOME} \
    -t ${GTF} \
    --utr5 ${UTR_BED} \
    -o uorf_results/WT/rep1/predicted_uorfs.txt

echo "Predicting uORFs in WT replicate 2..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam \
    -g ${GENOME} \
    -t ${GTF} \
    --utr5 ${UTR_BED} \
    -o uorf_results/WT/rep2/predicted_uorfs.txt

echo "Predicting uORFs in imb2 replicate 1..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam \
    -g ${GENOME} \
    -t ${GTF} \
    --utr5 ${UTR_BED} \
    -o uorf_results/imb2/rep1/predicted_uorfs.txt

echo "Predicting uORFs in imb2 replicate 2..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam \
    -g ${GENOME} \
    -t ${GTF} \
    --utr5 ${UTR_BED} \
    -o uorf_results/imb2/rep2/predicted_uorfs.txt

# Step 2: Differential analysis
echo "Performing differential analysis..."
ribotish tisdiff \
    --ctrl /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam \
    --treat /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam \
    -g ${GENOME} \
    -t ${GTF} \
    --utr5 ${UTR_BED} \
    -o uorf_results/differential_analysis.txt

# Create summary report
echo "Creating summary report..."
{
    echo "uORF Analysis Summary"
    echo "===================="
    echo
    echo "WT Replicate 1 uORFs:"
    wc -l uorf_results/WT/rep1/predicted_uorfs.txt
    echo
    echo "WT Replicate 2 uORFs:"
    wc -l uorf_results/WT/rep2/predicted_uorfs.txt
    echo
    echo "imb2 Replicate 1 uORFs:"
    wc -l uorf_results/imb2/rep1/predicted_uorfs.txt
    echo
    echo "imb2 Replicate 2 uORFs:"
    wc -l uorf_results/imb2/rep2/predicted_uorfs.txt
    echo
    echo "Differential Analysis Results:"
    echo "----------------------------"
    echo "Total differentially translated uORFs:"
    wc -l uorf_results/differential_analysis.txt
    echo
    echo "Significantly changed uORFs (p < 0.05):"
    awk '$5 < 0.05' uorf_results/differential_analysis.txt | wc -l
} > uorf_results/analysis_summary.txt

echo "Analysis complete. Check uorf_results/analysis_summary.txt for results."
