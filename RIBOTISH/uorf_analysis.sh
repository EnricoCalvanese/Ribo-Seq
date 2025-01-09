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

# Create output directories for both analyses
mkdir -p uorf_results/{AUG,nonAUG}/{WT,imb2}/{rep1,rep2}

# Define common paths
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"

echo "Starting AUG-initiated uORF analysis..."
echo "======================================="

# Step 1A: Predict AUG uORFs for each sample
echo "Predicting AUG uORFs in WT replicate 1..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/AUG/WT/rep1/predicted_uorfs.txt \
    --geneformat gtf \
    -p 24 \
    -v

echo "Predicting AUG uORFs in WT replicate 2..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/AUG/WT/rep2/predicted_uorfs.txt \
    --geneformat gtf \
    -p 24 \
    -v

echo "Predicting AUG uORFs in imb2 replicate 1..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/AUG/imb2/rep1/predicted_uorfs.txt \
    --geneformat gtf \
    -p 24 \
    -v

echo "Predicting AUG uORFs in imb2 replicate 2..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/AUG/imb2/rep2/predicted_uorfs.txt \
    --geneformat gtf \
    -p 24 \
    -v

echo "Starting non-AUG uORF analysis..."
echo "================================"

# Step 1B: Predict non-AUG uORFs for each sample
echo "Predicting non-AUG uORFs in WT replicate 1..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/nonAUG/WT/rep1/predicted_uorfs.txt \
    --geneformat gtf \
    --alt \
    --altcodons "CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT" \
    -p 24 \
    -v

echo "Predicting non-AUG uORFs in WT replicate 2..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/nonAUG/WT/rep2/predicted_uorfs.txt \
    --geneformat gtf \
    --alt \
    --altcodons "CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT" \
    -p 24 \
    -v

echo "Predicting non-AUG uORFs in imb2 replicate 1..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/nonAUG/imb2/rep1/predicted_uorfs.txt \
    --geneformat gtf \
    --alt \
    --altcodons "CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT" \
    -p 24 \
    -v

echo "Predicting non-AUG uORFs in imb2 replicate 2..."
ribotish predict \
    -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam \
    -g ${GTF} \
    -f ${GENOME} \
    -o uorf_results/nonAUG/imb2/rep2/predicted_uorfs.txt \
    --geneformat gtf \
    --alt \
    --altcodons "CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT" \
    -p 24 \
    -v

# Step 2A: Differential analysis for AUG uORFs
echo "Performing differential analysis for AUG uORFs..."
ribotish tisdiff \
    -1 "uorf_results/AUG/WT/rep1/predicted_uorfs.txt,uorf_results/AUG/WT/rep2/predicted_uorfs.txt" \
    -2 "uorf_results/AUG/imb2/rep1/predicted_uorfs.txt,uorf_results/AUG/imb2/rep2/predicted_uorfs.txt" \
    -a "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam" \
    -b "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam" \
    -g ${GENOME} \
    -o uorf_results/AUG/differential_analysis.txt \
    -p 24 \
    -v

# Step 2B: Differential analysis for non-AUG uORFs
echo "Performing differential analysis for non-AUG uORFs..."
ribotish tisdiff \
    -1 "uorf_results/nonAUG/WT/rep1/predicted_uorfs.txt,uorf_results/nonAUG/WT/rep2/predicted_uorfs.txt" \
    -2 "uorf_results/nonAUG/imb2/rep1/predicted_uorfs.txt,uorf_results/nonAUG/imb2/rep2/predicted_uorfs.txt" \
    -a "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam" \
    -b "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam" \
    -g ${GENOME} \
    -o uorf_results/nonAUG/differential_analysis.txt \
    -p 24 \
    -v

# Create summary report
echo "Creating summary report..."
{
    echo "uORF Analysis Summary"
    echo "===================="
    echo
    echo "AUG-initiated uORFs"
    echo "------------------"
    echo "Number of uORFs detected in each sample:"
    echo "WT Replicate 1:"
    wc -l uorf_results/AUG/WT/rep1/predicted_uorfs.txt
    echo "WT Replicate 2:"
    wc -l uorf_results/AUG/WT/rep2/predicted_uorfs.txt
    echo "imb2 Replicate 1:"
    wc -l uorf_results/AUG/imb2/rep1/predicted_uorfs.txt
    echo "imb2 Replicate 2:"
    wc -l uorf_results/AUG/imb2/rep2/predicted_uorfs.txt
    echo
    echo "Differential Analysis Results (AUG):"
    echo "Total differentially translated AUG uORFs:"
    wc -l uorf_results/AUG/differential_analysis.txt
    echo "Significantly changed AUG uORFs (p < 0.05):"
    awk '$5 < 0.05' uorf_results/AUG/differential_analysis.txt | wc -l
    echo
    echo "Non-AUG uORFs"
    echo "-------------"
    echo "Number of uORFs detected in each sample:"
    echo "WT Replicate 1:"
    wc -l uorf_results/nonAUG/WT/rep1/predicted_uorfs.txt
    echo "WT Replicate 2:"
    wc -l uorf_results/nonAUG/WT/rep2/predicted_uorfs.txt
    echo "imb2 Replicate 1:"
    wc -l uorf_results/nonAUG/imb2/rep1/predicted_uorfs.txt
    echo "imb2 Replicate 2:"
    wc -l uorf_results/nonAUG/imb2/rep2/predicted_uorfs.txt
    echo
    echo "Differential Analysis Results (non-AUG):"
    echo "Total differentially translated non-AUG uORFs:"
    wc -l uorf_results/nonAUG/differential_analysis.txt
    echo "Significantly changed non-AUG uORFs (p < 0.05):"
    awk '$5 < 0.05' uorf_results/nonAUG/differential_analysis.txt | wc -l
    
    echo
    echo "Top 10 most significant changes (AUG):"
    (echo -e "Gene\tPosition\tLog2FC\tP-value" && \
    sort -k5,5n uorf_results/AUG/differential_analysis.txt | head -10) | column -t
    
    echo
    echo "Top 10 most significant changes (non-AUG):"
    (echo -e "Gene\tPosition\tLog2FC\tP-value" && \
    sort -k5,5n uorf_results/nonAUG/differential_analysis.txt | head -10) | column -t
} > uorf_results/analysis_summary.txt

echo "Analysis complete. Check uorf_results/analysis_summary.txt for results."
