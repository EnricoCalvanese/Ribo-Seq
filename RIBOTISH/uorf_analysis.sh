#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=uorf_analysis
#SBATCH --output=uorf_analysis_%j.out
#SBATCH --error=uorf_analysis_%j.err

# Load required modules
module load bio/bedtools2/2.31.0-gcc-11.4.0

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create output directories
mkdir -p uorf_results/{AUG,nonAUG}/{WT,imb2}/{rep1,rep2}

# Define paths
GENOME="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf"
UTR_BED="reference/tair10_5utr.sorted.bed"

# Function to convert ribotish output to BED format
convert_to_bed() {
    input=$1
    output=$2
    awk 'NR>1 {
        split($5,a,":"); 
        split(a[2],b,"-");
        split(b[2],c,":");
        chr=a[1];
        start=b[1];
        end=c[1];
        strand=c[2];
        print chr"\t"start"\t"end"\t"$1"_"$2"\t.\t"strand
    }' "$input" > "$output"
}

# Function to filter predictions by 5' UTR regions
filter_by_utr() {
    predictions=$1
    output=$2
    tmp_bed="tmp.bed"
    
    # Convert predictions to BED
    convert_to_bed "$predictions" "$tmp_bed"
    
    # Intersect with 5' UTR regions
    bedtools intersect -a "$tmp_bed" -b "$UTR_BED" -wa > "intersect.bed"
    
    # Keep only predictions that overlap with 5' UTRs
    awk 'NR==1{print}' "$predictions" > "$output"  # Keep header
    while read -r line; do
        id=$(echo "$line" | cut -f4)
        gene=$(echo "$id" | cut -f1 -d"_")
        transcript=$(echo "$id" | cut -f2 -d"_")
        grep "^$gene[[:space:]].*$transcript" "$predictions" >> "$output"
    done < "intersect.bed"
    
    rm "$tmp_bed" "intersect.bed"
}

# Process each sample for AUG uORFs
echo "Processing AUG uORFs..."
for condition in WT imb2; do
    for rep in rep1 rep2; do
        sample=""
        case "${condition}_${rep}" in
            "WT_rep1") sample="LZT103-1" ;;
            "WT_rep2") sample="LZT103-2" ;;
            "imb2_rep1") sample="LZT104-1" ;;
            "imb2_rep2") sample="LZT104-2" ;;
        esac
        
        echo "Processing ${condition} ${rep} (${sample})..."
        
        # Run prediction
        ribotish predict \
            -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample}_uniq_sort.bam \
            -g ${GTF} \
            -f ${GENOME} \
            -o uorf_results/AUG/${condition}/${rep}/raw_predictions.txt \
            --geneformat gtf \
            -p 24 \
            -v
            
        # Filter predictions to 5' UTR regions
        filter_by_utr \
            "uorf_results/AUG/${condition}/${rep}/raw_predictions.txt" \
            "uorf_results/AUG/${condition}/${rep}/filtered_predictions.txt"
    done
done

# Process each sample for non-AUG uORFs
echo "Processing non-AUG uORFs..."
for condition in WT imb2; do
    for rep in rep1 rep2; do
        sample=""
        case "${condition}_${rep}" in
            "WT_rep1") sample="LZT103-1" ;;
            "WT_rep2") sample="LZT103-2" ;;
            "imb2_rep1") sample="LZT104-1" ;;
            "imb2_rep2") sample="LZT104-2" ;;
        esac
        
        echo "Processing ${condition} ${rep} (${sample})..."
        
        # Run prediction
        ribotish predict \
            -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample}_uniq_sort.bam \
            -g ${GTF} \
            -f ${GENOME} \
            -o uorf_results/nonAUG/${condition}/${rep}/raw_predictions.txt \
            --geneformat gtf \
            --alt \
            --altcodons "CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT" \
            -p 24 \
            -v
            
        # Filter predictions to 5' UTR regions
        filter_by_utr \
            "uorf_results/nonAUG/${condition}/${rep}/raw_predictions.txt" \
            "uorf_results/nonAUG/${condition}/${rep}/filtered_predictions.txt"
    done
done

# Run differential analysis with filtered predictions
echo "Performing differential analysis..."

# For AUG uORFs
ribotish tisdiff \
    -1 "uorf_results/AUG/WT/rep1/filtered_predictions.txt,uorf_results/AUG/WT/rep2/filtered_predictions.txt" \
    -2 "uorf_results/AUG/imb2/rep1/filtered_predictions.txt,uorf_results/AUG/imb2/rep2/filtered_predictions.txt" \
    -a "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam" \
    -b "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam" \
    -g ${GENOME} \
    -o uorf_results/AUG/differential_analysis.txt \
    -p 24 \
    -v

# For non-AUG uORFs
ribotish tisdiff \
    -1 "uorf_results/nonAUG/WT/rep1/filtered_predictions.txt,uorf_results/nonAUG/WT/rep2/filtered_predictions.txt" \
    -2 "uorf_results/nonAUG/imb2/rep1/filtered_predictions.txt,uorf_results/nonAUG/imb2/rep2/filtered_predictions.txt" \
    -a "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam" \
    -b "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam" \
    -g ${GENOME} \
    -o uorf_results/nonAUG/differential_analysis.txt \
    -p 24 \
    -v

# Create summary report
echo "Creating summary report..."
{
    echo "uORF Analysis Summary (5' UTR filtered)"
    echo "======================================"
    echo
    echo "AUG-initiated uORFs"
    echo "------------------"
    echo "Number of uORFs detected in each sample:"
    for condition in WT imb2; do
        for rep in rep1 rep2; do
            echo "${condition} ${rep}:"
            wc -l "uorf_results/AUG/${condition}/${rep}/filtered_predictions.txt"
        done
    done
    
    if [ -f "uorf_results/AUG/differential_analysis.txt" ]; then
        echo
        echo "Differential Analysis Results (AUG):"
        echo "Total differentially translated AUG uORFs:"
        wc -l "uorf_results/AUG/differential_analysis.txt"
        echo "Significantly changed AUG uORFs (p < 0.05):"
        awk '$5 < 0.05' "uorf_results/AUG/differential_analysis.txt" | wc -l
        
        echo
        echo "Top 10 most significant changes (AUG):"
        (echo -e "Gene\tPosition\tLog2FC\tP-value" && \
        sort -k5,5n "uorf_results/AUG/differential_analysis.txt" | head -10) | column -t
    fi
    
    echo
    echo "Non-AUG uORFs"
    echo "-------------"
    echo "Number of uORFs detected in each sample:"
    for condition in WT imb2; do
        for rep in rep1 rep2; do
            echo "${condition} ${rep}:"
            wc -l "uorf_results/nonAUG/${condition}/${rep}/filtered_predictions.txt"
        done
    done
    
    if [ -f "uorf_results/nonAUG/differential_analysis.txt" ]; then
        echo
        echo "Differential Analysis Results (non-AUG):"
        echo "Total differentially translated non-AUG uORFs:"
        wc -l "uorf_results/nonAUG/differential_analysis.txt"
        echo "Significantly changed non-AUG uORFs (p < 0.05):"
        awk '$5 < 0.05' "uorf_results/nonAUG/differential_analysis.txt" | wc -l
        
        echo
        echo "Top 10 most significant changes (non-AUG):"
        (echo -e "Gene\tPosition\tLog2FC\tP-value" && \
        sort -k5,5n "uorf_results/nonAUG/differential_analysis.txt" | head -10) | column -t
    fi
} > uorf_results/filtered_analysis_summary.txt

echo "Analysis complete. Check uorf_results/filtered_analysis_summary.txt for results."
