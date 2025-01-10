#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=2:00:00
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

# Function to filter predictions based on start positions in 5' UTRs
filter_predictions() {
    input=$1
    output=$2
    tmp_dir=$(mktemp -d)
    
    # Convert to BED format, keeping only start positions and ensuring proper field count
    awk -v OFS="\t" 'NR>1 {
        split($5,a,":");
        split(a[2],b,"-");
        split(b[2],c,":");
        chr=a[1];
        if(c[2] == "+") {
            print chr, b[1], b[1]+1, $1"_"$2, ".", c[2]
        } else {
            print chr, c[1]-1, c[1], $1"_"$2, ".", c[2]
        }
    }' "$input" > "${tmp_dir}/starts.bed"
    
    # Intersect with 5' UTRs
    bedtools intersect -a "${tmp_dir}/starts.bed" -b "$UTR_BED" -wa > "${tmp_dir}/overlapping.bed"
    
    # Create filtered output with original data
    head -n 1 "$input" > "$output"  # Keep header
    while read -r line; do
        id=$(echo "$line" | cut -f4)
        gene=$(echo "$id" | cut -f1 -d"_")
        transcript=$(echo "$id" | cut -f2 -d"_")
        grep -P "^${gene}\t${transcript}" "$input" >> "$output" || true
    done < "${tmp_dir}/overlapping.bed"
    
    # Cleanup
    rm -rf "${tmp_dir}"
}

# Process samples
process_sample() {
    local condition=$1
    local rep=$2
    local type=$3
    local sample

    case "${condition}_${rep}" in
        "WT_rep1") sample="LZT103-1" ;;
        "WT_rep2") sample="LZT103-2" ;;
        "imb2_rep1") sample="LZT104-1" ;;
        "imb2_rep2") sample="LZT104-2" ;;
    esac

    echo "Processing ${type} uORFs for ${condition} ${rep} (${sample})..."
    
    # Create output directory
    mkdir -p "uorf_results/${type}/${condition}/${rep}"
    
    # Base prediction command
    cmd="ribotish predict \
        -b /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/${sample}_uniq_sort.bam \
        -g ${GTF} \
        -f ${GENOME} \
        -o uorf_results/${type}/${condition}/${rep}/raw_predictions.txt \
        --geneformat gtf \
        -p 24 \
        -v"
    
    # Add non-AUG parameters if needed
    if [ "$type" == "nonAUG" ]; then
        cmd="$cmd --alt --altcodons CTG,GTG,TTG,ACG,AGG,AAG,ATC,ATA,ATT"
    fi
    
    # Run prediction
    eval $cmd
    
    # Filter predictions if raw file exists
    if [ -f "uorf_results/${type}/${condition}/${rep}/raw_predictions.txt" ]; then
        filter_predictions \
            "uorf_results/${type}/${condition}/${rep}/raw_predictions.txt" \
            "uorf_results/${type}/${condition}/${rep}/predictions.txt"
    else
        echo "Warning: No raw predictions found for ${condition} ${rep}"
        return 1
    fi
}

# Process all samples
echo "Processing AUG uORFs..."
for condition in WT imb2; do
    for rep in rep1 rep2; do
        process_sample "$condition" "$rep" "AUG"
    done
done

echo "Processing non-AUG uORFs..."
for condition in WT imb2; do
    for rep in rep1 rep2; do
        process_sample "$condition" "$rep" "nonAUG"
    done
done

# Function to check if files exist and have content
check_files() {
    local type=$1
    local error=0
    
    for condition in WT imb2; do
        for rep in rep1 rep2; do
            local file="uorf_results/${type}/${condition}/${rep}/predictions.txt"
            if [ ! -f "$file" ] || [ ! -s "$file" ]; then
                echo "Error: Missing or empty file: $file"
                error=1
            fi
        done
    done
    
    return $error
}

# Perform differential analysis if files exist
run_diff_analysis() {
    local type=$1
    
    if check_files "$type"; then
        echo "Running differential analysis for ${type} uORFs..."
        ribotish tisdiff \
            -1 "uorf_results/${type}/WT/rep1/predictions.txt,uorf_results/${type}/WT/rep2/predictions.txt" \
            -2 "uorf_results/${type}/imb2/rep1/predictions.txt,uorf_results/${type}/imb2/rep2/predictions.txt" \
            -a "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT103-2_uniq_sort.bam" \
            -b "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-1_uniq_sort.bam,/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/LZT104-2_uniq_sort.bam" \
            -g ${GENOME} \
            -o "uorf_results/${type}/differential_analysis.txt" \
            -p 24 \
            -v
    else
        echo "Skipping differential analysis for ${type} due to missing files"
    fi
}

# Run differential analysis
run_diff_analysis "AUG"
run_diff_analysis "nonAUG"

# Create summary report
echo "Creating summary report..."
{
    echo "uORF Analysis Summary (Start Position Filtered)"
    echo "=============================================="
    echo
    
    for type in AUG nonAUG; do
        echo "${type}-initiated uORFs"
        echo "----------------------"
        echo "Number of uORFs detected in each sample:"
        for condition in WT imb2; do
            for rep in rep1 rep2; do
                if [ -f "uorf_results/${type}/${condition}/${rep}/predictions.txt" ]; then
                    echo "${condition} ${rep}:"
                    wc -l "uorf_results/${type}/${condition}/${rep}/predictions.txt"
                else
                    echo "${condition} ${rep}: File not found"
                fi
            done
        done
        
        echo
        if [ -f "uorf_results/${type}/differential_analysis.txt" ]; then
            echo "Differential Analysis Results:"
            echo "Total differentially translated uORFs:"
            wc -l "uorf_results/${type}/differential_analysis.txt"
            echo "Significantly changed uORFs (p < 0.05):"
            awk '$5 < 0.05' "uorf_results/${type}/differential_analysis.txt" | wc -l
            
            echo
            echo "Top 10 most significant changes:"
            {
                echo -e "Gene\tPosition\tLog2FC\tP-value"
                sort -k5,5n "uorf_results/${type}/differential_analysis.txt" | head -10
            } | column -t
        else
            echo "No differential analysis results found"
        fi
        echo
    done
} > uorf_results/analysis_summary.txt

echo "Analysis complete. Check uorf_results/analysis_summary.txt for results."
