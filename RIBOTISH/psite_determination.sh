#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=ribotish_quality
#SBATCH --output=ribotish_quality_%j.out
#SBATCH --error=ribotish_quality_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Run the Python wrapper script
python3 ribotish_wrapper.py

# Create a summary report if the script completes successfully
if [ $? -eq 0 ]; then
    echo "Creating summary report..."
    {
        echo "Quality Control Summary"
        echo "======================"
        echo ""
        
        for sample in LZT103-1 LZT103-2 LZT104-1 LZT104-2; do
            echo "Sample: ${sample}"
            echo "-------------------------"
            if [ -f "quality_results/${sample}_qual.txt" ]; then
                echo "Quality control completed successfully"
                head -n 5 "quality_results/${sample}_qual.txt"
            else
                echo "Quality control failed or incomplete"
            fi
            echo ""
        done
    } > quality_results/summary.txt
fi
