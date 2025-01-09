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

# Create a detailed summary report if the script completes successfully
if [ $? -eq 0 ]; then
    echo "Creating summary report..."
    {
        echo "Quality Control Summary Report"
        echo "============================"
        echo "Generated: $(date)"
        echo ""
        
        for sample in LZT103-1 LZT103-2 LZT104-1 LZT104-2; do
            echo "Sample: ${sample}"
            echo "-------------------------"
            qual_file="quality_results/${sample}_qual.txt"
            if [ -f "$qual_file" ]; then
                echo "✓ Quality control completed successfully"
                echo "Key metrics:"
                echo "-------------------"
                head -n 10 "$qual_file"  # Show first 10 lines of metrics
                echo ""
                
                # Check if PDF was generated
                if [ -f "quality_results/${sample}_qual.pdf" ]; then
                    echo "✓ Quality control plots generated"
                else
                    echo "⚠ Quality control plots missing"
                fi
            else
                echo "⚠ Quality control failed or incomplete"
            fi
            echo "================================"
            echo ""
        done
    } > quality_results/detailed_summary.txt
    
    echo "Detailed summary available in quality_results/detailed_summary.txt"
fi
