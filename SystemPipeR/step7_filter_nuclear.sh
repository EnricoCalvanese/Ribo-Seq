#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1:00:00
#SBATCH --job-name=filter_nuclear
#SBATCH --output=filter_nuclear_%j.out
#SBATCH --error=filter_nuclear_%j.err

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/attempt2

# Create directory for filtered files
mkdir -p results/counts/nuclear

# Process each count file
for sample in LZT103-1 LZT103-2 LZT104-1 LZT104-2; do
    for type in uorf morf; do
        input_file="results/counts/${sample}_${type}_counts.txt"
        output_file="results/counts/nuclear/${sample}_${type}_counts.txt"
        
        echo "Processing ${sample} ${type} counts..."
        
        # Count total entries before filtering
        total=$(wc -l < "$input_file")
        
        # Filter out mitochondrial and chloroplast genes and save to new file
        grep -v 'M\|C' "$input_file" > "$output_file"
        
        # Count remaining entries
        remaining=$(wc -l < "$output_file")
        
        # Calculate removed entries
        removed=$((total - remaining))
        
        echo "${sample} ${type}: removed ${removed} entries out of ${total}"
    done
done

echo "Filtering complete. Nuclear genes are in results/counts/nuclear/"
