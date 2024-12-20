#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

# Create output directories
mkdir -p filtered_reads
mkdir -p alignment_logs
mkdir -p sam_files

# Set the index path
INDEX_PATH="/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/rRNA_tRNA_index"

# Process each sample pair
for R1 in /global/scratch/users/enricocalvane/riboseq/imb2/paired_reads/*R1.rmadaptor.clean.fq.paired.fq; do
    # Get the corresponding R2 file
    R2=${R1/R1/R2}
    
    # Extract sample name for output files
    SAMPLE=$(basename $R1 | sed 's/_L[0-9]_.*$//')
    
    echo "Processing sample: $SAMPLE"
    
    # Run bowtie2 with 24 cores
    bowtie2 -p 24 \
            -x $INDEX_PATH \
            -1 $R1 \
            -2 $R2 \
            --un-conc filtered_reads/${SAMPLE}.filtered.fq \
            > sam_files/${SAMPLE}.tRNA_rRNA.sam \
            2> alignment_logs/${SAMPLE}.bowtie2.log

    echo "Completed processing: $SAMPLE"
done

echo "All samples processed"
