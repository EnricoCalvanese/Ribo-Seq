#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=STAR_alignment.log

# Load required modules
module load bio/samtools/1.17-gcc-11.4.0

# Set working directories and paths
workdir="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast"
fastq_dir="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast"
star_index_dir="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast"

echo "Starting alignment process..."
echo "Time: $(date)"

# Step 2: Process single-end FASTQ files
cd $fastq_dir

# Get yeast sample names (SRR5008134 through SRR5008137)
samples=$(ls SRR500813[4-7].nonrRNA.fastq 2>/dev/null | sed 's/.nonrRNA.fastq//' | sort)

if [[ -z "$samples" ]]; then
    echo "Error: No matching fastq files found (SRR5008134-SRR5008137.nonrRNA.fastq)"
    exit 1
fi

for sample in $samples; do
    echo "Processing sample: $sample"
    echo "Time: $(date)"
    
    # Check if fastq file exists
    fastq_file="${sample}.nonrRNA.fastq"
    
    if [[ ! -f "$fastq_file" ]]; then
        echo "Warning: Missing file $fastq_file. Skipping..."
        continue
    fi
    
    # Output directory is the same as working directory
    output_dir="$workdir"
    
    echo "Aligning $sample with STAR..."
    
    # Step 3: Run STAR alignment (single-end)
    STAR --runThreadN 24 \
         --outFilterType Normal \
         --outWigType wiggle \
         --outWigStrand Stranded \
         --outWigNorm RPM \
         --alignEndsType EndToEnd \
         --outFilterMismatchNmax 1 \
         --outFilterMultimapNmax 1 \
         --genomeDir $star_index_dir \
         --readFilesIn $fastq_dir/$fastq_file \
         --outFileNamePrefix $output_dir/${sample}. \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMattributes All
    
    echo "Alignment completed for $sample at: $(date)"
    echo "Starting post-processing..."
    
    # Step 4: Sort and index transcriptome-mapped BAM
    if [[ -f "$output_dir/${sample}.Aligned.toTranscriptome.out.bam" ]]; then
        echo "Sorting transcriptome BAM for $sample..."
        samtools sort -T $output_dir/${sample}.Aligned.toTranscriptome.out.sorted \
                      -o $output_dir/${sample}.Aligned.toTranscriptome.out.sorted.bam \
                      $output_dir/${sample}.Aligned.toTranscriptome.out.bam
        
        echo "Indexing transcriptome BAM for $sample..."
        samtools index $output_dir/${sample}.Aligned.toTranscriptome.out.sorted.bam
        
        # Remove unsorted transcriptome BAM to save space
        rm $output_dir/${sample}.Aligned.toTranscriptome.out.bam
    fi
    
    # Step 5: Index genome-mapped BAM (already sorted by STAR due to SortedByCoordinate)
    if [[ -f "$output_dir/${sample}.Aligned.sortedByCoord.out.bam" ]]; then
        echo "Indexing genome BAM for $sample..."
        samtools index $output_dir/${sample}.Aligned.sortedByCoord.out.bam
    fi
    
    # Step 6: Print alignment statistics
    echo "=== Alignment Statistics for $sample ==="
    if [[ -f "$output_dir/${sample}.Log.final.out" ]]; then
        echo "STAR alignment summary:"
        grep -E "Number of input reads|Uniquely mapped reads|Number of reads mapped to multiple loci|Number of reads mapped to too many loci|% of reads mapped to multiple loci|% of reads mapped to too many loci|% of reads unmapped" $output_dir/${sample}.Log.final.out
    fi
    echo "=========================================="
    
    echo "Sample $sample processing completed at: $(date)"
    echo ""
done

echo "All samples processed successfully!"
echo "Final completion time: $(date)"

# Print overall summary
echo ""
echo "=== OVERALL SUMMARY ==="
echo "Total samples processed: $(echo $samples | wc -w)"
echo "Results located in: $workdir/"
echo "STAR index files located in: $star_index_dir"
echo ""
echo "Output files for each sample:"
echo "  - {sample}.Aligned.sortedByCoord.out.bam (genome-mapped reads)"
echo "  - {sample}.Aligned.sortedByCoord.out.bam.bai (genome BAM index)"
echo "  - {sample}.Aligned.toTranscriptome.out.sorted.bam (transcriptome-mapped reads)"
echo "  - {sample}.Aligned.toTranscriptome.out.sorted.bam.bai (transcriptome BAM index)"
echo "  - {sample}.ReadsPerGene.out.tab (gene counts)"
echo "  - {sample}.Log.final.out (alignment statistics)"
echo "  - {sample}.Signal.Unique.str1.out.wig (strand 1 wiggle)"
echo "  - {sample}.Signal.Unique.str2.out.wig (strand 2 wiggle)"
echo "======================="
