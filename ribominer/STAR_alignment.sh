#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=STAR_alignment_permissive.log

# Load required modules
module load bio/samtools/1.17-gcc-11.4.0

# Set working directories and paths
workdir="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
fastq_dir="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/filtered_reads"
genome_fasta="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
gtf_file="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_ribominer.gtf"
star_index_dir="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/STAR_index"

# Create necessary directories
mkdir -p $star_index_dir
mkdir -p $workdir/STAR_alignment_permissive

echo "Starting STAR genome indexing..."
echo "Time: $(date)"

# Step 1: Build STAR genome index
STAR --runMode genomeGenerate \
     --genomeDir $star_index_dir \
     --genomeFastaFiles $genome_fasta \
     --sjdbGTFfile $gtf_file \
     --runThreadN 24

echo "Genome indexing completed at: $(date)"
echo "Starting alignment process..."

# Step 2: Process paired-end FASTQ files
cd $fastq_dir

# Get unique sample names (remove .filtered.1.fq and .filtered.2.fq suffixes)
samples=$(ls *.filtered.1.fq | sed 's/.filtered.1.fq//' | sort -u)

for sample in $samples; do
    echo "Processing sample: $sample"
    echo "Time: $(date)"
    
    # Check if both paired files exist
    read1="${sample}.filtered.1.fq"
    read2="${sample}.filtered.2.fq"
    
    if [[ ! -f "$read1" || ! -f "$read2" ]]; then
        echo "Warning: Missing paired files for $sample. Skipping..."
        continue
    fi
    
    # Create output directory for this sample
    output_dir="$workdir/STAR_alignment_permissive"
    
    echo "Aligning $sample with STAR..."
    
    # Step 3: Run STAR alignment
    STAR --runThreadN 24 \
         --outFilterType Normal \
         --outWigType wiggle \
         --outWigStrand Stranded \
         --outWigNorm RPM \
         --outFilterMismatchNmax 3 \     
         --outFilterMultimapNmax 20 \     
         --outSAMmultNmax 1 \     
         --outMultimapperOrder Random \ 
         --genomeDir $star_index_dir \
         --readFilesIn $fastq_dir/$read1 $fastq_dir/$read2 \
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
echo "Results located in: $workdir/STAR_alignment_permissive/"
echo "Index located in: $star_index_dir"
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
