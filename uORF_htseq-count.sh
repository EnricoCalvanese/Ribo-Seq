#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

#riboseq must be activated
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR

# Create output directories for counts
mkdir -p uORF_counts/TEup uORF_counts/TEnc uORF_counts/TEdown

# For TEnc
for bam in /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/filtered_bams/TEnc/LZT10*_uniq_sort_TEnc.bam; do
    sample=$(basename $bam _uniq_sort_TEnc.bam)
    htseq-count -f bam -r name -s yes \
        -t sequence_feature \
        -i ID \
        -m union --nonunique none \
        $bam /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/uorf.gff > uORF_counts/TEnc/${sample}_uorf_counts.txt
done

# For TEup
for bam in /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/filtered_bams/TEup/LZT10*_uniq_sort_TEup.bam; do
    sample=$(basename $bam _uniq_sort_TEup.bam)
    htseq-count -f bam -r name -s yes \
        -t sequence_feature \
        -i ID \
        -m union --nonunique none \
        $bam /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/uorf.gff > uORF_counts/TEup/${sample}_uorf_counts.txt
done

# For TEdown
for bam in /global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/filtered_bams/TEdown/LZT10*_uniq_sort_TEdown.bam; do
    sample=$(basename $bam _uniq_sort_TEdown.bam)
    htseq-count -f bam -r name -s yes \
        -t sequence_feature \
        -i ID \
        -m union --nonunique none \
        $bam /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/uorf.gff > uORF_counts/TEdown/${sample}_uorf_counts.txt
done
