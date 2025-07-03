#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=03:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=indexing_alignment_hisat.log

cd /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads

mkdir unique_reads

cd unique_reads

module load bio/samtools/1.17-gcc-11.4.0

# Process LZT101-1
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT101-1.sam > H_LZT101-1.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT101-1.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT101-1_uniq.sam
cat H_LZT101-1.txt LZT101-1_uniq.sam > LZT101-1_uniq1.sam
mv LZT101-1_uniq1.sam LZT101-1_uniq.sam
samtools view -bS LZT101-1_uniq.sam > LZT101-1_uniq.bam

# Process LZT101-2
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT101-2.sam > H_LZT101-2.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT101-2.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT101-2_uniq.sam
cat H_LZT101-2.txt LZT101-2_uniq.sam > LZT101-2_uniq1.sam
mv LZT101-2_uniq1.sam LZT101-2_uniq.sam
samtools view -bS LZT101-2_uniq.sam > LZT101-2_uniq.bam

# Process LZT102-1
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT102-1.sam > H_LZT102-1.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT102-1.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT102-1_uniq.sam
cat H_LZT102-1.txt LZT102-1_uniq.sam > LZT102-1_uniq1.sam
mv LZT102-1_uniq1.sam LZT102-1_uniq.sam
samtools view -bS LZT102-1_uniq.sam > LZT102-1_uniq.bam

# Process LZT102-2
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT102-2.sam > H_LZT102-2.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT102-2.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT102-2_uniq.sam
cat H_LZT102-2.txt LZT102-2_uniq.sam > LZT102-2_uniq1.sam
mv LZT102-2_uniq1.sam LZT102-2_uniq.sam
samtools view -bS LZT102-2_uniq.sam > LZT102-2_uniq.bam

# Process LZT103-1
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT103-1.sam > H_LZT103-1.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT103-1.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT103-1_uniq.sam
cat H_LZT103-1.txt LZT103-1_uniq.sam > LZT103-1_uniq1.sam
mv LZT103-1_uniq1.sam LZT103-1_uniq.sam
samtools view -bS LZT103-1_uniq.sam > LZT103-1_uniq.bam

# Process LZT103-2
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT103-2.sam > H_LZT103-2.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT103-2.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT103-2_uniq.sam
cat H_LZT103-2.txt LZT103-2_uniq.sam > LZT103-2_uniq1.sam
mv LZT103-2_uniq1.sam LZT103-2_uniq.sam
samtools view -bS LZT103-2_uniq.sam > LZT103-2_uniq.bam

# Process LZT104-1
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT104-1.sam > H_LZT104-1.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT104-1.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT104-1_uniq.sam
cat H_LZT104-1.txt LZT104-1_uniq.sam > LZT104-1_uniq1.sam
mv LZT104-1_uniq1.sam LZT104-1_uniq.sam
samtools view -bS LZT104-1_uniq.sam > LZT104-1_uniq.bam

# Process LZT104-2
samtools view -H /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT104-2.sam > H_LZT104-2.txt
samtools view -S /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/reads/hisat_alignment/LZT104-2.sam | grep -P "\tNH:i:1\t|\tNH:i:1$" | grep "YT:Z:CP" > LZT104-2_uniq.sam
cat H_LZT104-2.txt LZT104-2_uniq.sam > LZT104-2_uniq1.sam
mv LZT104-2_uniq1.sam LZT104-2_uniq.sam
samtools view -bS LZT104-2_uniq.sam > LZT104-2_uniq.bam
