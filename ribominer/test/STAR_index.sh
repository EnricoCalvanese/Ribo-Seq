#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=STAR_alignment.log

# Set working directories and paths
genome_fasta="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast/Saccharomyces_genome.fa"
gtf_file="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast/Saccharomyces.gtf"

STAR --runMode genomeGenerate \
     --genomeDir . \
     --genomeFastaFiles $genome_fasta \
     --sjdbGTFfile $gtf_file \
     --runThreadN 24
