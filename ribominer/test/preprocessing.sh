#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=preprocessing.log

cd /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast

fastqc SRR5008135.fastq -o .

cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008134.trimmed.fastq SRR5008134.fastq > SRR5008134_trimmed.log
cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008135.trimmed.fastq SRR5008135.fastq > SRR5008135_trimmed.log
cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008136.trimmed.fastq SRR5008136.fastq > SRR5008136_trimmed.log
cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008137.trimmed.fastq SRR5008137.fastq > SRR5008137_trimmed.log

fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008134.trimmed.fastq -o SRR5008134.trimmed.Qfilter.fastq > SRR5008134.Qfilter.log
fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008135.trimmed.fastq -o SRR5008135.trimmed.Qfilter.fastq > SRR5008135.Qfilter.log
fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008136.trimmed.fastq -o SRR5008136.trimmed.Qfilter.fastq > SRR5008136.Qfilter.log
fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008137.trimmed.fastq -o SRR5008137.trimmed.Qfilter.fastq > SRR5008137.Qfilter.log

