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

module load bio/fastqc/0.12.1-gcc-11.4.0

fastqc SRR5008135.fastq -o .

cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008134.trimmed.fastq SRR5008134.fastq > SRR5008134_trimmed.log
cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008135.trimmed.fastq SRR5008135.fastq > SRR5008135_trimmed.log
cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008136.trimmed.fastq SRR5008136.fastq > SRR5008136_trimmed.log
cutadapt -m 15 -M 35 --match-read-wildcards -a CTGTAGGCACCATCAAT -o SRR5008137.trimmed.fastq SRR5008137.fastq > SRR5008137_trimmed.log

fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008134.trimmed.fastq -o SRR5008134.trimmed.Qfilter.fastq > SRR5008134.Qfilter.log
fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008135.trimmed.fastq -o SRR5008135.trimmed.Qfilter.fastq > SRR5008135.Qfilter.log
fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008136.trimmed.fastq -o SRR5008136.trimmed.Qfilter.fastq > SRR5008136.Qfilter.log
fastq_quality_filter -Q33 -v -q 25 -p 75 -i SRR5008137.trimmed.fastq -o SRR5008137.trimmed.Qfilter.fastq > SRR5008137.Qfilter.log

# Load Bowtie2 module
module load bio/bowtie2/2.5.1-gcc-11.4.0

# Use 24 threads and specify rRNA index
THREADS=24
INDEX="rRNA_index"

# Loop over all filtered FASTQ files
for fq in *.trimmed.Qfilter.fastq; do
    # Extract sample name prefix
    sample=$(basename "$fq" .trimmed.Qfilter.fastq)

    # Output files
    unaligned="${sample}.nonrRNA.fastq"
    aligned="${sample}.rRNA.sam"

    echo "Processing $sample ..."

    bowtie2 \
        -x $INDEX \
        -U "$fq" \
        -S "$aligned" \
        --un "$unaligned" \
        --norc \
        -p $THREADS \
        --all \
        --score-min "L,0,-0.15" \
        -L 15 \
        -D 20 \
        -R 3 \
        -N 0 \
        -i S,1,0.50

    echo "Finished $sample"
done
