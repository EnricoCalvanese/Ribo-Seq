#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
# Create output directory for processed files
mkdir -p processed_fastx

# Process R1 files
echo "Processing R1 files..."

# LZT101-1 R1
zcat LZT101-1_L1_Q18328W9619.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT101-1_L1_Q18328W9619.R1.rmadaptor.clean.fq

# LZT101-2 R1
zcat LZT101-2_L1_Q18329W9619.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT101-2_L1_Q18329W9619.R1.rmadaptor.clean.fq

# LZT102-1 R1
zcat LZT102-1_L5_Q18330W9620.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT102-1_L5_Q18330W9620.R1.rmadaptor.clean.fq

# LZT102-2 R1
zcat LZT102-2_L5_Q18331W9620.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT102-2_L5_Q18331W9620.R1.rmadaptor.clean.fq

# LZT103-1 R1
zcat LZT103-1_L7_Q18332W9621.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT103-1_L7_Q18332W9621.R1.rmadaptor.clean.fq

# LZT103-2 R1
zcat LZT103-2_L8_Q18333W9621.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT103-2_L8_Q18333W9621.R1.rmadaptor.clean.fq

# LZT104-1 R1
zcat LZT104-1_L8_Q18334W9622.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT104-1_L8_Q18334W9622.R1.rmadaptor.clean.fq

# LZT104-2 R1
zcat LZT104-2_L7_Q18335W9622.R1.clean.fastq.gz | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c > \
processed_fastx/LZT104-2_L7_Q18335W9622.R1.rmadaptor.clean.fq

# Process R2 files
echo "Processing R2 files..."

# LZT101-1 R2
zcat LZT101-1_L1_Q18328W9619.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT101-1_L1_Q18328W9619.R2.rmadaptor.clean.fq

# LZT101-2 R2
zcat LZT101-2_L1_Q18329W9619.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT101-2_L1_Q18329W9619.R2.rmadaptor.clean.fq

# LZT102-1 R2
zcat LZT102-1_L5_Q18330W9620.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT102-1_L5_Q18330W9620.R2.rmadaptor.clean.fq

# LZT102-2 R2
zcat LZT102-2_L5_Q18331W9620.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT102-2_L5_Q18331W9620.R2.rmadaptor.clean.fq

# LZT103-1 R2
zcat LZT103-1_L7_Q18332W9621.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT103-1_L7_Q18332W9621.R2.rmadaptor.clean.fq

# LZT103-2 R2
zcat LZT103-2_L8_Q18333W9621.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT103-2_L8_Q18333W9621.R2.rmadaptor.clean.fq

# LZT104-1 R2
zcat LZT104-1_L8_Q18334W9622.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT104-1_L8_Q18334W9622.R2.rmadaptor.clean.fq

# LZT104-2 R2
zcat LZT104-2_L7_Q18335W9622.R2.clean.fastq.gz | \
fastx_clipper -Q33 -a AGATCGGAAGAGCGT -l 17 | \
fastx_reverse_complement -Q33 | \
fastx_clipper -Q33 -a CTGTAGGCACCATCA -l 17 -c | \
fastx_reverse_complement -Q33 > \
processed_fastx/LZT104-2_L7_Q18335W9622.R2.rmadaptor.clean.fq

echo "All files have been processed"
