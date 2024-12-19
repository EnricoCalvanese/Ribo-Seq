#!/bin/bash

# Script to process RNA-seq data with fastp
# Processes paired-end reads and generates cleaned fastq files

# Create a directory for logs if it doesn't exist
mkdir -p fastp_logs

# Process LZT101-1
fastp -i LZT101-1_L1_Q18328W9619.R1.fastq.gz \
      -I LZT101-1_L1_Q18328W9619.R2.fastq.gz \
      -o LZT101-1_L1_Q18328W9619.R1.clean.fastq.gz \
      -O LZT101-1_L1_Q18328W9619.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT101-1_fastp.log

# Process LZT101-2
fastp -i LZT101-2_L1_Q18329W9619.R1.fastq.gz \
      -I LZT101-2_L1_Q18329W9619.R2.fastq.gz \
      -o LZT101-2_L1_Q18329W9619.R1.clean.fastq.gz \
      -O LZT101-2_L1_Q18329W9619.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT101-2_fastp.log

# Process LZT102-1
fastp -i LZT102-1_L5_Q18330W9620.R1.fastq.gz \
      -I LZT102-1_L5_Q18330W9620.R2.fastq.gz \
      -o LZT102-1_L5_Q18330W9620.R1.clean.fastq.gz \
      -O LZT102-1_L5_Q18330W9620.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT102-1_fastp.log

# Process LZT102-2
fastp -i LZT102-2_L5_Q18331W9620.R1.fastq.gz \
      -I LZT102-2_L5_Q18331W9620.R2.fastq.gz \
      -o LZT102-2_L5_Q18331W9620.R1.clean.fastq.gz \
      -O LZT102-2_L5_Q18331W9620.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT102-2_fastp.log

# Process LZT103-1
fastp -i LZT103-1_L7_Q18332W9621.R1.fastq.gz \
      -I LZT103-1_L7_Q18332W9621.R2.fastq.gz \
      -o LZT103-1_L7_Q18332W9621.R1.clean.fastq.gz \
      -O LZT103-1_L7_Q18332W9621.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT103-1_fastp.log

# Process LZT103-2
fastp -i LZT103-2_L8_Q18333W9621.R1.fastq.gz \
      -I LZT103-2_L8_Q18333W9621.R2.fastq.gz \
      -o LZT103-2_L8_Q18333W9621.R1.clean.fastq.gz \
      -O LZT103-2_L8_Q18333W9621.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT103-2_fastp.log

# Process LZT104-1
fastp -i LZT104-1_L8_Q18334W9622.R1.fastq.gz \
      -I LZT104-1_L8_Q18334W9622.R2.fastq.gz \
      -o LZT104-1_L8_Q18334W9622.R1.clean.fastq.gz \
      -O LZT104-1_L8_Q18334W9622.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT104-1_fastp.log

# Process LZT104-2
fastp -i LZT104-2_L7_Q18335W9622.R1.fastq.gz \
      -I LZT104-2_L7_Q18335W9622.R2.fastq.gz \
      -o LZT104-2_L7_Q18335W9622.R1.clean.fastq.gz \
      -O LZT104-2_L7_Q18335W9622.R2.clean.fastq.gz \
      -A -G -n 15 2>fastp_logs/LZT104-2_fastp.log

echo "All samples have been processed with fastp"
