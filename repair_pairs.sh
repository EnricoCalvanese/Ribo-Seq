#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

cd /global/scratch/users/enricocalvane/riboseq/imb2

# Create directory for logs
mkdir -p pair_repair_logs

# Create directory for paired output
mkdir -p paired_reads

echo "Starting paired-end read repair..."

# Process each sample pair
# LZT101-1
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT101-1_L1_Q18328W9619.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT101-1_L1_Q18328W9619.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT101-1_repair.log

# LZT101-2
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT101-2_L1_Q18329W9619.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT101-2_L1_Q18329W9619.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT101-2_repair.log

# LZT102-1
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT102-1_L5_Q18330W9620.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT102-1_L5_Q18330W9620.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT102-1_repair.log

# LZT102-2
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT102-2_L5_Q18331W9620.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT102-2_L5_Q18331W9620.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT102-2_repair.log

# LZT103-1
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT103-1_L7_Q18332W9621.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT103-1_L7_Q18332W9621.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT103-1_repair.log

# LZT103-2
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT103-2_L8_Q18333W9621.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT103-2_L8_Q18333W9621.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT103-2_repair.log

# LZT104-1
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT104-1_L8_Q18334W9622.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT104-1_L8_Q18334W9622.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT104-1_repair.log

# LZT104-2
fastq_pair \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT104-2_L7_Q18335W9622.R1.rmadaptor.clean.fq \
    /global/scratch/users/enricocalvane/riboseq/imb2/processed_fastx/LZT104-2_L7_Q18335W9622.R2.rmadaptor.clean.fq \
    2>/global/scratch/users/enricocalvane/riboseq/imb2/pair_repair_logs/LZT104-2_repair.log

# Move paired files to paired_reads directory
echo "Moving paired files to paired_reads directory..."
mv *paired.fq paired_reads/

echo "Repair process complete. Check pair_repair_logs directory for processing logs."
