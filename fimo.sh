#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-7
#SBATCH --job-name=seq_prep
#SBATCH --output=seq_prep_%A_%a.out
#SBATCH --error=seq_prep_%A_%a.err

cd /global/scratch/users/enricocalvane/riboseq/imb2/FIMO

# Map array index to category
case $SLURM_ARRAY_TASK_ID in
    1) category="uorf_down" ;;
    2) category="uorf_up" ;;
    3) category="uorf_nc" ;;
    4) category="translatome_down" ;;
    5) category="translatome_up" ;;
    6) category="translatome_nc" ;;
    7) category="full_transcriptome" ;;
esac

# Run the script for this category
/global/home/users/enricocalvane/.conda/envs/riboseq/bin/python prepare_for_MEME-FIMO.py $category
