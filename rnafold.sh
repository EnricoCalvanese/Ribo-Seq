#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=RNA_fold
#SBATCH --output=RNA_fold_%j.out
#SBATCH --error=RNA_fold_%j.err

# Run the analysis script
/global/home/users/enricocalvane/.conda/envs/riboseq/bin/python rnafold.py
