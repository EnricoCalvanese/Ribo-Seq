#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL

mkdir /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/translating_AUGs
cd /global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/translating_AUGs

module load python

/global/home/users/enricocalvane/.conda/envs/riboseq/bin/python calculate_translating_AUGs.py
