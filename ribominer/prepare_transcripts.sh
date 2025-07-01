#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=prepare_transcripts.out
#SBATCH --error=prepare_transcripts.err

set -euo pipefail

echo "Starting prepare_transcripts job on $(hostname)"
date

# Paths
FASTA="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
GTF_PATCHED="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_transposons.patched.gtf"
OUTDIR="/global/scratch/users/enricocalvane/riboseq/prepared_transcripts"
SIF="ribocode_ribominer_latest.sif"

mkdir -p "$OUTDIR"

echo "Running prepare_transcripts with Singularity..."
singularity exec "$SIF" \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$GTF_PATCHED" \
  -f "$FASTA" \
  -o "$OUTDIR"

echo "prepare_transcripts completed successfully."
date
