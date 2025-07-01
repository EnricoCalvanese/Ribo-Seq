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
GTF_ORIG="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_transposons.current.gtf"
FASTA="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
GTF_CLEANED="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/Araport11_GTF_cleaned.gtf"
GTF_PATCHED="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/Araport11_GTF_genes_transposons.patched.gtf"
OUTDIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts"
SIF="ribocode_ribominer_latest.sif"

mkdir -p "$OUTDIR"

echo "Cleaning GTF formatting..."
sed 's/;[[:space:]]*/;/g' "$GTF_ORIG" > "$GTF_CLEANED"

echo "Patching GTF with transcript_biotype..."
awk 'BEGIN{OFS="\t"} {
  attr = "";
  for (i=9; i<=NF; i++) attr = attr $i " ";
  tid = ""; gid = "";
  match(attr, /transcript_id[ ]*"([^"]+)"/, m1);
  match(attr, /gene_id[ ]*"([^"]+)"/, m2);
  if (m1[1] != "") tid = m1[1];
  if (m2[1] != "") gid = m2[1];
  $9 = "transcript_id \"" tid "\"; gene_id \"" gid "\"; transcript_biotype \"protein_coding\";";
  print $1, $2, $3, $4, $5, $6, $7, $8, $9;
}' "$GTF_CLEANED" > "$GTF_PATCHED"

sed 's/\tmRNA\t/\ttranscript\t/g' "$GTF_PATCHED" > "${GTF_PATCHED%.gtf}.transcript.gtf"

echo "Running prepare_transcripts with Singularity..."
singularity exec "$SIF" \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$GTF_PATCHED" \
  -f "$FASTA" \
  -o "$OUTDIR"

echo "prepare_transcripts completed successfully."
date
