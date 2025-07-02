#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=prepare_transcripts.log

set -euo pipefail

echo "Starting prepare_transcripts job on $(hostname)"
date

# Paths
GTF_ORIG="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_Araport11.gene.gtf"
GTF_PATCHED="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_Araport11.gene.patched.gtf"
FASTA="/global/scratch/users/enricocalvane/riboseq/Athaliana_447_TAIR10.fa"
OUTDIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts"
SIF="ribocode_ribominer_latest.sif"

mkdir -p "$OUTDIR"

awk '
BEGIN { OFS = "\t" }
($3 == "transcript" || $3 == "exon" || $3 == "CDS") {
  tid = ""; gid = ""; rest = "";
  match($9, /transcript_id "([^"]+)"/, t);
  match($9, /gene_id "([^"]+)"/, g);
  if (t[1] && g[1]) {
    tid = t[1]; gid = g[1];
    $9 = "transcript_id \"" tid "\"; gene_id \"" gid "\"; transcript_biotype \"protein_coding\";";
    print;
  }
}
' "$GTF_ORIG" > "$GTF_PATCHED"

echo "Running prepare_transcripts with Singularity..."
singularity exec "$SIF" \
  /root/miniconda3/bin/prepare_transcripts \
  -g "$GTF_PATCHED" \
  -f "$FASTA" \
  -o "$OUTDIR"

echo "prepare_transcripts completed successfully."
date
