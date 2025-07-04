#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=output_transcript_info.log

# Define variables
MODIFIED_GTF="/global/scratch/users/enricocalvane/riboseq/Araport11_GTF_genes_ribominer.gtf"
TRANSCRIPTS_CDS="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts/transcripts_cds.txt"
TRANSCRIPTS_FA="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts/transcripts_sequence.fa"
LONGEST_OUTPUT="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/longest.transcripts.info.txt"
ALL_OUTPUT="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/all.transcripts.info.txt"

# Run OutputTranscriptInfo
OutputTranscriptInfo \
  -c "$TRANSCRIPTS_CDS" \
  -g "$MODIFIED_GTF" \
  -f "$TRANSCRIPTS_FA" \
  -o "$LONGEST_OUTPUT" \
  -O "$ALL_OUTPUT"

echo "OutputTranscriptInfo completed successfully"
