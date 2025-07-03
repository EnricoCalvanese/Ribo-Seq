#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=get_protein_coding_sequence.log

# Define variables
SIF="ribocode_ribominer_latest.sif"
TRANSCRIPTS_FA="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/prepared_transcripts/transcripts_sequence.fa"
LONGEST_TRANSCRIPTS_INFO="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/longest.transcripts.info.txt"
OUTPUT_PREFIX="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/longest"

# Run GetProteinCodingSequence
singularity exec "$SIF" \
  /root/miniconda3/bin/GetProteinCodingSequence \
  -i "$TRANSCRIPTS_FA" \
  -c "$LONGEST_TRANSCRIPTS_INFO" \
  -o "$OUTPUT_PREFIX" \
  --mode whole \
  --table 1 \
  -S

echo "GetProteinCodingSequence completed successfully"
echo "Generated files:"
echo "  ${OUTPUT_PREFIX}_aa_sequences.fa (amino acid sequences)"
echo "  ${OUTPUT_PREFIX}_transcript_sequences.fa (transcript sequences)"
echo "  ${OUTPUT_PREFIX}_cds_sequences.fa (CDS sequences)"
