#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=data_prep.log

set -euo pipefail

cd /global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer/yeast

prepare_transcripts -g Saccharomyces.gtf -f Saccharomyces_genome.fa -o RiboCode_annot

OutputTranscriptInfo -c RiboCode_annot/transcripts_cds.txt -g Saccharomyces.gtf -f RiboCode_annot/transcripts_sequence.fa -o longest.transcripts.info.txt -O all.transcripts.info.txt

GetProteinCodingSequence -i RiboCode_annot/transcripts_sequence.fa  -c longest.transcripts.info.txt -o <output_prefix> --mode whole --table 1
