#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=meme_fimo
#SBATCH --output=meme_fimo_%A_%a.out
#SBATCH --error=meme_fimo_%A_%a.err
#SBATCH --array=1-6


# Set base directory
BASE_DIR="/global/scratch/users/enricocalvane/riboseq/imb2/FIMO"
cd $BASE_DIR

# First, run MEME on TEdown uORF leaders (this runs only once, not in parallel)
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    echo "Running MEME analysis..."
    meme leader_sequences/uORF_down_leaders.fa -dna -nmotifs 5 -oc meme_output
fi

# Wait for MEME to complete before starting FIMO
sleep 30

# Map array index to dataset for FIMO analysis
case $SLURM_ARRAY_TASK_ID in
    1) dataset="uorf_up" ;;
    2) dataset="uorf_nc" ;;
    3) dataset="translatome_down" ;;
    4) dataset="translatome_up" ;;
    5) dataset="translatome_nc" ;;
    6) dataset="full_transcriptome" ;;
esac

# Convert dataset name to corresponding input file name
if [[ $dataset == uorf_* ]]; then
    input_file="leader_sequences/uORF_${dataset#uorf_}_leaders.fa"
elif [[ $dataset == "full_transcriptome" ]]; then
    input_file="leader_sequences/full_transcriptome_leaders.fa"
else
    input_file="leader_sequences/translatome_${dataset#translatome_}_leaders.fa"
fi

# Run FIMO analysis
echo "Running FIMO analysis for $dataset..."
fimo --oc fimo_output_${dataset} meme_output/meme.txt $input_file

echo "Completed FIMO analysis for $dataset"
