#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=metagene_analysis_whole_region.log

# Define variables
WORK_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
LONGEST_TRANSCRIPTS_INFO="${WORK_DIR}/longest.transcripts.info.txt"
ATTRIBUTES_FILE="${WORK_DIR}/attributes_transcriptome.txt"
SELECT_TRANS_LIST="${WORK_DIR}/total_TEs.txt"
OUTPUT_DIR="${WORK_DIR}/metagene_plots"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Change to working directory
cd ${WORK_DIR}

# Run MetageneAnalysisForTheWholeRegions using the attributes file
MetageneAnalysisForTheWholeRegions \
    -f ${ATTRIBUTES_FILE} \
    -c ${LONGEST_TRANSCRIPTS_INFO} \
    -o ${OUTPUT_DIR}/metagene_analysis \
    -b 15,90,60 \
    -S ${SELECT_TRANS_LIST} \
    -l 150 \
    -n 128 \
    -m 64 \
    -e 30 \
    --id-type=gene_id \
    --plot=yes

echo "MetageneAnalysisForTheWholeRegions completed successfully!"
echo "Output files are located in: ${OUTPUT_DIR}"
