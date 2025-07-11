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
ATTRIBUTES_FILE="${WORK_DIR}/attributes.txt"
SELECT_TRANS_LIST="${WORK_DIR}/total_TEs.txt"
OUTPUT_DIR="${WORK_DIR}/metagene_plots"
OUTPUT_PREFIX="${OUTPUT_DIR}/metagene_analysis"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Change to working directory
cd ${WORK_DIR}

# Run MetageneAnalysisForTheWholeRegions using the attributes file
echo "Running MetageneAnalysisForTheWholeRegions..."
MetageneAnalysisForTheWholeRegions \
    -f ${ATTRIBUTES_FILE} \
    -c ${LONGEST_TRANSCRIPTS_INFO} \
    -o ${OUTPUT_PREFIX} \
    -b 15,90,60 \
    -S ${SELECT_TRANS_LIST} \
    -l 150 \
    -n 128 \
    -m 64 \
    -e 30 \
    --id-type=gene_id \
    --plot=yes

echo "MetageneAnalysisForTheWholeRegions completed successfully!"

# Check if the output file exists before proceeding to plotting
DENSITY_FILE="${OUTPUT_PREFIX}_scaled_density_dataframe.txt"
if [ -f "$DENSITY_FILE" ]; then
    echo "Found density dataframe file: $DENSITY_FILE"
    echo "Running PlotMetageneAnalysisForTheWholeRegions..."
    
    # Run PlotMetageneAnalysisForTheWholeRegions
    PlotMetageneAnalysisForTheWholeRegions \
        -i ${DENSITY_FILE} \
        -o ${OUTPUT_PREFIX}_plot \
        -g WT,imb2 \
        -r "WT 1,WT 2__imb2 1,imb2 2" \
        -b 15,90,60 \
        --mode all \
        --xlabel-loc -0.4
    
    echo "PlotMetageneAnalysisForTheWholeRegions completed successfully!"
    echo "Plot files are located in: ${OUTPUT_DIR}"
else
    echo "Error: Density dataframe file not found: $DENSITY_FILE"
    echo "Cannot proceed with plotting."
    exit 1
fi

echo "All analyses completed successfully!"
echo "Output files are located in: ${OUTPUT_DIR}"
