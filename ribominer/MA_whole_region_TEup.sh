#!/bin/bash
#SBATCH --account=fc_rnaseq
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:45:00
#SBATCH --mail-user=enrico_calvane@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --output=metagene_analysis_TEup_genes.log

# Define variables
WORK_DIR="/global/scratch/users/enricocalvane/riboseq/metagene_plot_ribominer"
LONGEST_TRANSCRIPTS_INFO="${WORK_DIR}/longest.transcripts.info.txt"
ATTRIBUTES_FILE="${WORK_DIR}/attributes.txt"
SELECT_TRANS_LIST="${WORK_DIR}/TEup0.5_genes.txt"
OUTPUT_DIR="${WORK_DIR}/metagene_plots"
OUTPUT_PREFIX="${OUTPUT_DIR}/TE_up"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p ${WORK_DIR}/metagene_plots/common_genes

# Change to working directory
cd ${WORK_DIR}/metagene_plots

#Diagnostic ribominer run

MetageneAnalysisForTheWholeRegions \
    -f ${ATTRIBUTES_FILE} \
    -c ${LONGEST_TRANSCRIPTS_INFO} \
    -o ${OUTPUT_PREFIX} \
    -S ${SELECT_TRANS_LIST} \
    -b 15,90,60 \
    -l 30 \
    -n 1 \
    -m 1 \
    -e 0 \
    --id-type=gene_id \
    --plot=no

# Extract kept transcripts from each diagnostic output
awk '{print $1}' TE_up_WT-1_bins_density.txt > common_genes/TE_up_WT-1_transcripts.txt
awk '{print $1}' TE_up_WT-2_bins_density.txt > common_genes/TE_up_WT-2_transcripts.txt
awk '{print $1}' TE_up_imb2-1_bins_density.txt > common_genes/TE_up_imb2-1_transcripts.txt
awk '{print $1}' TE_up_imb2-2_bins_density.txt > common_genes/TE_up_imb2-2_transcripts.txt

cd ${WORK_DIR}/metagene_plots/common_genes

grep -Fxf TE_up_WT-1_transcripts.txt TE_up_WT-2_transcripts.txt > TE_up_WT_common.txt
grep -Fxf TE_up_imb2-1_transcripts.txt TE_up_imb2-2_transcripts.txt > TE_up_imb2_common.txt
grep -Fxf TE_up_WT_common.txt TE_up_imb2_common.txt > TE_up_common.txt

SELECT_TRANS_LIST="${WORK_DIR}/metagene_plots/common_genes/TE_up_common.txt"


# 1. Run MetageneAnalysisForTheWholeRegions
echo "1. Running MetageneAnalysisForTheWholeRegions..."
MetageneAnalysisForTheWholeRegions \
    -f ${ATTRIBUTES_FILE} \
    -c ${LONGEST_TRANSCRIPTS_INFO} \
    -o ${OUTPUT_PREFIX} \
    -S ${SELECT_TRANS_LIST} \
    -b 60,90,15 \
    -l 3 \
    -n 1 \
    -m 1 \
    -e 0 \
    --plot=yes

echo "MetageneAnalysisForTheWholeRegions completed successfully!"

# 2. Plot MetageneAnalysisForTheWholeRegions
DENSITY_FILE="${OUTPUT_PREFIX}_scaled_density_dataframe.txt"
if [ -f "$DENSITY_FILE" ]; then
    echo "2. Running PlotMetageneAnalysisForTheWholeRegions..."
    
    PlotMetageneAnalysisForTheWholeRegions \
        -i ${DENSITY_FILE} \
        -o ${OUTPUT_PREFIX}_plot \
        -g WT,imb2 \
        -r "WT-1,WT-2__imb2-1,imb2-2" \
        -b 60,90,15 \
        --mode all \
        --xlabel-loc -0.4
    
    echo "PlotMetageneAnalysisForTheWholeRegions completed successfully!"
else
    echo "Error: Density dataframe file not found: $DENSITY_FILE"
    echo "Cannot proceed with plotting whole regions."
fi

# 3. Polarity Calculation
echo "3. Running PolarityCalculation..."
PolarityCalculation \
    -f ${ATTRIBUTES_FILE} \
    -c ${LONGEST_TRANSCRIPTS_INFO} \
    -o ${OUTPUT_PREFIX} \
    -S ${SELECT_TRANS_LIST} \
    --id-type=gene_id \
    -n 0

echo "PolarityCalculation completed successfully!"

# 4. Plot Polarity
POLARITY_FILE="${OUTPUT_PREFIX}_polarity_dataframe.txt"
if [ -f "$POLARITY_FILE" ]; then
    echo "4. Running PlotPolarity..."
    
    PlotPolarity \
        -i ${POLARITY_FILE} \
        -o ${OUTPUT_PREFIX}_polarity_plot \
        -g WT,imb2 \
        -r "WT-1,WT-2__imb2-1,imb2-2" \
        --mode mean \
        -y 5
    
    echo "PlotPolarity completed successfully!"
else
    echo "Error: Polarity dataframe file not found: $POLARITY_FILE"
    echo "Cannot proceed with polarity plotting."
fi

# 5. Metagene Analysis for CDS
echo "5. Running MetageneAnalysis for CDS..."
MetageneAnalysis \
    -f ${ATTRIBUTES_FILE} \
    -c ${LONGEST_TRANSCRIPTS_INFO} \
    -o ${OUTPUT_PREFIX}_CDS \
    -S ${SELECT_TRANS_LIST} \
    -U codon \
    -M RPKM \
    -u 0 \
    -d 500 \
    -l 3 \
    -n 1 \
    -m 1 \
    -e 0 \
    --norm yes \
    -y 100 \
    --CI 0.95 \
    --type CDS

echo "MetageneAnalysis for CDS completed successfully!"

# 6. Plot Metagene Analysis for CDS
CDS_FILE="${OUTPUT_PREFIX}_CDS_dataframe.txt"
if [ -f "$CDS_FILE" ]; then
    echo "6. Running PlotMetageneAnalysis for CDS..."
    
    PlotMetageneAnalysis \
        -i ${CDS_FILE} \
        -o ${OUTPUT_PREFIX}_CDS_grouped_plot \
        -u 0 \
        -d 500 \
        -g WT,imb2 \
        -r "WT-1,WT-2__imb2-1,imb2-2" \
        -U codon \
        --CI 0.95 \
        --mode mean
    
    echo "PlotMetageneAnalysis for CDS completed successfully!"
else
    echo "Cannot find CDS dataframe file for plotting."
fi
