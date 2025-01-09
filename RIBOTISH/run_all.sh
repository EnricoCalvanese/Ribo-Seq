#!/bin/bash

# Set working directory
cd /global/scratch/users/enricocalvane/riboseq/imb2/ribotish

# Create all necessary directories
mkdir -p {reference,psite_results,results/{AUG_uORFs,nonAUG_uORFs}}

# Submit 5' UTR extraction job
UTR_JOB=$(sbatch extract_utr5.sh | awk '{print $4}')

# Submit P-site determination job (depends on UTR extraction)
PSITE_JOB=$(sbatch --dependency=afterok:${UTR_JOB} psite_determination.sh | awk '{print $4}')

# Submit uORF analysis job (depends on P-site determination)
ANALYSIS_JOB=$(sbatch --dependency=afterok:${PSITE_JOB} uorf_analysis.sh | awk '{print $4}')

echo "All jobs submitted:"
echo "UTR extraction job ID: ${UTR_JOB}"
echo "P-site determination job ID: ${PSITE_JOB}"
echo "uORF analysis job ID: ${ANALYSIS_JOB}"
