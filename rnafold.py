"""
RNA Secondary Structure Analysis Pipeline
This script analyzes RNA secondary structures in 5' leader sequences,
comparing different translational efficiency categories using ViennaRNA.
"""

import RNA
from Bio import SeqIO
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import re
import os
from multiprocessing import Pool
from functools import partial
import glob
from datetime import datetime

# Define input and output paths
INPUT_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2/FIMO/leader_sequences"
OUTPUT_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2/RNAfold"
NUM_CORES = 24  # Matching your SLURM configuration

def setup_output_directory():
    """Create output directory structure with timestamp"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(OUTPUT_DIR, f"run_{timestamp}")
    subdirs = ['plots', 'data', 'statistics']
    
    for subdir in subdirs:
        os.makedirs(os.path.join(run_dir, subdir), exist_ok=True)
    
    return run_dir

def analyze_single_sequence(sequence_record):
    """
    Analyze a single RNA sequence to determine its structural properties.
    
    This function calculates several important metrics:
    1. Minimum free energy (MFE) structure - the most energetically favorable fold
    2. Ensemble diversity - how much the structure varies across different possible folds
    3. Base pairing patterns - the distribution of paired and unpaired bases
    
    Parameters:
        sequence_record: A BioPython SeqRecord object containing the sequence
        
    Returns:
        dict: A dictionary containing all calculated structural metrics
    """
    # Convert DNA to RNA sequence
    rna_seq = str(sequence_record.seq).replace('T', 'U')
    seq_length = len(rna_seq)
    
    # Calculate minimum free energy structure
    # fold() returns a tuple of (structure_string, mfe)
    (structure, mfe) = RNA.fold(rna_seq)
    
    # Calculate partition function for ensemble properties
    # pf_fold() returns (ensemble_structure, ensemble_energy)
    (ensemble_struct, ensemble_energy) = RNA.pf_fold(rna_seq)
    
    # Calculate ensemble diversity using sequence length
    ensemble_diversity = RNA.mean_bp_distance(seq_length)
    
    # Analyze base pairing using the structure string
    stem_count = structure.count('(')  # Count base pairs
    loop_count = structure.count('.')  # Count unpaired bases
    
    # Calculate structure-based metrics
    pairing_ratio = stem_count / seq_length if seq_length > 0 else 0
    
    # Find structural motifs in the MFE structure
    motifs = find_structural_motifs(structure)
    
    return {
        'sequence_id': sequence_record.id,
        'length': seq_length,
        'mfe': mfe,
        'ensemble_energy': ensemble_energy,
        'ensemble_diversity': ensemble_diversity,
        'stem_density': pairing_ratio,
        'loop_density': loop_count / seq_length,
        'structural_complexity': ensemble_diversity / seq_length,
        **motifs  # Include all identified motifs
    }    
def find_structural_motifs(structure):
    """
    Identify regulatory structural motifs using pattern matching
    """
    # Initialize counters for different motif types
    stem_loops = len(re.findall(r'\([.]+\)', structure))
    internal_loops = len(re.findall(r'\([.]+\([.]+\)[.]+\)', structure))
    bulges = len(re.findall(r'\([.]+\(+\)+[.]+\)', structure))
    
    # Calculate complexity metrics
    total_motifs = stem_loops + internal_loops + bulges
    motif_density = total_motifs / len(structure) if len(structure) > 0 else 0
    
    return {
        'stem_loops': stem_loops,
        'internal_loops': internal_loops,
        'bulges': bulges,
        'motif_density': motif_density
    }

def analyze_fasta_file_parallel(fasta_file, num_cores):
    """
    Process a FASTA file using parallel computing
    """
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    with Pool(num_cores) as pool:
        results = pool.map(analyze_single_sequence, sequences)
    
    return pd.DataFrame(results)

def generate_visualizations(dataframes, categories, output_dir):
    """
    Create comprehensive visualizations comparing structural features
    """
    features = ['mfe', 'ensemble_diversity', 'stem_density', 'loop_density',
                'structural_complexity', 'motif_density', 'stem_loops',
                'internal_loops', 'bulges']
    
    for feature in features:
        plt.figure(figsize=(12, 6))
        data = pd.concat([df.assign(Category=cat) for df, cat in zip(dataframes, categories)])
        
        # Create violin plot
        sns.violinplot(x='Category', y=feature, data=data)
        plt.title(f'Distribution of {feature} across categories')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        # Save plot
        plt.savefig(os.path.join(output_dir, 'plots', f'{feature}_comparison.png'))
        plt.close()

def perform_statistical_analysis(dataframes, categories, output_dir):
    """
    Perform comprehensive statistical analysis and save results
    """
    features = dataframes[0].columns.drop(['sequence_id'])
    stats_results = []
    
    for feature in features:
        # Kruskal-Wallis test
        h_stat, p_val = stats.kruskal(*[df[feature] for df in dataframes])
        
        # Pairwise Mann-Whitney U tests
        pairwise_tests = []
        for i, cat1 in enumerate(categories):
            for j, cat2 in enumerate(categories[i+1:], i+1):
                stat, p = stats.mannwhitneyu(dataframes[i][feature],
                                           dataframes[j][feature])
                pairwise_tests.append({
                    'comparison': f'{cat1}_vs_{cat2}',
                    'p_value': p
                })
        
        stats_results.append({
            'feature': feature,
            'kruskal_wallis_p': p_val,
            'pairwise_tests': pairwise_tests
        })
    
    # Save results
    with open(os.path.join(output_dir, 'statistics', 'statistical_analysis.txt'), 'w') as f:
        for result in stats_results:
            f.write(f"\nFeature: {result['feature']}\n")
            f.write(f"Kruskal-Wallis p-value: {result['kruskal_wallis_p']}\n")
            f.write("Pairwise comparisons:\n")
            for test in result['pairwise_tests']:
                f.write(f"{test['comparison']}: p = {test['p_value']}\n")

def main():
    # Create output directory
    run_dir = setup_output_directory()
    
    # Define input files and categories
    input_files = {
        'translatome_down': 'translatome_down_leaders.fa',
        'translatome_nc': 'translatome_nc_leaders.fa',
        'translatome_up': 'translatome_up_leaders.fa',
        'uORF_down': 'uORF_down_leaders.fa',
        'uORF_nc': 'uORF_nc_leaders.fa',
        'uORF_up': 'uORF_up_leaders.fa'
    }
    
    # Process each file
    results = {}
    for category, filename in input_files.items():
        print(f"Processing {category}...")
        filepath = os.path.join(INPUT_DIR, filename)
        results[category] = analyze_fasta_file_parallel(filepath, NUM_CORES)
        
        # Save individual results
        results[category].to_csv(os.path.join(run_dir, 'data', f'{category}_analysis.csv'))
    
    # Generate visualizations
    generate_visualizations(list(results.values()), list(results.keys()), run_dir)
    
    # Perform statistical analysis
    perform_statistical_analysis(list(results.values()), list(results.keys()), run_dir)
    
    print(f"Analysis complete. Results saved in {run_dir}")

if __name__ == "__main__":
    main()
