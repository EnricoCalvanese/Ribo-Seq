"""
RNA Secondary Structure Analysis Pipeline

This script analyzes RNA secondary structures in 5' leader sequences,
comparing different translational efficiency categories using ViennaRNA.
It includes robust error handling, progress tracking, and statistical analysis.

Author: Claude
Date: January 2024
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
from datetime import datetime

# Define input and output paths
INPUT_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2/FIMO/leader_sequences"
OUTPUT_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2/RNAfold"
NUM_CORES = 24  # Matching SLURM configuration

def setup_output_directory():
    """
    Create an organized output directory structure with timestamp.
    Returns the path to the newly created directory.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(OUTPUT_DIR, f"run_{timestamp}")
    
    # Create subdirectories for different output types
    subdirs = ['plots', 'data', 'statistics']
    for subdir in subdirs:
        os.makedirs(os.path.join(run_dir, subdir), exist_ok=True)
    
    return run_dir

def get_local_structure_complexity(structure, window_size=30):
    """
    Calculate structural complexity using a sliding window approach.
    This method is more stable for analyzing long sequences.
    
    Parameters:
        structure: The dot-bracket notation of the RNA structure
        window_size: Size of the sliding window for local analysis
    
    Returns:
        float: Average local structural complexity
    """
    if len(structure) <= window_size:
        return structure.count('(') / len(structure)
    
    complexities = []
    for i in range(len(structure) - window_size + 1):
        window = structure[i:i + window_size]
        complexities.append(window.count('(') / window_size)
    
    return np.mean(complexities)

def find_structural_motifs(structure):
    """
    Identify regulatory structural motifs using pattern matching.
    
    Parameters:
        structure: The dot-bracket notation of the RNA structure
    
    Returns:
        dict: Counts of different structural motifs
    """
    stem_loops = len(re.findall(r'\([.]+\)', structure))
    internal_loops = len(re.findall(r'\([.]+\([.]+\)[.]+\)', structure))
    bulges = len(re.findall(r'\([.]+\(+\)+[.]+\)', structure))
    
    total_motifs = stem_loops + internal_loops + bulges
    motif_density = total_motifs / len(structure) if len(structure) > 0 else 0
    
    return {
        'stem_loops': stem_loops,
        'internal_loops': internal_loops,
        'bulges': bulges,
        'motif_density': motif_density
    }

def analyze_single_sequence(sequence_record):
    """
    Analyze a single RNA sequence with improved handling of long sequences
    and numerical stability.
    
    Parameters:
        sequence_record: A BioPython SeqRecord object
    
    Returns:
        dict: Structural analysis results
    """
    rna_seq = str(sequence_record.seq).replace('T', 'U')
    seq_length = len(rna_seq)
    
    try:
        # Adjust scaling parameter for long sequences
        if seq_length > 1000:
            RNA.cvar.pf_scale = 2.0
        
        # Calculate MFE structure
        (structure, mfe) = RNA.fold(rna_seq)
        
        try:
            # Attempt partition function calculation
            (ensemble_struct, ensemble_energy) = RNA.pf_fold(rna_seq)
            ensemble_diversity = RNA.mean_bp_distance(seq_length)
        except Exception as e:
            print(f"Warning: Partition function calculation failed for {sequence_record.id}: {str(e)}")
            ensemble_energy = float('nan')
            ensemble_diversity = float('nan')
        
        # Basic structural analysis
        stem_count = structure.count('(')
        loop_count = structure.count('.')
        gc_content = (rna_seq.count('G') + rna_seq.count('C')) / seq_length
        local_structure_complexity = get_local_structure_complexity(structure)
        
        return {
            'sequence_id': sequence_record.id,
            'length': seq_length,
            'mfe': mfe,
            'ensemble_energy': ensemble_energy,
            'ensemble_diversity': ensemble_diversity,
            'stem_density': stem_count / seq_length,
            'loop_density': loop_count / seq_length,
            'gc_content': gc_content,
            'local_complexity': local_structure_complexity,
            **find_structural_motifs(structure)
        }
        
    except Exception as e:
        print(f"Error processing sequence {sequence_record.id}: {str(e)}")
        return None

def analyze_fasta_file_parallel(fasta_file, num_cores):
    """
    Process a FASTA file using parallel computing with progress tracking.
    
    Parameters:
        fasta_file: Path to the FASTA file
        num_cores: Number of CPU cores to use
    
    Returns:
        pandas.DataFrame: Analysis results for all sequences
    """
    # Count total sequences for progress tracking
    total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    print(f"Processing {len(sequences)} sequences from {os.path.basename(fasta_file)}")
    
    # Process sequences in parallel with progress tracking
    with Pool(num_cores) as pool:
        chunk_size = max(1, min(1000, total_sequences // (num_cores * 4)))
        results = []
        
        for i, result in enumerate(pool.imap(analyze_single_sequence, 
                                           sequences, 
                                           chunksize=chunk_size)):
            if result is not None:
                results.append(result)
            
            if (i + 1) % 100 == 0:
                print(f"Processed {i + 1}/{total_sequences} sequences "
                      f"({((i + 1)/total_sequences)*100:.1f}%)")
    
    return pd.DataFrame(results)

def generate_visualizations(dataframes, categories, output_dir):
    """
    Create comprehensive visualizations comparing structural features.
    
    Parameters:
        dataframes: List of pandas DataFrames containing analysis results
        categories: List of category names
        output_dir: Directory to save the plots
    """
    features = [col for col in dataframes[0].columns 
               if col not in ['sequence_id']]
    
    for feature in features:
        plt.figure(figsize=(12, 6))
        data = pd.concat([df.assign(Category=cat) for df, cat 
                         in zip(dataframes, categories)])
        
        sns.violinplot(x='Category', y=feature, data=data)
        plt.title(f'Distribution of {feature} across categories')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_dir, 'plots', f'{feature}_comparison.png'))
        plt.close()

def perform_statistical_analysis(dataframes, categories, output_dir):
    """
    Perform statistical analysis with improved robustness and error handling.
    
    Parameters:
        dataframes: List of pandas DataFrames containing analysis results
        categories: List of category names
        output_dir: Directory to save the results
    """
    features = [col for col in dataframes[0].columns 
               if col not in ['sequence_id']]
    
    stats_results = []
    
    for feature in features:
        try:
            valid_data = [df[feature].dropna() for df in dataframes]
            
            if all(len(data) > 0 for data in valid_data):
                if len(set(np.concatenate(valid_data))) > 1:
                    h_stat, p_val = stats.kruskal(*valid_data)
                    
                    pairwise_tests = []
                    for i, cat1 in enumerate(categories):
                        for j, cat2 in enumerate(categories[i+1:], i+1):
                            try:
                                stat, p = stats.mannwhitneyu(valid_data[i],
                                                           valid_data[j])
                                pairwise_tests.append({
                                    'comparison': f'{cat1}_vs_{cat2}',
                                    'p_value': p
                                })
                            except Exception as e:
                                print(f"Warning: Pairwise comparison failed for "
                                      f"{cat1} vs {cat2} on {feature}: {str(e)}")
                    
                    stats_results.append({
                        'feature': feature,
                        'kruskal_wallis_p': p_val,
                        'pairwise_tests': pairwise_tests
                    })
                else:
                    print(f"Warning: All values identical for feature {feature}")
            
        except Exception as e:
            print(f"Warning: Statistical analysis failed for feature {feature}: {str(e)}")
    
    # Save detailed results
    with open(os.path.join(output_dir, 'statistics', 'statistical_analysis.txt'), 'w') as f:
        f.write("Statistical Analysis Results\n")
        f.write("==========================\n\n")
        
        for result in stats_results:
            f.write(f"\nFeature: {result['feature']}\n")
            f.write("-" * (len(result['feature']) + 9) + "\n")
            f.write(f"Kruskal-Wallis p-value: {result['kruskal_wallis_p']:.6f}\n")
            f.write("\nPairwise comparisons:\n")
            for test in result['pairwise_tests']:
                f.write(f"{test['comparison']}: p = {test['p_value']:.6f}\n")

def main():
    """
    Main execution function that coordinates the analysis pipeline.
    """
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
    
    # Process each file and store results
    results = {}
    for category, filename in input_files.items():
        print(f"Processing {category}...")
        filepath = os.path.join(INPUT_DIR, filename)
        results[category] = analyze_fasta_file_parallel(filepath, NUM_CORES)
        
        # Save individual results
        results[category].to_csv(
            os.path.join(run_dir, 'data', f'{category}_analysis.csv'))
    
    # Generate visualizations
    print("Generating visualizations...")
    generate_visualizations(list(results.values()), 
                          list(results.keys()), 
                          run_dir)
    
    # Perform statistical analysis
    print("Performing statistical analysis...")
    perform_statistical_analysis(list(results.values()), 
                               list(results.keys()), 
                               run_dir)
    
    print(f"Analysis complete. Results saved in {run_dir}")

if __name__ == "__main__":
    main()
