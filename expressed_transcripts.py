import os
import re
from collections import defaultdict
from typing import List, Dict, Set, Tuple

import pysam
import pandas as pd
import numpy as np
import gffutils

def get_chromosome_from_transcript(transcript_id: str) -> str:
    """
    Extracts chromosome number from an Arabidopsis transcript ID.
    Handles special cases for organellar genomes (chloroplast and mitochondria).
    
    Args:
        transcript_id: TAIR-style transcript ID (e.g., 'AT1G05730.1')
        
    Returns:
        Chromosome identifier (e.g., '1', 'Pt', 'Mt')
    """
    # Remove 'transcript:' prefix if present
    if transcript_id.startswith('transcript:'):
        transcript_id = transcript_id[len('transcript:'):]
    
    # Extract chromosome identifier between 'AT' and 'G'
    match = re.match(r'AT([0-9CM])G', transcript_id)
    if match:
        chr_id = match.group(1)
        # Convert special chromosomes
        chr_map = {
            'C': 'Pt',  # Chloroplast -> Plastid
            'M': 'Mt'   # Mitochondria
        }
        return chr_map.get(chr_id, str(chr_id))
    return None

def calculate_unified_background_threshold(background_results: pd.DataFrame) -> float:
    """
    Calculates a unified background threshold from normalized counts across all samples.
    This represents the overall background noise level in the experiment.
    
    Args:
        background_results: DataFrame containing background analysis results
        
    Returns:
        Single Q3 threshold value representing background across all samples
    """
    # Pool all normalized counts together
    all_normalized_counts = background_results['normalized_count'].values
    
    # Calculate Q3 from pooled data
    unified_q3 = float(np.percentile(all_normalized_counts, 75))
    
    # Print detailed statistics for transparency
    print("\nBackground Threshold Analysis:")
    print(f"Total background measurements: {len(all_normalized_counts):,}")
    print(f"Sample-specific Q3 values:")
    
    for sample in background_results['sample'].unique():
        sample_data = background_results[background_results['sample'] == sample]
        sample_q3 = np.percentile(sample_data['normalized_count'], 75)
        print(f"  {sample}: {sample_q3:.2f}")
    
    print(f"\nUnified Q3 threshold (pooled data): {unified_q3:.2f}")
    
    return unified_q3

def count_maug_reads(bam_file: str,
                    db: gffutils.FeatureDB,
                    transcripts: Set[str],
                    window_size: int = 3) -> Dict[str, Dict[str, float]]:
    """
    Counts and normalizes reads at mAUG sites for each transcript.
    Uses a window around the mAUG to capture translation initiation events.
    
    Args:
        bam_file: Path to BAM file
        db: GFF database
        transcripts: Set of transcript IDs to analyze
        window_size: Size of window around mAUG (default: 3nt each side)
        
    Returns:
        Dictionary with raw and normalized counts for each transcript
    """
    raw_counts = defaultdict(int)
    normalized_counts = defaultdict(float)
    processed = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get total mapped reads for normalization
        total_reads = bam.count()
        print(f"Total mapped reads: {total_reads:,}")
        
        for transcript_id in transcripts:
            try:
                # Get chromosome and feature information
                feature = db[transcript_id]
                chromosome = get_chromosome_from_transcript(transcript_id)
                
                if not chromosome:
                    continue
                
                # Find first CDS for mAUG position
                cds_features = list(db.children(feature, featuretype='CDS'))
                if not cds_features:
                    continue
                
                first_cds = min(cds_features, key=lambda x: x.start)
                maug_pos = first_cds.start
                
                # Count reads in window around mAUG
                count = 0
                for read in bam.fetch(chromosome, 
                                    max(0, maug_pos - window_size),
                                    maug_pos + window_size):
                    count += 1
                
                # Store both raw and normalized counts
                raw_counts[transcript_id] = count
                normalized_counts[transcript_id] = (count / total_reads) * 1e6  # RPM normalization
                
                processed += 1
                if processed % 1000 == 0:
                    print(f"Processed {processed:,} transcripts...")
                
            except Exception as e:
                print(f"Error processing {transcript_id}: {str(e)}")
                continue
    
    return {'raw': raw_counts, 'normalized': normalized_counts}

def identify_expressed_transcripts(ribo_seq_files: List[str],
                                 db: gffutils.FeatureDB,
                                 all_transcripts: Set[str],
                                 background_q3: float,
                                 min_raw_count: int = 10) -> Set[str]:
    """
    Identifies transcripts with reliable translation initiation at mAUGs.
    Applies both normalized count and raw count filters across all samples.
    
    Args:
        ribo_seq_files: List of BAM file paths
        db: GFF database
        all_transcripts: Set of all transcript IDs to consider
        background_q3: Q3 threshold from background analysis
        min_raw_count: Minimum required raw read count
        
    Returns:
        Set of transcript IDs passing all filters
    """
    # Get counts for each sample
    sample_counts = {}
    for bam_file in ribo_seq_files:
        print(f"\nProcessing {os.path.basename(bam_file)}...")
        counts = count_maug_reads(bam_file, db, all_transcripts)
        sample_counts[bam_file] = counts
    
    # Find transcripts passing thresholds in all samples
    passing_transcripts = set()
    for transcript_id in all_transcripts:
        passes_all = True
        for bam_file in ribo_seq_files:
            raw_count = sample_counts[bam_file]['raw'][transcript_id]
            norm_count = sample_counts[bam_file]['normalized'][transcript_id]
            
            if raw_count < min_raw_count or norm_count < background_q3:
                passes_all = False
                break
        
        if passes_all:
            passing_transcripts.add(transcript_id)
    
    return passing_transcripts

def main():
    """
    Main function implementing the complete analysis protocol:
    1. Load background analysis results
    2. Calculate unified background threshold
    3. Count reads at mAUGs
    4. Identify expressed transcripts using both raw and normalized count filters
    """
    # File paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    
    RIBO_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT103-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT103-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-2_uniq_sort.bam")
    ]
    
    # Load background analysis results
    print("Loading background analysis results...")
    background_results = pd.read_csv(
        os.path.join(BASE_DIR, "systemPipeR/maug_upstream_analysis.tsv"),
        sep='\t'
    )
    
    # Calculate unified background threshold
    background_threshold = calculate_unified_background_threshold(background_results)
    
    # Load GFF database
    print("\nLoading GFF database...")
    db = gffutils.FeatureDB(GFF3_FILE + '.db')
    
    # Get all transcripts
    print("Getting all transcripts...")
    all_transcripts = set(feature.id for feature in db.features_of_type('mRNA'))
    print(f"Found {len(all_transcripts):,} total transcripts")
    
    # Find expressed transcripts
    expressed_transcripts = identify_expressed_transcripts(
        RIBO_SEQ_FILES,
        db,
        all_transcripts,
        background_threshold,
        min_raw_count=10
    )
    
    # Print results
    print(f"\nAnalysis Results:")
    print(f"Found {len(expressed_transcripts):,} expressed transcripts")
    print(f"These transcripts have:")
    print(f"- Normalized counts ≥ {background_threshold:.2f} (background Q3)")
    print(f"- Raw counts ≥ 10 in all samples")
    
    # Save results
    output_file = os.path.join(BASE_DIR, "systemPipeR/translating_AUGs/expressed_transcripts.txt")
    with open(output_file, 'w') as f:
        for transcript_id in sorted(expressed_transcripts):
            f.write(f"{transcript_id}\n")
    
    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    main()
