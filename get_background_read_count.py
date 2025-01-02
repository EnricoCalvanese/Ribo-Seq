import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
from typing import List, Dict, Set, Tuple

# Previous helper functions remain the same...
# [load_uorf_positions, get_transcript_lengths, get_eligible_transcripts remain unchanged]

def count_upstream_reads(bam_file: str,
                        transcripts: Set[str],
                        upstream_distance: int = 50) -> Dict[str, int]:
    """
    Count reads in regions upstream of mAUGs using pysam's fetch method.
    
    Args:
        bam_file: Path to BAM file
        transcripts: Set of eligible transcript IDs
        upstream_distance: Distance upstream of mAUG to analyze
        
    Returns:
        Dictionary mapping transcript IDs to read counts
    """
    counts = defaultdict(int)
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for transcript_id in transcripts:
            try:
                # Get all reads for this transcript first
                reads = list(bam.fetch(transcript_id))
                if not reads:
                    continue
                
                # Find the CDS start position (assuming first read position as reference)
                cds_start = reads[0].reference_start
                
                # Count reads in upstream region
                for read in reads:
                    # Check if read falls within our upstream window
                    if cds_start - upstream_distance <= read.reference_start < cds_start:
                        counts[transcript_id] += 1
                        
            except (ValueError, StopIteration) as e:
                print(f"Warning: Could not process transcript {transcript_id}: {str(e)}")
                continue
                
    return counts

def calculate_statistics(normalized_counts: Dict[str, float]) -> Dict[str, float]:
    """
    Calculate statistical measures for normalized counts with proper error handling.
    
    Args:
        normalized_counts: Dictionary of normalized read counts
        
    Returns:
        Dictionary containing statistical measures
    """
    stats = {}
    
    # Convert dictionary values to list for numpy operations
    counts_list = list(normalized_counts.values())
    
    if not counts_list:
        # Return None for all statistics if no data
        stats['q3'] = None
        stats['median'] = None
        stats['mean'] = None
        return stats
    
    # Calculate statistics
    try:
        stats['q3'] = float(np.percentile(counts_list, 75))
        stats['median'] = float(np.median(counts_list))
        stats['mean'] = float(np.mean(counts_list))
    except Exception as e:
        print(f"Warning: Error calculating statistics: {str(e)}")
        stats['q3'] = None
        stats['median'] = None
        stats['mean'] = None
    
    return stats

def main():
    # File paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")
    
    RIBO_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT103-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT103-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-2_uniq_sort.bam")
    ]
    
    print("Loading GFF database...")
    db = gffutils.FeatureDB(GFF3_FILE + '.db')
    
    print("Loading uORF positions...")
    uorf_positions = load_uorf_positions(UORF_GFF)
    
    print("Getting transcript leader lengths...")
    leader_lengths = get_transcript_lengths(db)
    
    print("Identifying eligible transcripts...")
    eligible_transcripts = get_eligible_transcripts(leader_lengths, uorf_positions)
    print(f"Found {len(eligible_transcripts)} eligible transcripts")
    
    # Process each Ribo-seq file
    results = []
    
    for bam_file in RIBO_SEQ_FILES:
        print(f"Processing {os.path.basename(bam_file)}...")
        
        # Get total reads for normalization
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = bam.count()
        
        # Count upstream reads
        counts = count_upstream_reads(bam_file, eligible_transcripts)
        
        # Load transcript abundances (you'll need to implement this based on your RNA-seq data)
        transcript_abundances = {tid: 1.0 for tid in eligible_transcripts}
        
        # Normalize counts
        normalized_counts = normalize_counts(counts, total_reads, transcript_abundances)
        
        # Calculate statistics
        stats = calculate_statistics(normalized_counts)
        
        # Store results
        sample_result = {
            'sample': os.path.basename(bam_file),
            'normalized_counts': normalized_counts,
        }
        # Add statistics to results
        sample_result.update(stats)
        results.append(sample_result)
        
        # Print statistics
        print(f"Statistics for {os.path.basename(bam_file)}:")
        print(f"  Q3: {stats['q3']:.2f if stats['q3'] is not None else 'N/A'}")
        print(f"  Median: {stats['median']:.2f if stats['median'] is not None else 'N/A'}")
        print(f"  Mean: {stats['mean']:.2f if stats['mean'] is not None else 'N/A'}")
    
    # Save results
    output_file = os.path.join(BASE_DIR, "systemPipeR/maug_upstream_analysis.tsv")
    
    # Prepare results for output
    output_data = []
    for result in results:
        for transcript_id, count in result['normalized_counts'].items():
            output_data.append({
                'sample': result['sample'],
                'transcript_id': transcript_id,
                'normalized_count': count,
                'q3': result['q3'],
                'median': result['median'],
                'mean': result['mean']
            })
    
    # Create DataFrame and save results
    if output_data:
        pd.DataFrame(output_data).to_csv(output_file, sep='\t', index=False)
        print(f"Results saved to {output_file}")
    else:
        print("Warning: No data to save")

if __name__ == "__main__":
    main()
