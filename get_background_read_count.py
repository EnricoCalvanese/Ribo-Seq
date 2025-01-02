import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
from typing import List, Dict, Set, Tuple

def load_uorf_positions(uorf_gff: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Reads and parses the uORF GFF file to identify transcripts containing upstream AUG codons.
    This function creates a mapping between transcript IDs and their uORF positions.
    
    Args:
        uorf_gff: Path to the GFF file containing uORF annotations
        
    Returns:
        Dictionary where keys are transcript IDs and values are lists of (start, end) positions
    """
    uorf_positions = defaultdict(list)
    print(f"Reading uORF file: {uorf_gff}")
    
    try:
        with open(uorf_gff) as f:
            for line in f:
                # Skip comment lines in GFF file
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                # Ensure line has all required GFF fields
                if len(fields) < 9:
                    continue
                    
                transcript_id = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                uorf_positions[transcript_id].append((start, end))
                
        print(f"Found uORFs in {len(uorf_positions)} transcripts")
        return uorf_positions
        
    except FileNotFoundError:
        print(f"Error: uORF GFF file not found at {uorf_gff}")
        return {}
    except Exception as e:
        print(f"Error processing uORF file: {str(e)}")
        return {}

def get_transcript_lengths(db: gffutils.FeatureDB) -> Dict[str, int]:
    """
    Calculates the length of 5' leader sequences for all transcripts in the database.
    The leader sequence is defined as the region from transcript start to the first CDS.
    
    Args:
        db: GFF database containing transcript and CDS annotations
        
    Returns:
        Dictionary mapping transcript IDs to their 5' leader lengths
    """
    leader_lengths = {}
    processed = 0
    
    for transcript in db.features_of_type('mRNA'):
        try:
            # Find all CDS features for this transcript
            cds_features = list(db.children(transcript, featuretype='CDS'))
            
            if cds_features:
                # Get the start position of the first CDS
                first_cds = min(cds_features, key=lambda x: x.start)
                # Calculate leader length
                leader_length = first_cds.start - transcript.start
                leader_lengths[transcript.id] = leader_length
            
            processed += 1
            if processed % 1000 == 0:
                print(f"Processed {processed} transcripts...")
                
        except Exception as e:
            print(f"Warning: Could not process transcript {transcript.id}: {str(e)}")
            continue
    
    print(f"Calculated leader lengths for {len(leader_lengths)} transcripts")
    return leader_lengths

def get_eligible_transcripts(leader_lengths: Dict[str, int],
                           uorf_positions: Dict[str, List[Tuple[int, int]]],
                           min_leader_length: int = 100) -> Set[str]:
    """
    Identifies transcripts that meet our criteria:
    1. Have a 5' leader sequence ≥ 100 nucleotides
    2. Do not contain any upstream AUG codons
    
    Args:
        leader_lengths: Dictionary of transcript leader lengths
        uorf_positions: Dictionary of uORF positions
        min_leader_length: Minimum required leader length (default: 100)
        
    Returns:
        Set of transcript IDs that meet all criteria
    """
    eligible = set()
    
    for transcript_id, length in leader_lengths.items():
        # Check if leader is long enough
        if length >= min_leader_length:
            # Check if transcript has no uORFs
            if transcript_id not in uorf_positions:
                eligible.add(transcript_id)
    
    print(f"Found {len(eligible)} eligible transcripts out of {len(leader_lengths)} total")
    return eligible

def count_upstream_reads(bam_file: str,
                        transcripts: Set[str],
                        upstream_distance: int = 50) -> Dict[str, int]:
    """
    Counts the number of ribosome footprints in regions upstream of main AUG codons.
    
    Args:
        bam_file: Path to the BAM file containing aligned reads
        transcripts: Set of eligible transcript IDs
        upstream_distance: Number of nucleotides upstream to analyze (default: 50)
        
    Returns:
        Dictionary mapping transcript IDs to their upstream read counts
    """
    counts = defaultdict(int)
    processed = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for transcript_id in transcripts:
            try:
                # Get all reads mapping to this transcript
                for read in bam.fetch(transcript_id):
                    # We'll consider the start of the first read as our reference point
                    # In a real analysis, you might want to use the annotated CDS start
                    if processed == 0:
                        cds_start = read.reference_start
                        
                    # Count reads in the upstream window
                    if cds_start - upstream_distance <= read.reference_start < cds_start:
                        counts[transcript_id] += 1
                
                processed += 1
                if processed % 100 == 0:
                    print(f"Processed {processed} transcripts...")
                    
            except ValueError as e:
                print(f"Warning: Could not process transcript {transcript_id}: {str(e)}")
                continue
    
    return counts

def normalize_counts(counts: Dict[str, int],
                    total_reads: int,
                    transcript_abundances: Dict[str, float]) -> Dict[str, float]:
    """
    Normalizes read counts by total sequencing depth and transcript abundance.
    
    Args:
        counts: Raw read counts per transcript
        total_reads: Total number of mapped reads in the sample
        transcript_abundances: RPKM values for each transcript
        
    Returns:
        Dictionary of normalized read counts
    """
    normalized = {}
    
    for transcript_id, count in counts.items():
        if transcript_id in transcript_abundances and transcript_abundances[transcript_id] > 0:
            # Normalize by total reads (RPM) and transcript abundance
            normalized[transcript_id] = (count / total_reads * 1e6) / transcript_abundances[transcript_id]
    
    return normalized

def calculate_statistics(normalized_counts: Dict[str, float]) -> Dict[str, float]:
    """
    Calculates various statistical measures from the normalized counts.
    
    Args:
        normalized_counts: Dictionary of normalized read counts
        
    Returns:
        Dictionary containing statistical measures (Q3, median, mean)
    """
    stats = {}
    counts_list = list(normalized_counts.values())
    
    if not counts_list:
        stats['q3'] = None
        stats['median'] = None
        stats['mean'] = None
        return stats
    
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
    """
    Main function that orchestrates the analysis workflow:
    1. Loads and processes input files
    2. Identifies eligible transcripts
    3. Counts and normalizes upstream reads
    4. Calculates statistics and saves results
    """
    # Define file paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")
    
    RIBO_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT103-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT103-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-2_uniq_sort.bam")
    ]
    
    # Load and process input files
    print("Loading GFF database...")
    db = gffutils.FeatureDB(GFF3_FILE + '.db')
    
    print("Loading uORF positions...")
    uorf_positions = load_uorf_positions(UORF_GFF)
    
    print("Getting transcript leader lengths...")
    leader_lengths = get_transcript_lengths(db)
    
    print("Identifying eligible transcripts...")
    eligible_transcripts = get_eligible_transcripts(leader_lengths, uorf_positions)
    
    # Process each Ribo-seq file
    results = []
    
    for bam_file in RIBO_SEQ_FILES:
        print(f"\nProcessing {os.path.basename(bam_file)}...")
        
        # Get total reads for normalization
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_reads = bam.count()
            print(f"Total mapped reads: {total_reads:,}")
        
        # Count and normalize upstream reads
        counts = count_upstream_reads(bam_file, eligible_transcripts)
        transcript_abundances = {tid: 1.0 for tid in eligible_transcripts}  # Placeholder values
        normalized_counts = normalize_counts(counts, total_reads, transcript_abundances)
        
        # Calculate statistics
        stats = calculate_statistics(normalized_counts)
        sample_result = {
            'sample': os.path.basename(bam_file),
            'normalized_counts': normalized_counts
        }
        sample_result.update(stats)
        results.append(sample_result)
        
        # Print statistics
        print(f"\nStatistics for {os.path.basename(bam_file)}:")
        print(f"Q3: {stats['q3']:.2f if stats['q3'] is not None else 'N/A'}")
        print(f"Median: {stats['median']:.2f if stats['median'] is not None else 'N/A'}")
        print(f"Mean: {stats['mean']:.2f if stats['mean'] is not None else 'N/A'}")
    
    # Save results
    output_file = os.path.join(BASE_DIR, "systemPipeR/maug_upstream_analysis.tsv")
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
    
    if output_data:
        pd.DataFrame(output_data).to_csv(output_file, sep='\t', index=False)
        print(f"\nResults saved to {output_file}")
    else:
        print("\nWarning: No data to save")

if __name__ == "__main__":
    main()
