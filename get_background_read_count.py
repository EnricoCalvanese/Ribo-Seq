import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
from typing import List, Dict, Set, Tuple
from statistics import quantiles

def load_uorf_positions(uorf_gff: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Load uORF positions from the GFF file to identify transcripts with uAUGs.
    
    Args:
        uorf_gff: Path to the uORF GFF file
        
    Returns:
        Dictionary mapping transcript IDs to lists of uORF positions
    """
    uorf_positions = defaultdict(list)
    
    with open(uorf_gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:  # Skip malformed lines
                continue
                
            transcript_id = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            uorf_positions[transcript_id].append((start, end))
    
    return uorf_positions

def get_transcript_lengths(db: gffutils.FeatureDB) -> Dict[str, int]:
    """
    Get the length of 5' leader sequences for all transcripts.
    
    Args:
        db: GFF database
        
    Returns:
        Dictionary mapping transcript IDs to leader lengths
    """
    leader_lengths = {}
    
    for transcript in db.features_of_type('mRNA'):
        try:
            # Get CDS start position
            cds = next(db.children(transcript, featuretype='CDS'))
            # Calculate leader length (assuming transcript coordinates start at 1)
            leader_length = cds.start - transcript.start
            leader_lengths[transcript.id] = leader_length
        except StopIteration:
            continue  # Skip transcripts without CDS
            
    return leader_lengths

def get_eligible_transcripts(leader_lengths: Dict[str, int],
                           uorf_positions: Dict[str, List[Tuple[int, int]]],
                           min_leader_length: int = 100) -> Set[str]:
    """
    Identify transcripts with long enough leaders and no uAUGs.
    
    Args:
        leader_lengths: Dictionary of transcript leader lengths
        uorf_positions: Dictionary of uORF positions
        min_leader_length: Minimum required leader length
        
    Returns:
        Set of eligible transcript IDs
    """
    eligible = set()
    
    for transcript_id, length in leader_lengths.items():
        # Check leader length
        if length < min_leader_length:
            continue
            
        # Check for absence of uORFs
        if transcript_id not in uorf_positions:
            eligible.add(transcript_id)
            
    return eligible

def count_upstream_reads(bam_file: str,
                        transcripts: Set[str],
                        upstream_distance: int = 50) -> Dict[str, int]:
    """
    Count reads in regions upstream of mAUGs.

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
            print(f"Checking transcript: {transcript_id}")  # Debugging statement
            try:
                read_found = False  # Flag to check if any reads are found
                for read in bam.fetch(transcript_id):
                    if read.reference_start < upstream_distance:
                        counts[transcript_id] += 1
                        read_found = True
                if not read_found:
                    print(f"No reads found upstream of {transcript_id}")  # Debugging statement
            except ValueError:
                print(f"Error processing transcript: {transcript_id}")  # Debugging statement
                continue
                
    print(f"Counts: {counts}")  # Debugging statement
    return counts

def normalize_counts(counts: Dict[str, int],
                    total_reads: int,
                    transcript_abundances: Dict[str, float]) -> Dict[str, float]:
    """
    Normalize read counts by total reads and transcript abundance.
    
    Args:
        counts: Raw read counts
        total_reads: Total number of mapped reads
        transcript_abundances: RPKM values for each transcript
        
    Returns:
        Dictionary of normalized counts
    """
    normalized = {}
    
    for transcript_id, count in counts.items():
        if transcript_id in transcript_abundances and transcript_abundances[transcript_id] > 0:
            normalized[transcript_id] = (count / total_reads) / transcript_abundances[transcript_id]
            
    return normalized

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
            print(f"Total reads in {os.path.basename(bam_file)}: {total_reads}") # Debugging statement
        
        # Count upstream reads
        counts = count_upstream_reads(bam_file, eligible_transcripts)
        
        # Load transcript abundances (you'll need to implement this based on your RNA-seq data)
        # For now, using placeholder values
        transcript_abundances = {tid: 1.0 for tid in eligible_transcripts}
        print(f"Transcript Abundances: {transcript_abundances}") # Debugging statement
        # Normalize counts
        normalized_counts = normalize_counts(counts, total_reads, transcript_abundances)
        print(f"Normalized Counts: {normalized_counts}")  # Debugging statement
        
        # Calculate Q3
        if normalized_counts:
            q3 = np.percentile(list(normalized_counts.values()), 75)
            print(f"Q3 for {os.path.basename(bam_file)}: {q3:.2f}")
        
        results.append({
            'sample': os.path.basename(bam_file),
            'normalized_counts': normalized_counts,
            'q3': q3
        })
    
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
                'q3': result['q3']
            })
    
    pd.DataFrame(output_data).to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
