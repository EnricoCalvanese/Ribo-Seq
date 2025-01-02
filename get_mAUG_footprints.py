import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
from typing import List, Dict, Tuple

def create_gff_db(gff3_file: str) -> gffutils.FeatureDB:
    """
    Create or connect to a GFF database for efficient querying.
    
    Args:
        gff3_file: Path to the GFF3 file
        
    Returns:
        gffutils database object
    """
    db_path = gff3_file + '.db'
    if not os.path.exists(db_path):
        gffutils.create_db(gff3_file, db_path)
    return gffutils.FeatureDB(db_path)

def load_detected_transcripts(filepath: str) -> set:
    """
    Load the filtered transcript IDs from file.
    
    Args:
        filepath: Path to the detected transcripts file
        
    Returns:
        Set of transcript IDs
    """
    with open(filepath, 'r') as f:
        return set(line.strip() for line in f)

def get_total_read_counts(bam_files: List[str]) -> Dict[str, int]:
    """
    Calculate total mapped reads for each BAM file.
    
    Args:
        bam_files: List of paths to BAM files
        
    Returns:
        Dictionary mapping file paths to total read counts
    """
    counts = {}
    for bam_file in bam_files:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            counts[bam_file] = bam.count()
    return counts

def find_mAUG_positions(db: gffutils.FeatureDB, transcript_ids: set) -> Dict[str, List[int]]:
    """
    Find all mAUG positions in the filtered transcripts.
    
    Args:
        db: GFF database
        transcript_ids: Set of filtered transcript IDs
        
    Returns:
        Dictionary mapping transcript IDs to lists of mAUG positions
    """
    mAUG_positions = defaultdict(list)
    
    for transcript_id in transcript_ids:
        transcript = db[transcript_id]
        seq = transcript.sequence  # Assuming sequence is available in GFF
        
        # Find all AUG positions
        for i in range(len(seq)-2):
            if seq[i:i+3].upper() == 'AUG':
                mAUG_positions[transcript_id].append(i)
    
    return mAUG_positions

def count_footprints(bam_file: str, 
                    mAUG_positions: Dict[str, List[int]], 
                    window_size: int = 3) -> Dict[str, Dict[int, int]]:
    """
    Count ribosome footprints around mAUG sites.
    
    Args:
        bam_file: Path to BAM file
        mAUG_positions: Dictionary of mAUG positions per transcript
        window_size: Size of the window around mAUG to count reads
        
    Returns:
        Nested dictionary of counts per transcript and mAUG position
    """
    counts = defaultdict(lambda: defaultdict(int))
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for transcript_id, positions in mAUG_positions.items():
            for pos in positions:
                # Count reads in window around mAUG
                for read in bam.fetch(transcript_id, 
                                    pos - window_size,
                                    pos + window_size):
                    counts[transcript_id][pos] += 1
                    
    return counts

def calculate_transcript_abundance(rna_seq_files: List[str],
                                transcript_ids: set) -> Dict[str, float]:
    """
    Calculate transcript abundance (RPKM) from RNA-seq data.
    
    Args:
        rna_seq_files: List of RNA-seq BAM files
        transcript_ids: Set of filtered transcript IDs
        
    Returns:
        Dictionary of transcript abundances
    """
    abundances = defaultdict(float)
    total_counts = get_total_read_counts(rna_seq_files)
    
    for bam_file in rna_seq_files:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for transcript_id in transcript_ids:
                count = sum(1 for _ in bam.fetch(transcript_id))
                length = bam.get_reference_length(transcript_id)
                rpkm = (count * 1e9) / (total_counts[bam_file] * length)
                abundances[transcript_id] += rpkm / len(rna_seq_files)
    
    return abundances

def main():
    # File paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    
    RNA_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT101-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT101-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT102-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT102-2_uniq_sort.bam")
    ]
    
    RIBO_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT103-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT103-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-2_uniq_sort.bam")
    ]
    
    # Load detected transcripts
    detected_transcripts = load_detected_transcripts(
        os.path.join(BASE_DIR, "systemPipeR/translating_AUGs/detected_transcripts.txt")
    )
    
    # Create GFF database
    db = create_gff_db(GFF3_FILE)
    
    # Find mAUG positions
    mAUG_positions = find_mAUG_positions(db, detected_transcripts)
    
    # Calculate transcript abundances
    abundances = calculate_transcript_abundance(RNA_SEQ_FILES, detected_transcripts)
    
    # Process each Ribo-seq file
    results = []
    total_counts = get_total_read_counts(RIBO_SEQ_FILES)
    
    for bam_file in RIBO_SEQ_FILES:
        # Count footprints at mAUG sites
        footprint_counts = count_footprints(bam_file, mAUG_positions)
        
        # Normalize counts
        for transcript_id, positions in footprint_counts.items():
            for pos, count in positions.items():
                # Normalize by total reads and transcript abundance
                normalized_count = (count / total_counts[bam_file]) / abundances[transcript_id]
                
                results.append({
                    'transcript_id': transcript_id,
                    'mAUG_position': pos,
                    'sample': os.path.basename(bam_file),
                    'raw_count': count,
                    'normalized_count': normalized_count
                })
    
    # Save results
    results_df = pd.DataFrame(results)
    output_file = os.path.join(BASE_DIR, "systemPipeR/translating_AUGs/mAUG_footprints.tsv")
    results_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
