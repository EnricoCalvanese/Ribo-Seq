import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
from typing import List, Dict, Set, Tuple

def examine_bam_references(bam_file: str, num_examples: int = 5) -> Set[str]:
    """
    Examines the reference names in a BAM file and prints examples.
    
    Args:
        bam_file: Path to BAM file
        num_examples: Number of example references to print
        
    Returns:
        Set of all reference names in the BAM file
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        references = set(bam.references)
        print(f"\nFound {len(references)} unique references in BAM file")
        print("Example reference names:")
        for ref in list(references)[:num_examples]:
            print(f"  {ref}")
    return references

def examine_transcript_ids(transcripts: Set[str], num_examples: int = 5) -> None:
    """
    Prints examples of transcript IDs from the GFF database.
    
    Args:
        transcripts: Set of transcript IDs
        num_examples: Number of examples to print
    """
    print("\nExample transcript IDs from GFF:")
    for transcript_id in list(transcripts)[:num_examples]:
        print(f"  {transcript_id}")

def clean_transcript_id(transcript_id: str) -> str:
    """
    Standardizes transcript ID format by removing prefixes and cleaning up the ID.
    
    Args:
        transcript_id: Raw transcript ID
        
    Returns:
        Cleaned transcript ID
    """
    # First, print the original ID for debugging
    print(f"Cleaning transcript ID: {transcript_id}")
    
    # Remove common prefixes
    prefixes_to_remove = ['transcript:', 'gene:', 'mRNA:']
    cleaned_id = transcript_id
    for prefix in prefixes_to_remove:
        if cleaned_id.startswith(prefix):
            cleaned_id = cleaned_id[len(prefix):]
    
    # Print the cleaned ID
    print(f"Cleaned to: {cleaned_id}")
    return cleaned_id

def verify_bam_references(bam_file: str, transcripts: Set[str]) -> Set[str]:
    """
    Verifies which transcripts are present in the BAM file.
    Now includes detailed debugging information.
    
    Args:
        bam_file: Path to BAM file
        transcripts: Set of transcript IDs to verify
        
    Returns:
        Set of valid transcript IDs
    """
    valid_transcripts = set()
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_references = set(bam.references)
        print(f"\nBAM file contains {len(bam_references)} references")
        
        # Print some example BAM references
        print("\nExample BAM references:")
        for ref in list(bam_references)[:5]:
            print(f"  {ref}")
        
        # Check each transcript
        print("\nChecking transcript IDs against BAM references...")
        for transcript_id in transcripts:
            clean_id = clean_transcript_id(transcript_id)
            if clean_id in bam_references:
                valid_transcripts.add(clean_id)
            elif len(valid_transcripts) == 0:  # Only print first few failures
                print(f"Failed to find: {clean_id}")
    
    return valid_transcripts

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
    
    # Load GFF database
    print("Loading GFF database...")
    db = gffutils.FeatureDB(GFF3_FILE + '.db')
    
    # Examine first BAM file to understand reference format
    print("\nExamining first BAM file format...")
    bam_references = examine_bam_references(RIBO_SEQ_FILES[0])
    
    # Process a few transcripts to understand the format differences
    print("\nProcessing sample transcripts...")
    sample_transcripts = set(list(db.features_of_type('mRNA'))[:5])
    examine_transcript_ids(sample_transcripts)
    
    # Print full processing information for debugging
    print("\nFull processing information:")
    for transcript in sample_transcripts:
        print(f"\nOriginal transcript ID: {transcript.id}")
        cleaned_id = clean_transcript_id(transcript.id)
        print(f"Cleaned ID: {cleaned_id}")
        print(f"Found in BAM references: {cleaned_id in bam_references}")

if __name__ == "__main__":
    main()
