###Transcript Filtering Script Based on RPKM Thresholds
import os
import pysam
import numpy as np
import pandas as pd
from typing import List, Dict
from collections import defaultdict

def calculate_rpkm(bam_file: str, regions: Dict[str, List[tuple]], 
                  region_type: str = 'CDS') -> Dict[str, float]:
    """
    Calculate RPKM values for given genomic regions from a BAM file.
    
    Args:
        bam_file: Path to the BAM file
        regions: Dictionary of transcript_id -> list of (start, end) coordinates
        region_type: String indicating region type ('CDS' or 'exon')
    
    Returns:
        Dictionary of transcript_id -> RPKM value
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_mapped_reads = bam.mapped / 1_000_000
    
    rpkm_values = {}
    
    for transcript_id, coordinates in regions.items():
        read_count = 0
        region_length = 0
        
        for start, end in coordinates:
            read_count += bam.count(transcript_id, start, end)
            region_length += end - start
        
        if region_length > 0 and total_mapped_reads > 0:
            rpkm = (read_count * 1000 * 1_000_000) / (region_length * total_mapped_reads)
        else:
            rpkm = 0
            
        rpkm_values[transcript_id] = rpkm
    
    bam.close()
    return rpkm_values

def filter_transcripts_strict(rna_seq_files: List[str], 
                            ribo_seq_files: List[str],
                            cds_regions: Dict[str, List[tuple]],
                            exon_regions: Dict[str, List[tuple]],
                            rpkm_threshold: float = 1.0) -> set:
    """
    Filter transcripts that meet RPKM thresholds in ALL RNA-seq and Ribo-seq samples.
    This implements the first step of the protocol to identify transcripts with
    detectable translation initiation from mAUGs.
    
    Args:
        rna_seq_files: List of RNA-seq BAM files
        ribo_seq_files: List of Ribo-seq BAM files
        cds_regions: Dictionary of transcript_id -> list of CDS coordinates
        exon_regions: Dictionary of transcript_id -> list of exon coordinates
        rpkm_threshold: RPKM threshold for filtering (default 1.0)
    
    Returns:
        Set of transcript IDs that pass the filtering criteria in ALL samples
    """
    # Calculate RPKM for all samples
    rna_rpkms = []
    ribo_rpkms = []
    
    # Process RNA-seq samples
    for rna_file in rna_seq_files:
        rpkm_values = calculate_rpkm(rna_file, exon_regions, 'exon')
        rna_rpkms.append(rpkm_values)
    
    # Process Ribo-seq samples
    for ribo_file in ribo_seq_files:
        rpkm_values = calculate_rpkm(ribo_file, cds_regions, 'CDS')
        ribo_rpkms.append(rpkm_values)
    
    # Find transcripts that pass thresholds in ALL samples
    passing_transcripts = set()
    
    for transcript_id in cds_regions.keys():
        # Check if transcript passes threshold in ALL samples
        rna_pass = all(rpkm[transcript_id] >= rpkm_threshold for rpkm in rna_rpkms)
        ribo_pass = all(rpkm[transcript_id] >= rpkm_threshold for rpkm in ribo_rpkms)
        
        if rna_pass and ribo_pass:
            passing_transcripts.add(transcript_id)
    
    return passing_transcripts

# Your file paths
BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
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

def main():
    """
    Main function implementing the first step of the protocol:
    1. Filter transcripts with RPKM ≥ 1 in ALL samples
    
    Subsequent steps (not implemented here) will:
    2. Calculate normalized ribosome footprints at mAUGs
    3. Calculate background from Q3 of normalized counts upstream of mAUGs
    4. Apply both normalized (≥23.17) and raw (≥10) count thresholds
    """
    # Load your CDS and exon regions
    cds_regions = {}  # Load from your annotation file
    exon_regions = {}  # Load from your annotation file
    
    # Filter transcripts with strict criteria
    passing_transcripts = filter_transcripts_strict(
        RNA_SEQ_FILES,
        RIBO_SEQ_FILES,
        cds_regions,
        exon_regions,
        rpkm_threshold=1.0
    )
    
    # Save results for next steps
    with open("detected_transcripts.txt", "w") as f:
        for transcript_id in sorted(passing_transcripts):
            f.write(f"{transcript_id}\n")
    
    print(f"Number of detected transcripts: {len(passing_transcripts)}")

if __name__ == "__main__":
    main()
