import os
import pysam
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from collections import defaultdict
import logging
from dataclasses import dataclass

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

@dataclass
class StartCodonPosition:
    """Store information about a start codon location and its read counts."""
    chrom: str
    position: int  # Position of the first base of start codon
    counts: Dict[str, int]  # Sample -> raw count mapping
    normalized_counts: Dict[str, float]  # Sample -> normalized count mapping

def calculate_p_site_offset(read_length: int) -> int:
    """
    Calculate P-site offset based on read length.
    For most organisms:
    - 28-nt reads: 12-nt offset
    - 29-nt reads: 12-nt offset
    - 30-nt reads: 13-nt offset
    """
    if read_length <= 29:
        return 12
    else:
        return 13

def get_start_codon_positions(gff3_file: str) -> Dict[str, StartCodonPosition]:
    """
    Extract start codon positions from GFF3 file.
    """
    start_codons = {}
    
    with open(gff3_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
                
            chrom = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            attrs = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
            
            # Get parent transcript ID
            parent = attrs.get('Parent', '')
            if not parent:
                continue
                
            # For positive strand, start codon is at start position
            # For negative strand, start codon is at end position - 2
            start_pos = start if strand == '+' else end - 2
            start_codons[parent] = StartCodonPosition(
                chrom=chrom,
                position=start_pos,
                counts={},
                normalized_counts={}
            )
    
    logging.info(f"Found {len(start_codons)} start codon positions")
    return start_codons

def count_p_site_reads(bam_file: str, 
                      start_codons: Dict[str, StartCodonPosition],
                      window_size: int = 1) -> None:
    """
    Count P-site mapped reads around start codons.
    Updates the counts in the StartCodonPosition objects.
    """
    sample_name = os.path.basename(bam_file)
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_mapped_reads = bam.mapped
    
    for transcript_id, start_codon in start_codons.items():
        read_count = 0
        
        # Define region to check (window around start codon)
        region_start = start_codon.position - window_size
        region_end = start_codon.position + window_size + 3  # Include full start codon
        
        # Get all reads in region
        try:
            for read in bam.fetch(start_codon.chrom, region_start, region_end):
                # Calculate P-site position for this read
                p_site_offset = calculate_p_site_offset(read.query_length)
                p_site_pos = read.reference_start + p_site_offset
                
                # Check if P-site falls within our window around start codon
                if region_start <= p_site_pos <= region_end:
                    read_count += 1
        except ValueError:
            logging.warning(f"Could not fetch reads for {transcript_id} at {start_codon.chrom}:{region_start}-{region_end}")
            
        # Store raw count
        start_codon.counts[sample_name] = read_count
        
        # Calculate normalized count (per million mapped reads)
        normalized_count = (read_count * 1_000_000) / total_mapped_reads
        start_codon.normalized_counts[sample_name] = normalized_count
    
    bam.close()

def filter_by_start_codon_usage(start_codons: Dict[str, StartCodonPosition],
                              ribo_seq_files: List[str],
                              raw_count_threshold: int = 10,
                              norm_count_threshold: float = 23.17) -> set:
    """
    Filter transcripts based on start codon usage criteria.
    """
    logging.info("Starting start codon usage analysis")
    
    # Calculate P-site reads for all samples
    for bam_file in ribo_seq_files:
        logging.info(f"Processing {bam_file}")
        count_p_site_reads(bam_file, start_codons)
    
    # Filter based on criteria
    passing_transcripts = set()
    
    for transcript_id, start_codon in start_codons.items():
        # Check raw count threshold in all samples
        raw_counts_pass = all(count >= raw_count_threshold 
                            for count in start_codon.counts.values())
        
        # Check normalized count threshold in all samples
        norm_counts_pass = all(count >= norm_count_threshold 
                             for count in start_codon.normalized_counts.values())
        
        if raw_counts_pass and norm_counts_pass:
            passing_transcripts.add(transcript_id)
    
    logging.info(f"Found {len(passing_transcripts)} transcripts with significant start codon usage")
    return passing_transcripts

def main():
    # File paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    
    RIBO_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT103-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT103-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-2_uniq_sort.bam")
    ]
    
    # Load previously detected transcripts
    detected_transcripts = set()
    with open("/global/scratch/users/enricocalvane/riboseq/imb2/systemPipeR/translating_AUGs/detected_transcripts.txt") as f:
        for line in f:
            detected_transcripts.add(line.strip())
    
    logging.info(f"Loaded {len(detected_transcripts)} previously detected transcripts")
    
    # Get start codon positions
    start_codons = get_start_codon_positions(GFF3_FILE)
    
    # Filter to previously detected transcripts
    start_codons = {tid: pos for tid, pos in start_codons.items() 
                   if tid in detected_transcripts}
    
    # Filter based on start codon usage
    passing_transcripts = filter_by_start_codon_usage(
        start_codons,
        RIBO_SEQ_FILES,
        raw_count_threshold=10,
        norm_count_threshold=23.17
    )
    
    # Save results
    output_file = os.path.join(BASE_DIR, "expressed_transcripts.txt")
    with open(output_file, "w") as f:
        for transcript_id in sorted(passing_transcripts):
            f.write(f"{transcript_id}\n")
    
    # Save detailed metrics
    metrics_file = os.path.join(BASE_DIR, "start_codon_metrics.tsv")
    with open(metrics_file, "w") as f:
        # Write header
        samples = sorted(next(iter(start_codons.values())).counts.keys())
        header = ["transcript_id"] + \
                [f"{s}_raw_count" for s in samples] + \
                [f"{s}_normalized_count" for s in samples]
        f.write("\t".join(header) + "\n")
        
        # Write metrics for each transcript
        for tid, pos in start_codons.items():
            raw_counts = [str(pos.counts[s]) for s in samples]
            norm_counts = [f"{pos.normalized_counts[s]:.2f}" for s in samples]
            row = [tid] + raw_counts + norm_counts
            f.write("\t".join(row) + "\n")
    
    logging.info(f"Results saved to {output_file}")
    logging.info(f"Detailed metrics saved to {metrics_file}")

if __name__ == "__main__":
    main()
