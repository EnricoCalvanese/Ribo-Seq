###Transcript Filtering Script Based on RPKM Thresholds
import os
import re
import pysam
import numpy as np
from typing import List, Dict, Tuple, Set
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def parse_gff3_attributes(attr_string: str) -> dict:
    """
    Parse GFF3 attribute string into a dictionary.
    Example: ID=transcript:AT1G01010.1;Parent=gene:AT1G01010 -> {'ID': 'transcript:AT1G01010.1', ...}
    """
    attrs = {}
    for attr in attr_string.strip(';').split(';'):
        if '=' not in attr:
            continue
        key, value = attr.split('=')
        attrs[key] = value
    return attrs

def load_regions_from_gff3(gff3_file: str) -> Tuple[Dict[str, List[Tuple[str, int, int]]], 
                                                    Dict[str, List[Tuple[str, int, int]]]]:
    """
    Parse TAIR10 GFF3 file to extract CDS and exon regions.
    Returns regions as (chromosome, start, end) to handle chromosome mapping.
    """
    cds_regions = defaultdict(list)
    exon_regions = defaultdict(list)
    parent_map = {}  # Map to link CDS/exon features to their transcript IDs
    
    logging.info(f"Loading regions from {gff3_file}")
    
    with open(gff3_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attributes = parse_gff3_attributes(fields[8])
            
            # Handle different feature types
            if feature_type == 'mRNA':
                # Store transcript ID mapping
                transcript_id = attributes.get('ID', '')
                parent_map[transcript_id] = transcript_id
            elif feature_type == 'CDS':
                parent = attributes.get('Parent', '')
                if parent in parent_map:
                    transcript_id = parent_map[parent]
                    cds_regions[transcript_id].append((chrom, start, end))
            elif feature_type == 'exon':
                parent = attributes.get('Parent', '')
                if parent in parent_map:
                    transcript_id = parent_map[parent]
                    exon_regions[transcript_id].append((chrom, start, end))
    
    logging.info(f"Loaded {len(cds_regions)} transcripts with CDS regions")
    logging.info(f"Loaded {len(exon_regions)} transcripts with exon regions")
    
    return dict(cds_regions), dict(exon_regions)

def calculate_rpkm(bam_file: str, 
                  regions: Dict[str, List[Tuple[str, int, int]]], 
                  region_type: str = 'CDS') -> Dict[str, float]:
    """
    Calculate RPKM values for genomic regions from a BAM file.
    Handles chromosome names from the TAIR10 annotation.
    """
    logging.info(f"Calculating {region_type} RPKM values for {bam_file}")
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_mapped_reads = bam.mapped / 1_000_000
    
    rpkm_values = {}
    
    for transcript_id, coordinates in regions.items():
        read_count = 0
        region_length = 0
        
        # Sum up reads and lengths across all regions for this transcript
        for chrom, start, end in coordinates:
            try:
                read_count += bam.count(chrom, start, end)
                region_length += end - start
            except ValueError as e:
                logging.warning(f"Error counting reads for {transcript_id} at {chrom}:{start}-{end}: {e}")
                continue
        
        # Calculate RPKM
        if region_length > 0 and total_mapped_reads > 0:
            rpkm = (read_count * 1000 * 1_000_000) / (region_length * total_mapped_reads)
        else:
            rpkm = 0
            
        rpkm_values[transcript_id] = rpkm
    
    bam.close()
    return rpkm_values

def filter_transcripts_strict(rna_seq_files: List[str], 
                            ribo_seq_files: List[str],
                            cds_regions: Dict[str, List[Tuple[str, int, int]]],
                            exon_regions: Dict[str, List[Tuple[str, int, int]]],
                            rpkm_threshold: float = 1.0) -> Set[str]:
    """
    Filter transcripts that meet RPKM thresholds in ALL RNA-seq and Ribo-seq samples.
    """
    logging.info("Starting transcript filtering")
    
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
    all_transcripts = set(cds_regions.keys()) & set(exon_regions.keys())
    
    for transcript_id in all_transcripts:
        # Check if transcript passes threshold in ALL samples
        rna_pass = all(rpkm[transcript_id] >= rpkm_threshold for rpkm in rna_rpkms)
        ribo_pass = all(rpkm[transcript_id] >= rpkm_threshold for rpkm in ribo_rpkms)
        
        if rna_pass and ribo_pass:
            passing_transcripts.add(transcript_id)
    
    logging.info(f"Found {len(passing_transcripts)} passing transcripts")
    return passing_transcripts

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
    
    # Load regions from GFF3
    cds_regions, exon_regions = load_regions_from_gff3(GFF3_FILE)
    
    # Filter transcripts
    passing_transcripts = filter_transcripts_strict(
        RNA_SEQ_FILES,
        RIBO_SEQ_FILES,
        cds_regions,
        exon_regions,
        rpkm_threshold=1.0
    )
    
    # Save results
    output_file = os.path.join(BASE_DIR, "detected_transcripts.txt")
    with open(output_file, "w") as f:
        for transcript_id in sorted(passing_transcripts):
            f.write(f"{transcript_id}\n")
    
    logging.info(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
