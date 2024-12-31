import pandas as pd
import numpy as np
import pysam
import gffutils
from collections import defaultdict
import os
import sys

# File paths remain the same
BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")

def parse_attributes(attr_string):
    """Parse GFF attribute string into a dictionary"""
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
    return attrs

def parse_transcript_id(feature_by):
    """Extract transcript ID from feature_by attribute"""
    if not feature_by:
        return None
    parts = feature_by.split(':')
    if len(parts) >= 2:
        return parts[1]
    return None

def parse_uorf_data():
    """Parse uORF GFF file and create transcript annotation database"""
    print("\nParsing uORF GFF file...")
    
    uorf_data = defaultdict(list)
    transcript_info = {
        'exon_coords': defaultdict(list),
        'cds_coords': defaultdict(list),
        'exon_lengths': {},
        'cds_lengths': {},
        'mAUG_positions': {},
        'transcript_ids': set(),
        'five_prime_leaders': {}
    }
    
    feature_counts = defaultdict(int)
    
    try:
        with open(UORF_GFF, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) != 9:
                    continue
                
                chrom, source, feature_type, start, end, score, strand, phase, attributes = fields
                
                # Parse attributes
                attrs = parse_attributes(attributes)
                
                if feature_type == 'sequence_feature':
                    if 'feature_by' in attrs and 'featuretype' in attrs:
                        if attrs['featuretype'] == 'uORF':
                            transcript_id = parse_transcript_id(attrs['feature_by'])
                            if transcript_id:
                                uorf_data[transcript_id].append({
                                    'start': int(start),
                                    'end': int(end),
                                    'strand': strand,
                                    'score': float(score) if score != '.' else 0,
                                    'coords': f"{start}-{end}"
                                })
                                transcript_info['transcript_ids'].add(transcript_id)
                
                feature_counts[feature_type] += 1
    
    except Exception as e:
        print(f"Error processing GFF file: {str(e)}")
        sys.exit(1)
    
    print("\nFeature counts by type:")
    for ftype, count in feature_counts.items():
        print(f"{ftype}: {count}")
    
    print(f"\nFound {len(uorf_data)} transcripts with uORFs")
    print("Sample of transcripts with uORFs:")
    sample_transcripts = list(uorf_data.keys())[:5]
    for transcript in sample_transcripts:
        print(f"\nTranscript: {transcript}")
        print(f"Number of uORFs: {len(uorf_data[transcript])}")
        print("uORF coordinates:")
        for uorf in uorf_data[transcript]:
            print(f"  {uorf['coords']} ({uorf['strand']})")
    
    return uorf_data, transcript_info

def load_transcript_annotations(gtf_file):
    """Load transcript annotations from GTF file"""
    # This function would load the basic transcript structure (exons, CDS)
    # We'll implement this if needed
    pass

def calculate_rpkm(bam_file, feature_coordinates, feature_lengths):
    """Calculate RPKM for given genomic features"""
    try:
        # Open BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
        
        # Get total mapped reads
        total_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        
        # Initialize counts
        counts = defaultdict(int)
        
        # Count reads per feature
        for feature_id, regions in feature_coordinates.items():
            for start, end in regions:
                try:
                    # Get chromosome/contig name
                    chrom = feature_id.split('.')[0]
                    
                    # Count reads in feature region
                    for read in bam.fetch(chrom, start, end):
                        if not read.is_unmapped:
                            counts[feature_id] += 1
                except ValueError:
                    continue
        
        # Calculate RPKM
        rpkm = {}
        for feature_id in counts:
            if feature_id in feature_lengths and feature_lengths[feature_id] > 0:
                rpkm[feature_id] = (counts[feature_id] * 1e9) / (total_reads * feature_lengths[feature_id])
        
        bam.close()
        return rpkm
    
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
        return {}

def main():
    """Main analysis workflow"""
    print("Starting analysis...")
    
    # Parse uORF data
    uorf_data, transcript_info = parse_uorf_data()
    
    # Here we would need to load the basic transcript structure (exons, CDS)
    # from a GTF/GFF file to proceed with the RPKM calculations
    
    print("\nNOTE: To proceed with the full analysis, we need:")
    print("1. Path to the Arabidopsis transcript annotation file (GTF/GFF)")
    print("2. Confirmation of the chromosome naming convention in your BAM files")

if __name__ == "__main__":
    main()
