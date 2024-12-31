import pandas as pd
import numpy as np
import pysam
import gffutils
from collections import defaultdict
import os
import sys

# Base directory structure and file paths remain the same
BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
GENOME_GFF = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")
OUTPUT_DIR = os.path.join(BASE_DIR, "systemPipeR/translating_AUGs")
DB_PATH = os.path.join(OUTPUT_DIR, "genome.db")

def get_chromosome_id(transcript_id):
    """Extract chromosome number from transcript ID and format for BAM"""
    try:
        if transcript_id.startswith('transcript:AT'):
            chr_num = transcript_id.split('G')[0][-1]
            # Just return the number instead of "Chr{chr_num}"
            return chr_num
        print(f"Unmatched transcript ID format: {transcript_id}")
    except Exception as e:
        print(f"Error processing transcript ID {transcript_id}: {str(e)}")
    return None

def calculate_rpkm(bam_file, feature_coords, feature_lengths):
    """Calculate RPKM values for genomic features with debugging output"""
    rpkm_values = {}
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        total_mapped_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        print(f"\nTotal mapped reads: {total_mapped_reads}")
        
        # Print available chromosomes in BAM
        valid_chromosomes = set(bam.references)
        print(f"Valid chromosomes in BAM: {valid_chromosomes}")
        
        # Debug counter
        processed = 0
        matched_chroms = 0
        
        for transcript_id, regions in feature_coords.items():
            processed += 1
            if processed <= 5:  # Print first 5 transcripts for debugging
                print(f"\nProcessing transcript: {transcript_id}")
            
            count = 0
            chrom = get_chromosome_id(transcript_id)
            
            if processed <= 5:
                print(f"Mapped to chromosome: {chrom}")
            
            if chrom and chrom in valid_chromosomes:
                matched_chroms += 1
                for start, end in regions:
                    try:
                        region_count = sum(1 for read in bam.fetch(chrom, start-1, end)
                                         if not read.is_unmapped)
                        count += region_count
                        if processed <= 5:
                            print(f"Region {start}-{end}: {region_count} reads")
                    except ValueError as e:
                        if processed <= 5:
                            print(f"Error fetching region {start}-{end}: {str(e)}")
                        continue
                
                if transcript_id in feature_lengths and feature_lengths[transcript_id] > 0:
                    rpkm = (count * 1e9) / (total_mapped_reads * feature_lengths[transcript_id])
                    rpkm_values[transcript_id] = rpkm
                    if processed <= 5:
                        print(f"RPKM: {rpkm}")
            
            if processed % 1000 == 0:
                print(f"Processed {processed} transcripts...")
        
        print(f"\nProcessing summary:")
        print(f"Total transcripts processed: {processed}")
        print(f"Chromosomes matched: {matched_chroms}")
        print(f"Transcripts with RPKM values: {len(rpkm_values)}")
        
        bam.close()
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
    
    return rpkm_values

# [Rest of the code remains the same as in the previous artifact]
