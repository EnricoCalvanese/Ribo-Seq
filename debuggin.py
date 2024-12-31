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
            chr_id = f"Chr{chr_num}"
            return chr_id
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

def main():
    """Main analysis workflow with additional debugging"""
    print("Starting translation analysis...")
    
    # Create/load genome database
    db = gffutils.FeatureDB(DB_PATH)
    
    # Load transcript annotations
    print("\nLoading transcript structures...")
    transcript_info = {
        'exon_coords': defaultdict(list),
        'cds_coords': defaultdict(list),
        'exon_lengths': {},
        'cds_lengths': {},
        'mAUG_positions': {},
        'transcript_ids': set()
    }
    
    # Debug counters
    processed_transcripts = 0
    transcripts_with_exons = 0
    transcripts_with_cds = 0
    
    # Load mRNA features
    for transcript in db.features_of_type('mRNA'):
        processed_transcripts += 1
        transcript_id = transcript.id
        transcript_info['transcript_ids'].add(transcript_id)
        
        if processed_transcripts <= 5:  # Print details for first 5 transcripts
            print(f"\nProcessing transcript: {transcript_id}")
            print(f"Strand: {transcript.strand}")
        
        # Get exons and CDS regions
        exon_count = 0
        for exon in db.children(transcript, featuretype='exon'):
            exon_count += 1
            transcript_info['exon_coords'][transcript_id].append((exon.start, exon.end))
        
        if exon_count > 0:
            transcripts_with_exons += 1
            
        cds_count = 0
        for cds in db.children(transcript, featuretype='CDS'):
            cds_count += 1
            transcript_info['cds_coords'][transcript_id].append((cds.start, cds.end))
        
        if cds_count > 0:
            transcripts_with_cds += 1
            
        if processed_transcripts <= 5:
            print(f"Exons: {exon_count}")
            print(f"CDS regions: {cds_count}")
        
        # Calculate lengths
        if transcript_id in transcript_info['exon_coords']:
            transcript_info['exon_lengths'][transcript_id] = sum(
                end - start + 1 for start, end in transcript_info['exon_coords'][transcript_id]
            )
        if transcript_id in transcript_info['cds_coords']:
            transcript_info['cds_lengths'][transcript_id] = sum(
                end - start + 1 for start, end in transcript_info['cds_coords'][transcript_id]
            )
            # Set mAUG position
            if transcript.strand == '+':
                transcript_info['mAUG_positions'][transcript_id] = min(start for start, _ in transcript_info['cds_coords'][transcript_id])
            else:
                transcript_info['mAUG_positions'][transcript_id] = max(end for _, end in transcript_info['cds_coords'][transcript_id])
    
    print(f"\nTranscript loading summary:")
    print(f"Total transcripts processed: {processed_transcripts}")
    print(f"Transcripts with exons: {transcripts_with_exons}")
    print(f"Transcripts with CDS: {transcripts_with_cds}")
    
    # Continue with the rest of the analysis...
    # [Rest of the main function remains the same]

if __name__ == "__main__":
    main()
