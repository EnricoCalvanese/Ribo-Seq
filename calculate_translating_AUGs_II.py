import pandas as pd
import numpy as np
import pysam
import gffutils
from collections import defaultdict
import os
import sys

# Base directory structure
BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"

# Input files
GENOME_GFF = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")

# BAM file paths
RNA_SEQ_FILES = [
    os.path.join(BASE_DIR, "unique_reads", f"LZT101-{i}_uniq_sort.bam") for i in [1, 2]
] + [
    os.path.join(BASE_DIR, "unique_reads", f"LZT102-{i}_uniq_sort.bam") for i in [1, 2]
]

RIBO_SEQ_FILES = [
    os.path.join(BASE_DIR, "unique_reads", f"LZT103-{i}_uniq_sort.bam") for i in [1, 2]
] + [
    os.path.join(BASE_DIR, "unique_reads", f"LZT104-{i}_uniq_sort.bam") for i in [1, 2]
]

# Output directory and files
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

def get_expressed_transcripts(transcript_info, min_rpkm=1.0):
    """Identify transcripts with RPKM ≥ 1 in all samples"""
    expressed = set(transcript_info['transcript_ids'])
    print(f"Starting with {len(expressed)} transcripts")
    
    # Check RNA-seq samples
    for rna_file in RNA_SEQ_FILES:
        print(f"Processing RNA-seq file: {os.path.basename(rna_file)}")
        rpkm = calculate_rpkm(rna_file, transcript_info['exon_coords'],
                            transcript_info['exon_lengths'])
        expressed_in_sample = {tid for tid, value in rpkm.items() if value >= min_rpkm}
        expressed &= expressed_in_sample
        print(f"Transcripts with RPKM ≥ {min_rpkm}: {len(expressed_in_sample)}")
        print(f"Remaining expressed transcripts: {len(expressed)}")
    
    # Check Ribo-seq samples
    for ribo_file in RIBO_SEQ_FILES:
        print(f"Processing Ribo-seq file: {os.path.basename(ribo_file)}")
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        expressed_in_sample = {tid for tid, value in rpkm.items() if value >= min_rpkm}
        expressed &= expressed_in_sample
        print(f"Transcripts with RPKM ≥ {min_rpkm}: {len(expressed_in_sample)}")
        print(f"Remaining expressed transcripts: {len(expressed)}")
    
    return expressed

def calculate_normalized_counts(bam_file, positions, transcript_rpkm, min_raw_count=1):
    """
    Calculate normalized read counts with improved counting strategy
    Returns both normalized and raw counts that meet minimum threshold
    """
    normalized_counts = {}
    raw_counts = {}
    
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        total_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        valid_chromosomes = set(bam.references)
        
        for transcript_id, position in positions.items():
            chrom = get_chromosome_id(transcript_id)
            if chrom and chrom in valid_chromosomes:
                try:
                    # Extend the window to capture reads spanning the position
                    window_start = max(0, position - 15)  # Adjust window size as needed
                    window_end = position + 15
                    
                    count = sum(1 for read in bam.fetch(chrom, window_start, window_end)
                              if not read.is_unmapped and
                              min(read.reference_end, window_end) - max(read.reference_start, window_start) > 0)
                    
                    if count >= min_raw_count:
                        raw_counts[transcript_id] = count
                        if transcript_id in transcript_rpkm and transcript_rpkm[transcript_id] > 0:
                            # Normalize by total reads and transcript expression
                            norm_count = (count * 1e6) / (total_reads * transcript_rpkm[transcript_id])
                            normalized_counts[transcript_id] = norm_count
                
                except ValueError as e:
                    continue
        
        bam.close()
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
    
    return normalized_counts, raw_counts

def calculate_background_cutoff(expressed_transcripts, transcript_info, ribo_files):
    """
    Calculate background cutoff from regions upstream of mAUGs
    """
    print("\nCalculating background cutoff...")
    background_counts = []
    raw_counts_dist = []
    
    for ribo_file in ribo_files:
        print(f"Processing {os.path.basename(ribo_file)}")
        
        # Calculate RPKM for normalization
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        
        # Get positions 50nt upstream of mAUGs
        upstream_positions = {}
        for tid in expressed_transcripts:
            if tid in transcript_info['mAUG_positions']:
                mAUG_pos = transcript_info['mAUG_positions'][tid]
                # Ensure position is at least 50nt upstream
                upstream_positions[tid] = max(1, mAUG_pos - 50)
        
        # Calculate counts at upstream positions
        norm_counts, raw_counts = calculate_normalized_counts(
            ribo_file, upstream_positions, rpkm, min_raw_count=1
        )
        
        background_counts.extend(norm_counts.values())
        raw_counts_dist.extend(raw_counts.values())
        
        print(f"Found {len(norm_counts)} normalized counts in this sample")
        if len(norm_counts) > 0:
            print(f"Range: {min(norm_counts.values()):.2f} - {max(norm_counts.values()):.2f}")
    
    if not background_counts:
        print("WARNING: No background counts found!")
        return 0.0
    
    # Calculate and return the 75th percentile (Q3)
    background_cutoff = np.percentile(background_counts, 75)
    
    print("\nBackground calculation summary:")
    print(f"Total background measurements: {len(background_counts)}")
    print(f"Raw counts statistics:")
    print(f"  Min: {min(raw_counts_dist)}")
    print(f"  Max: {max(raw_counts_dist)}")
    print(f"  Mean: {np.mean(raw_counts_dist):.2f}")
    print(f"Normalized counts statistics:")
    print(f"  Min: {min(background_counts):.2f}")
    print(f"  Max: {max(background_counts):.2f}")
    print(f"  Mean: {np.mean(background_counts):.2f}")
    print(f"  Q3 (75th percentile): {background_cutoff:.2f}")
    
    return background_cutoff

def main():
    """Main analysis workflow with improved background calculation"""
    print("Starting translation analysis...")
    
    # Load database and transcript info (previous code remains the same)
    db = gffutils.FeatureDB(DB_PATH)
    transcript_info = load_transcript_annotations(db)
    
    # Get expressed transcripts
    expressed_transcripts = get_expressed_transcripts(transcript_info)
    print(f"\nFound {len(expressed_transcripts)} expressed transcripts")
    
    if len(expressed_transcripts) == 0:
        print("No expressed transcripts found. Check RPKM calculations and BAM files.")
        return
    
    # Calculate background cutoff with new function
    background_cutoff = calculate_background_cutoff(
        expressed_transcripts, transcript_info, RIBO_SEQ_FILES
    )
    
    # Save results with more detailed information
    output_file = os.path.join(OUTPUT_DIR, "analysis_results.txt")
    print(f"\nSaving results to {output_file}")
    
    with open(output_file, "w") as f:
        f.write("Analysis Results\n")
        f.write("===============\n\n")
        f.write(f"Total transcripts analyzed: {len(transcript_info['transcript_ids'])}\n")
        f.write(f"Expressed transcripts (RPKM ≥ 1): {len(expressed_transcripts)}\n")
        f.write(f"Background cutoff: {background_cutoff:.2f}\n\n")
        f.write("Details:\n")
        f.write("--------\n")
        f.write("Criteria used:\n")
        f.write("- RPKM ≥ 1 in all RNA-seq and Ribo-seq samples\n")
        f.write("- Background calculated from 50nt upstream of mAUGs\n")
        f.write("- Q3 (75th percentile) used as cutoff\n\n")
        f.write("Expressed transcript IDs:\n")
        f.write("----------------------\n")
        f.write("\n".join(sorted(expressed_transcripts)))

if __name__ == "__main__":
    main()
