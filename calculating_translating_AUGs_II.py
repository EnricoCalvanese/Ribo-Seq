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
OUTPUT_DIR = os.path.join(BASE_DIR, "translation_analysis")
os.makedirs(OUTPUT_DIR, exist_ok=True)
DB_PATH = os.path.join(OUTPUT_DIR, "genome.db")

def create_genome_db():
    """Create SQLite database from genome annotation"""
    print("\nProcessing genome annotation...")
    
    if os.path.exists(DB_PATH):
        print("Using existing genome database")
        return gffutils.FeatureDB(DB_PATH)
    
    print("Creating genome database (this may take a few minutes)...")
    db = gffutils.create_db(
        GENOME_GFF,
        DB_PATH,
        merge_strategy='create_unique',
        sort_attribute_values=True,
        disable_infer_genes=True
    )
    return db

def load_transcript_annotations(db):
    """Load transcript structures and identify 5' leaders"""
    print("Loading transcript structures...")
    
    transcript_info = {
        'exon_coords': defaultdict(list),
        'cds_coords': defaultdict(list),
        'exon_lengths': {},
        'cds_lengths': {},
        'mAUG_positions': {},
        'five_prime_leaders': {},
        'transcript_ids': set()
    }
    
    for transcript in db.features_of_type('transcript'):
        transcript_id = transcript.id
        transcript_info['transcript_ids'].add(transcript_id)
        
        # Get exons and CDS regions
        for exon in db.children(transcript, featuretype='exon'):
            transcript_info['exon_coords'][transcript_id].append((exon.start, exon.end))
        
        cds_regions = list(db.children(transcript, featuretype='CDS'))
        if cds_regions:
            for cds in cds_regions:
                transcript_info['cds_coords'][transcript_id].append((cds.start, cds.end))
            
            # Identify mAUG position
            if transcript.strand == '+':
                transcript_info['mAUG_positions'][transcript_id] = min(start for start, _ in transcript_info['cds_coords'][transcript_id])
            else:
                transcript_info['mAUG_positions'][transcript_id] = max(end for _, end in transcript_info['cds_coords'][transcript_id])
    
    # Calculate lengths
    for transcript_id in transcript_info['transcript_ids']:
        if transcript_id in transcript_info['exon_coords']:
            transcript_info['exon_lengths'][transcript_id] = sum(
                end - start + 1 for start, end in transcript_info['exon_coords'][transcript_id]
            )
        if transcript_id in transcript_info['cds_coords']:
            transcript_info['cds_lengths'][transcript_id] = sum(
                end - start + 1 for start, end in transcript_info['cds_coords'][transcript_id]
            )
    
    return transcript_info

def calculate_rpkm(bam_file, feature_coords, feature_lengths):
    """Calculate RPKM values for genomic features"""
    rpkm_values = {}
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        total_mapped_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        
        for transcript_id, regions in feature_coords.items():
            count = 0
            for start, end in regions:
                try:
                    chrom = f"Chr{transcript_id[2]}" if transcript_id[2].isdigit() else None
                    if chrom:
                        count += sum(1 for read in bam.fetch(chrom, start-1, end)
                                   if not read.is_unmapped)
                except ValueError:
                    continue
            
            if transcript_id in feature_lengths and feature_lengths[transcript_id] > 0:
                rpkm = (count * 1e9) / (total_mapped_reads * feature_lengths[transcript_id])
                rpkm_values[transcript_id] = rpkm
        
        bam.close()
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
    
    return rpkm_values

def get_expressed_transcripts(transcript_info, min_rpkm=1.0):
    """Identify transcripts with RPKM ≥ 1 in all samples"""
    expressed = set(transcript_info['transcript_ids'])
    
    # Check RNA-seq samples
    for rna_file in RNA_SEQ_FILES:
        print(f"Processing RNA-seq file: {os.path.basename(rna_file)}")
        rpkm = calculate_rpkm(rna_file, transcript_info['exon_coords'],
                            transcript_info['exon_lengths'])
        expressed &= {tid for tid, value in rpkm.items() if value >= min_rpkm}
    
    # Check Ribo-seq samples
    for ribo_file in RIBO_SEQ_FILES:
        print(f"Processing Ribo-seq file: {os.path.basename(ribo_file)}")
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        expressed &= {tid for tid, value in rpkm.items() if value >= min_rpkm}
    
    return expressed

def calculate_normalized_counts(bam_file, positions, transcript_rpkm):
    """Calculate normalized read counts at specific positions"""
    normalized_counts = {}
    raw_counts = {}
    
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        total_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        
        for transcript_id, position in positions.items():
            try:
                chrom = f"Chr{transcript_id[2]}" if transcript_id[2].isdigit() else None
                if chrom:
                    count = sum(1 for read in bam.fetch(chrom, position-1, position+1)
                              if not read.is_unmapped)
                    raw_counts[transcript_id] = count
                    if transcript_id in transcript_rpkm and transcript_rpkm[transcript_id] > 0:
                        normalized_counts[transcript_id] = (count / total_reads) / transcript_rpkm[transcript_id]
            except ValueError:
                continue
        
        bam.close()
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
    
    return normalized_counts, raw_counts

def main():
    """Main analysis workflow"""
    print("Starting translation analysis...")
    
    # Create/load genome database
    db = create_genome_db()
    
    # Load transcript annotations
    transcript_info = load_transcript_annotations(db)
    print(f"Loaded {len(transcript_info['transcript_ids'])} transcripts")
    
    # Get expressed transcripts (RPKM ≥ 1)
    expressed_transcripts = get_expressed_transcripts(transcript_info)
    print(f"Found {len(expressed_transcripts)} expressed transcripts")
    
    # Calculate background from regions 50nt upstream of mAUGs
    background_counts = []
    for ribo_file in RIBO_SEQ_FILES:
        print(f"Calculating background from {os.path.basename(ribo_file)}")
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        
        # Get positions 50nt upstream of mAUGs
        upstream_positions = {
            tid: pos - 50 for tid, pos in transcript_info['mAUG_positions'].items()
            if tid in expressed_transcripts
        }
        
        counts, _ = calculate_normalized_counts(ribo_file, upstream_positions, rpkm)
        background_counts.extend(counts.values())
    
    # Calculate Q3 cutoff
    background_cutoff = np.percentile(background_counts, 75)
    print(f"Background cutoff (Q3): {background_cutoff:.2f}")
    
    # Save results
    with open(os.path.join(OUTPUT_DIR, "analysis_results.txt"), "w") as f:
        f.write(f"Total transcripts analyzed: {len(transcript_info['transcript_ids'])}\n")
        f.write(f"Expressed transcripts (RPKM ≥ 1): {len(expressed_transcripts)}\n")
        f.write(f"Background cutoff: {background_cutoff:.2f}\n")
        f.write("\nExpressed transcript IDs:\n")
        f.write("\n".join(sorted(expressed_transcripts)))

if __name__ == "__main__":
    main()
