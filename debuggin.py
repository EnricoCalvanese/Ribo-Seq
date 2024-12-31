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

def create_genome_db():
    """Create SQLite database from genome annotation with debug info"""
    print("\nProcessing genome annotation...")
    
    if os.path.exists(DB_PATH):
        print("Using existing genome database")
        db = gffutils.FeatureDB(DB_PATH)
        
        # Debug: Check database contents
        feature_counts = defaultdict(int)
        print("\nChecking database contents...")
        for feature in db.all_features():
            feature_counts[feature.featuretype] += 1
        
        print("\nFeature counts in database:")
        for ftype, count in feature_counts.items():
            print(f"{ftype}: {count}")
        
        return db
    
    print("Creating new genome database...")
    db = gffutils.create_db(
        GENOME_GFF,
        DB_PATH,
        merge_strategy='create_unique',
        sort_attribute_values=True,
        disable_infer_genes=True,
        verbose=True  # Add verbose output
    )
    return db

def load_transcript_annotations(db):
    """Load transcript structures with additional debug information"""
    print("\nLoading transcript structures...")
    
    transcript_info = {
        'exon_coords': defaultdict(list),
        'cds_coords': defaultdict(list),
        'exon_lengths': {},
        'cds_lengths': {},
        'mAUG_positions': {},
        'five_prime_leaders': {},
        'transcript_ids': set()
    }
    
    # Debug: Print available feature types
    print("\nAvailable feature types in database:")
    for ftype in db.featuretypes():
        print(f"- {ftype}")
    
    # Try different transcript identifiers
    transcript_types = ['transcript', 'mRNA']
    transcripts_found = False
    
    for transcript_type in transcript_types:
        print(f"\nTrying to load features of type: {transcript_type}")
        count = 0
        
        try:
            for transcript in db.features_of_type(transcript_type):
                count += 1
                transcript_id = transcript.id
                transcript_info['transcript_ids'].add(transcript_id)
                
                # Get exons
                exon_count = 0
                for exon in db.children(transcript, featuretype='exon'):
                    exon_count += 1
                    transcript_info['exon_coords'][transcript_id].append((exon.start, exon.end))
                
                # Get CDS regions
                cds_count = 0
                cds_regions = list(db.children(transcript, featuretype='CDS'))
                for cds in cds_regions:
                    cds_count += 1
                    transcript_info['cds_coords'][transcript_id].append((cds.start, cds.end))
                
                if cds_regions:
                    # Identify mAUG position
                    if transcript.strand == '+':
                        transcript_info['mAUG_positions'][transcript_id] = min(start for start, _ in transcript_info['cds_coords'][transcript_id])
                    else:
                        transcript_info['mAUG_positions'][transcript_id] = max(end for _, end in transcript_info['cds_coords'][transcript_id])
                
                if count % 1000 == 0:
                    print(f"Processed {count} {transcript_type}s...")
                
                # Debug output for first few transcripts
                if count <= 5:
                    print(f"\nExample transcript {count}:")
                    print(f"ID: {transcript_id}")
                    print(f"Exons: {exon_count}")
                    print(f"CDS regions: {cds_count}")
        
            print(f"\nTotal {transcript_type}s processed: {count}")
            if count > 0:
                transcripts_found = True
                break
                
        except Exception as e:
            print(f"Error processing {transcript_type}: {str(e)}")
    
    if not transcripts_found:
        print("\nWARNING: No transcripts were loaded! Database might be empty or incorrectly formatted.")
        print("Checking database file...")
        print(f"Database path: {DB_PATH}")
        print(f"Database exists: {os.path.exists(DB_PATH)}")
        print(f"Database size: {os.path.getsize(DB_PATH) if os.path.exists(DB_PATH) else 'N/A'} bytes")
    
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

def main():
    """Main analysis workflow with debug output"""
    print("Starting translation analysis...")
    
    # Create/load genome database
    db = create_genome_db()
    
    # Load transcript annotations
    transcript_info = load_transcript_annotations(db)
    
    print("\nSummary statistics:")
    print(f"Total transcripts: {len(transcript_info['transcript_ids'])}")
    print(f"Transcripts with exons: {len(transcript_info['exon_coords'])}")
    print(f"Transcripts with CDS: {len(transcript_info['cds_coords'])}")
    print(f"Transcripts with mAUG positions: {len(transcript_info['mAUG_positions'])}")
    
    if len(transcript_info['transcript_ids']) == 0:
        print("\nNo transcripts loaded. Stopping analysis.")
        return
        
    # Rest of the analysis code...

if __name__ == "__main__":
    main()
