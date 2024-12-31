import pandas as pd
import numpy as np
import pysam
import gffutils
from collections import defaultdict
import os
import sys

# File paths
BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")
GENOME_GFF = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"

def create_genome_db():
    """Create SQLite database from genome annotation"""
    print("\nProcessing genome annotation...")
    db_path = "genome.db"
    
    if os.path.exists(db_path):
        print("Using existing genome database")
        return gffutils.FeatureDB(db_path)
    
    print("Creating genome database (this may take a few minutes)...")
    db = gffutils.create_db(
        GENOME_GFF,
        db_path,
        merge_strategy='create_unique',
        sort_attribute_values=True,
        disable_infer_genes=True
    )
    return db

def load_transcript_annotations(db):
    """Load transcript structures from genome annotation"""
    print("Loading transcript structures...")
    
    transcript_info = {
        'exon_coords': defaultdict(list),
        'cds_coords': defaultdict(list),
        'exon_lengths': {},
        'cds_lengths': {},
        'mAUG_positions': {},
        'transcript_ids': set()
    }
    
    feature_counts = defaultdict(int)
    
    for feature in db.all_features():
        feature_counts[feature.featuretype] += 1
        
        if feature.featuretype in ['exon', 'CDS']:
            # Get transcript ID
            transcript_id = None
            if 'Parent' in feature.attributes:
                parent = feature.attributes['Parent'][0]
                if 'transcript:' in parent:
                    transcript_id = parent.split('transcript:')[1]
                else:
                    transcript_id = parent
            
            if transcript_id:
                if feature.featuretype == 'exon':
                    transcript_info['exon_coords'][transcript_id].append((feature.start, feature.end))
                    length = sum(end - start for start, end in transcript_info['exon_coords'][transcript_id])
                    transcript_info['exon_lengths'][transcript_id] = length
                else:  # CDS
                    transcript_info['cds_coords'][transcript_id].append((feature.start, feature.end))
                    if feature.strand == '+':
                        transcript_info['mAUG_positions'][transcript_id] = min(start for start, _ in transcript_info['cds_coords'][transcript_id])
                    else:
                        transcript_info['mAUG_positions'][transcript_id] = max(end for _, end in transcript_info['cds_coords'][transcript_id])
                    
                    length = sum(end - start for start, end in transcript_info['cds_coords'][transcript_id])
                    transcript_info['cds_lengths'][transcript_id] = length
                
                transcript_info['transcript_ids'].add(transcript_id)
    
    print("\nFeature counts from genome annotation:")
    for ftype, count in feature_counts.items():
        print(f"{ftype}: {count}")
    
    print(f"\nLoaded {len(transcript_info['transcript_ids'])} transcripts")
    return transcript_info

def calculate_rpkm(bam_file, feature_coordinates, feature_lengths):
    """Calculate RPKM for given genomic features"""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        
        # Get total mapped reads
        total_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        
        # Initialize counts
        counts = defaultdict(int)
        
        # Count reads per feature
        for transcript_id, regions in feature_coordinates.items():
            for start, end in regions:
                try:
                    # Get chromosome from transcript ID (e.g., AT1G01020.1 -> 1)
                    chrom = transcript_id[2]  # Gets the chromosome number
                    if not chrom.isdigit():
                        continue
                        
                    # Count reads in feature region
                    for read in bam.fetch(chrom, start, end):
                        if not read.is_unmapped:
                            counts[transcript_id] += 1
                except ValueError:
                    continue
        
        # Calculate RPKM
        rpkm = {}
        for transcript_id in counts:
            if transcript_id in feature_lengths and feature_lengths[transcript_id] > 0:
                rpkm[transcript_id] = (counts[transcript_id] * 1e9) / (total_reads * feature_lengths[transcript_id])
        
        bam.close()
        return rpkm
    
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
        return {}

def analyze_data_distribution(transcript_info):
    """Analyze data distribution to determine appropriate cutoffs"""
    print("\nAnalyzing data distribution...")
    
    rna_files = [
        os.path.join(BASE_DIR, "unique_reads", f"LZT101-{i}_uniq_sort.bam") for i in [1, 2]
    ] + [
        os.path.join(BASE_DIR, "unique_reads", f"LZT102-{i}_uniq_sort.bam") for i in [1, 2]
    ]
    
    ribo_files = [
        os.path.join(BASE_DIR, "unique_reads", f"LZT103-{i}_uniq_sort.bam") for i in [1, 2]
    ] + [
        os.path.join(BASE_DIR, "unique_reads", f"LZT104-{i}_uniq_sort.bam") for i in [1, 2]
    ]
    
    # Calculate RPKM distributions
    print("\nProcessing RNA-seq files...")
    exon_rpkm_values = []
    for rna_file in rna_files:
        print(f"Processing {os.path.basename(rna_file)}")
        rpkm = calculate_rpkm(rna_file, transcript_info['exon_coords'], 
                            transcript_info['exon_lengths'])
        exon_rpkm_values.extend(rpkm.values())
    
    print("\nProcessing Ribo-seq files...")
    cds_rpkm_values = []
    for ribo_file in ribo_files:
        print(f"Processing {os.path.basename(ribo_file)}")
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        cds_rpkm_values.extend(rpkm.values())
    
    # Calculate percentiles
    percentiles = [10, 25, 50, 75, 90, 95]
    distribution_stats = {
        'RNA_exon_percentiles': {p: np.percentile(exon_rpkm_values, p) 
                                for p in percentiles},
        'Ribo_CDS_percentiles': {p: np.percentile(cds_rpkm_values, p) 
                                for p in percentiles}
    }
    
    return distribution_stats

def main():
    """Main analysis workflow"""
    print("Starting analysis...")
    
    # Create/load genome database
    genome_db = create_genome_db()
    
    # Load transcript annotations
    transcript_info = load_transcript_annotations(genome_db)
    
    # Analyze distributions
    distribution_stats = analyze_data_distribution(transcript_info)
    
    # Print statistics
    print("\nData Distribution Statistics:")
    print("\nRNA-seq RPKM percentiles:")
    for p, value in distribution_stats['RNA_exon_percentiles'].items():
        print(f"{p}th percentile: {value:.2f}")
    
    print("\nRibo-seq RPKM percentiles:")
    for p, value in distribution_stats['Ribo_CDS_percentiles'].items():
        print(f"{p}th percentile: {value:.2f}")

if __name__ == "__main__":
    main()
