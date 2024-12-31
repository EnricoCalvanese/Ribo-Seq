import pandas as pd
import numpy as np
import pysam
import gffutils
from collections import defaultdict
import os

# File paths
BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")

# Sample information
SAMPLE_INFO = {
    'LZT101-1': {'type': 'RNA-Seq', 'condition': 'WT', 'replicate': 1},
    'LZT101-2': {'type': 'RNA-Seq', 'condition': 'WT', 'replicate': 2},
    'LZT102-1': {'type': 'RNA-Seq', 'condition': 'imb2', 'replicate': 1},
    'LZT102-2': {'type': 'RNA-Seq', 'condition': 'imb2', 'replicate': 2},
    'LZT103-1': {'type': 'Ribo-Seq', 'condition': 'WT', 'replicate': 1},
    'LZT103-2': {'type': 'Ribo-Seq', 'condition': 'WT', 'replicate': 2},
    'LZT104-1': {'type': 'Ribo-Seq', 'condition': 'imb2', 'replicate': 1},
    'LZT104-2': {'type': 'Ribo-Seq', 'condition': 'imb2', 'replicate': 2}
}

def get_bam_files():
    """Get sorted lists of RNA-seq and Ribo-seq BAM files"""
    bam_dir = os.path.join(BASE_DIR, "unique_reads")
    rna_files = []
    ribo_files = []
    
    for sample_id, info in SAMPLE_INFO.items():
        bam_file = os.path.join(bam_dir, f"{sample_id}_uniq_sort.bam")
        if info['type'] == 'RNA-Seq':
            rna_files.append(bam_file)
        else:
            ribo_files.append(bam_file)
    
    return sorted(rna_files), sorted(ribo_files)

def parse_uorf_data():
    """Parse uORF GFF file and create transcript annotation database"""
    # Create temporary SQLite database
    db_path = "temp_uorf.db"
    if os.path.exists(db_path):
        os.remove(db_path)
    
    # Parse GFF file
    db = gffutils.create_db(
        UORF_GFF,
        db_path,
        force=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )
    
    uorf_data = defaultdict(list)
    transcript_info = {
        'exon_coords': {},
        'cds_coords': {},
        'exon_lengths': {},
        'cds_lengths': {},
        'mAUG_positions': {},
        'transcript_ids': set(),
        'five_prime_leaders': {}
    }
    
    # Extract uORF and transcript information
    for feature in db.all_features():
        if feature.featuretype == 'uORF':
            transcript_id = feature.attributes.get('Parent', [None])[0]
            if transcript_id:
                uorf_data[transcript_id].append({
                    'start': feature.start,
                    'end': feature.end,
                    'strand': feature.strand,
                    'score': float(feature.score) if feature.score != '.' else 0
                })
        
        # Extract transcript features
        elif feature.featuretype in ['exon', 'CDS']:
            transcript_id = feature.attributes.get('Parent', [None])[0]
            if transcript_id:
                if feature.featuretype == 'exon':
                    transcript_info['exon_coords'][transcript_id] = (feature.start, feature.end)
                    transcript_info['exon_lengths'][transcript_id] = feature.end - feature.start
                else:
                    transcript_info['cds_coords'][transcript_id] = (feature.start, feature.end)
                    transcript_info['cds_lengths'][transcript_id] = feature.end - feature.start
                    # Assuming mAUG is at CDS start
                    transcript_info['mAUG_positions'][transcript_id] = feature.start
                
                transcript_info['transcript_ids'].add(transcript_id)
    
    os.remove(db_path)
    return uorf_data, transcript_info

def analyze_data_distribution():
    """Analyze data distribution to determine appropriate cutoffs"""
    rna_files, ribo_files = get_bam_files()
    uorf_data, transcript_info = parse_uorf_data()
    
    # Calculate RPKM distributions
    exon_rpkm = defaultdict(list)
    cds_rpkm = defaultdict(list)
    
    for rna_file in rna_files:
        rpkm = calculate_rpkm(rna_file, transcript_info['exon_coords'], 
                            transcript_info['exon_lengths'])
        for value in rpkm.values():
            exon_rpkm['RNA'].append(value)
    
    for ribo_file in ribo_files:
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        for value in rpkm.values():
            cds_rpkm['Ribo'].append(value)
    
    # Calculate percentiles
    percentiles = [25, 50, 75, 90, 95]
    distribution_stats = {
        'RNA_exon_percentiles': {p: np.percentile(exon_rpkm['RNA'], p) 
                                for p in percentiles},
        'Ribo_CDS_percentiles': {p: np.percentile(cds_rpkm['Ribo'], p) 
                                for p in percentiles}
    }
    
    return distribution_stats

def calculate_rpkm(bam_file, feature_coordinates, feature_lengths):
    """Calculate RPKM for given genomic features"""
    # Implementation remains the same as in previous version
    pass

def main():
    """Main analysis workflow"""
    # Get file paths
    rna_files, ribo_files = get_bam_files()
    
    # Parse uORF and transcript data
    uorf_data, transcript_info = parse_uorf_data()
    
    # Analyze data distribution for cutoff determination
    distribution_stats = analyze_data_distribution()
    
    # Print initial statistics
    print("Data Distribution Statistics:")
    print("RNA-seq RPKM percentiles:", distribution_stats['RNA_exon_percentiles'])
    print("Ribo-seq RPKM percentiles:", distribution_stats['Ribo_CDS_percentiles'])
    
    # Based on distribution, we can adjust these cutoffs
    rpkm_cutoff = 1  # We might want to adjust this based on distribution
    background_cutoff = None  # Will be calculated based on data
    min_raw_count = 10
    
    # Run main analysis
    results = analyze_augs(rna_files, ribo_files, transcript_info, uorf_data)
    
    return results

if __name__ == "__main__":
    main()
