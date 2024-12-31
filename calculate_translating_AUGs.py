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

def calculate_rpkm(bam_file, feature_coordinates, feature_lengths):
    """
    Calculate RPKM for given genomic features
    
    Parameters:
    bam_file (str): Path to BAM file
    feature_coordinates (dict): Dictionary of feature coordinates {feature_id: (start, end)}
    feature_lengths (dict): Dictionary of feature lengths {feature_id: length}
    
    Returns:
    dict: Dictionary of RPKM values for each feature
    """
    try:
        # Open BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
        
        # Get total mapped reads
        total_reads = float(sum(1 for read in bam.fetch() if not read.is_unmapped))
        
        # Initialize counts
        counts = defaultdict(int)
        
        # Count reads per feature
        for feature_id, (start, end) in feature_coordinates.items():
            try:
                # Get chromosome/contig name from the feature_id if needed
                chrom = feature_id.split(':')[0] if ':' in feature_id else feature_id
                
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
                # RPKM = (reads * 10^9) / (total_reads * feature_length)
                rpkm[feature_id] = (counts[feature_id] * 1e9) / (total_reads * feature_lengths[feature_id])
        
        bam.close()
        return rpkm
    
    except Exception as e:
        print(f"Error processing {bam_file}: {str(e)}")
        return {}

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
    print("Parsing uORF GFF file...")
    
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
    
    print("Extracting features...")
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
                    transcript_info['mAUG_positions'][transcript_id] = feature.start
                
                transcript_info['transcript_ids'].add(transcript_id)
    
    print(f"Found {len(uorf_data)} transcripts with uORFs")
    print(f"Found {len(transcript_info['transcript_ids'])} total transcripts")
    
    os.remove(db_path)
    return uorf_data, transcript_info

def analyze_data_distribution():
    """Analyze data distribution to determine appropriate cutoffs"""
    print("Analyzing data distribution...")
    rna_files, ribo_files = get_bam_files()
    uorf_data, transcript_info = parse_uorf_data()
    
    # Calculate RPKM distributions
    exon_rpkm = defaultdict(list)
    cds_rpkm = defaultdict(list)
    
    print("Processing RNA-seq files...")
    for rna_file in rna_files:
        print(f"Processing {os.path.basename(rna_file)}")
        rpkm = calculate_rpkm(rna_file, transcript_info['exon_coords'], 
                            transcript_info['exon_lengths'])
        exon_rpkm['RNA'].extend(rpkm.values())
    
    print("Processing Ribo-seq files...")
    for ribo_file in ribo_files:
        print(f"Processing {os.path.basename(ribo_file)}")
        rpkm = calculate_rpkm(ribo_file, transcript_info['cds_coords'],
                            transcript_info['cds_lengths'])
        cds_rpkm['Ribo'].extend(rpkm.values())
    
    # Calculate percentiles
    percentiles = [25, 50, 75, 90, 95]
    distribution_stats = {
        'RNA_exon_percentiles': {p: np.percentile(exon_rpkm['RNA'], p) 
                                for p in percentiles},
        'Ribo_CDS_percentiles': {p: np.percentile(cds_rpkm['Ribo'], p) 
                                for p in percentiles}
    }
    
    return distribution_stats

def main():
    """Main analysis workflow"""
    print("Starting analysis...")
    
    # Get file paths
    rna_files, ribo_files = get_bam_files()
    print(f"Found {len(rna_files)} RNA-seq files and {len(ribo_files)} Ribo-seq files")
    
    # Analyze data distribution for cutoff determination
    distribution_stats = analyze_data_distribution()
    
    # Print initial statistics
    print("\nData Distribution Statistics:")
    print("RNA-seq RPKM percentiles:", distribution_stats['RNA_exon_percentiles'])
    print("Ribo-seq RPKM percentiles:", distribution_stats['Ribo_CDS_percentiles'])

if __name__ == "__main__":
    main()
