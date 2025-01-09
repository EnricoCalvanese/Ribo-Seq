#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GTF attribute string to dictionary"""
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            key, value = attr.strip().split(' "', 1)
            attrs[key.strip()] = value.strip('"')
    return attrs

def main():
    # Store features by transcript
    transcripts = defaultdict(lambda: {'exons': [], 'cds': [], 'strand': None})
    
    # Read GTF file
    with open('/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
                
            chrom, _, feature, start, end, _, strand, _, attr_str = fields
            start = int(start)
            end = int(end)
            
            # Parse attributes
            attrs = parse_attributes(attr_str)
            if 'transcript_id' not in attrs:
                continue
                
            transcript_id = attrs['transcript_id']
            transcripts[transcript_id]['strand'] = strand
            
            if feature == 'exon':
                transcripts[transcript_id]['exons'].append((chrom, start, end))
            elif feature == 'CDS':
                transcripts[transcript_id]['cds'].append((chrom, start, end))

    # Output file
    with open('reference/tair10_5utr.bed', 'w') as out:
        for transcript_id, features in transcripts.items():
            if not features['exons'] or not features['cds']:
                continue
                
            strand = features['strand']
            chrom = features['exons'][0][0]
            
            if strand == '+':
                # Sort by start position
                first_exon_start = min(start for _, start, _ in features['exons'])
                first_cds_start = min(start for _, start, _ in features['cds'])
                
                if first_exon_start < first_cds_start:
                    # BED format: 0-based start, 1-based end
                    out.write(f'{chrom}\t{first_exon_start-1}\t{first_cds_start-1}\t{transcript_id}\t.\t{strand}\n')
            
            else:  # strand == '-'
                # Sort by end position
                last_exon_end = max(end for _, _, end in features['exons'])
                last_cds_end = max(end for _, _, end in features['cds'])
                
                if last_cds_end < last_exon_end:
                    # BED format: 0-based start, 1-based end
                    out.write(f'{chrom}\t{last_cds_end-1}\t{last_exon_end}\t{transcript_id}\t.\t{strand}\n')

if __name__ == '__main__':
    main()
