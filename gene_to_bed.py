#!/usr/bin/env python3

import sys
import argparse

def parse_gff(gff_file):
    """Parse GFF file and extract gene coordinates"""
    gene_coords = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
                
            # Extract gene ID from the attributes field
            attributes = dict(attr.split('=') for attr in parts[8].split(';') if '=' in attr)
            gene_id = attributes.get('ID', '').replace('gene:', '')
            
            if gene_id:
                gene_coords[gene_id] = {
                    'chr': parts[0],
                    'start': parts[3],
                    'end': parts[4],
                    'strand': parts[6]
                }
    return gene_coords

def create_bed(gene_list_file, gff_file, output_file):
    """Create BED file from gene list using coordinates from GFF"""
    # Parse GFF file
    gene_coords = parse_gff(gff_file)
    
    # Read gene list
    with open(gene_list_file) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    # Write BED file
    with open(output_file, 'w') as out:
        for gene in genes:
            if gene in gene_coords:
                coords = gene_coords[gene]
                bed_line = f"{coords['chr']}\t{coords['start']}\t{coords['end']}\t{gene}\t.\t{coords['strand']}\n"
                out.write(bed_line)
            else:
                print(f"Warning: Gene {gene} not found in GFF file", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Convert TAIR gene IDs to BED format')
    parser.add_argument('gene_list', help='File containing list of TAIR gene IDs')
    parser.add_argument('gff_file', help='TAIR GFF annotation file')
    parser.add_argument('output_bed', help='Output BED file name')
    
    args = parser.parse_args()
    
    create_bed(args.gene_list, args.gff_file, args.output_bed)

if __name__ == '__main__':
    main()
