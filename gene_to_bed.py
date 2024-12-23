#!/usr/bin/env python3

import sys
import argparse

def parse_gtf(gtf_file):
    """Parse GTF file and extract gene coordinates using gene_id field"""
    gene_coords = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            # Extract gene ID from the GTF attributes field
            attributes = parts[8]
            gene_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    # Extract the gene ID between quotes and remove "gene:" prefix
                    gene_id = attr.split('"')[1].replace("gene:", "")
                    break
                    
            if gene_id and gene_id not in gene_coords:
                gene_coords[gene_id] = {
                    'chr': parts[0],
                    'start': parts[3],
                    'end': parts[4],
                    'strand': parts[6]
                }

    # Print statistics for debugging
    print(f"Debug - Total genes found in GTF: {len(gene_coords)}", file=sys.stderr)
    if gene_coords:
        print(f"Debug - Example gene IDs in GTF: {list(gene_coords.keys())[:5]}", file=sys.stderr)
    
    return gene_coords

def create_bed(gene_list_file, gtf_file, output_file):
    """Create BED file from gene list using coordinates from GTF"""
    # Parse GTF file
    gene_coords = parse_gtf(gtf_file)
    
    # Read gene list
    with open(gene_list_file) as f:
        genes = [line.strip() for line in f if line.strip()]
    
    # Print first few genes from input list for debugging
    print(f"Debug - First few genes from input list: {genes[:5]}", file=sys.stderr)
    
    # Count matches and mismatches
    matches = 0
    mismatches = 0
    
    # Write BED file
    with open(output_file, 'w') as out:
        for gene in genes:
            if gene in gene_coords:
                coords = gene_coords[gene]
                bed_line = f"{coords['chr']}\t{coords['start']}\t{coords['end']}\t{gene}\t.\t{coords['strand']}\n"
                out.write(bed_line)
                matches += 1
            else:
                print(f"Warning: Gene {gene} not found in GTF file", file=sys.stderr)
                mismatches += 1
    
    print(f"Summary: {matches} genes found, {mismatches} genes not found", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Convert TAIR gene IDs to BED format')
    parser.add_argument('gene_list', help='File containing list of TAIR gene IDs')
    parser.add_argument('gtf_file', help='TAIR GTF annotation file')
    parser.add_argument('output_bed', help='Output BED file name')
    
    args = parser.parse_args()
    
    create_bed(args.gene_list, args.gtf_file, args.output_bed)

if __name__ == '__main__':
    main()
