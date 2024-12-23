#!/usr/bin/env python3

import sys
import argparse

def parse_gtf(gtf_file):
    """Parse GTF file and extract gene coordinates"""
    gene_coords = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
                
            # Extract gene ID from the GTF attributes field
            attributes = parts[8]
            # Print first few attributes for debugging
            if len(gene_coords) == 0:
                print(f"Debug - First gene attributes: {attributes}", file=sys.stderr)
                
            # Handle different possible formats of gene IDs in GTF
            gene_id = None
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1]
                    break
                    
            if gene_id:
                # Remove any prefix/suffix to match the basic TAIR ID format
                gene_id = gene_id.replace("gene:", "").strip()
                gene_coords[gene_id] = {
                    'chr': parts[0],
                    'start': parts[3],
                    'end': parts[4],
                    'strand': parts[6]
                }
                
                # Print first gene found for debugging
                if len(gene_coords) == 1:
                    print(f"Debug - First gene parsed: {gene_id}", file=sys.stderr)
    
    # Print number of genes found
    print(f"Debug - Total genes found in GTF: {len(gene_coords)}", file=sys.stderr)
    # Print a few example gene IDs
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
            clean_gene = gene.strip().replace("gene:", "")
            if clean_gene in gene_coords:
                coords = gene_coords[clean_gene]
                bed_line = f"{coords['chr']}\t{coords['start']}\t{coords['end']}\t{clean_gene}\t.\t{coords['strand']}\n"
                out.write(bed_line)
                matches += 1
            else:
                print(f"Warning: Gene {clean_gene} not found in GTF file", file=sys.stderr)
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
