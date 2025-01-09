#!/usr/bin/env python3

import sys
import os

def clean_fasta(input_path, output_path):
    """Clean FASTA file to ensure proper formatting."""
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        is_first = True
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Simplify header to just chromosome number/name
                header = line.split()[0]  # Take first part of header
                if not is_first:
                    outfile.write('\n')
                outfile.write(header + '\n')
                is_first = False
            else:
                outfile.write(line)

def main():
    input_genome = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
    output_genome = "/global/scratch/users/enricocalvane/riboseq/imb2/ribotish/reference/cleaned_genome.fa"
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_genome), exist_ok=True)
    
    print("Cleaning FASTA file...")
    clean_fasta(input_genome, output_genome)
    print("FASTA cleaning complete")

    # Print first few lines of cleaned file
    print("\nFirst few lines of cleaned FASTA:")
    with open(output_genome, 'r') as f:
        for i, line in enumerate(f):
            if i < 10:
                print(line.strip())
            else:
                break

if __name__ == '__main__':
    main()
