#!/usr/bin/env python3

import sys
import subprocess
import os

def run_ribotish_quality(bam_path, output_prefix):
    """
    Run ribotish quality control with proper parameter handling.
    
    Parameters:
    bam_path: Path to the input BAM file
    output_prefix: Prefix for output files
    """
    # Define base command with parameters we know work
    cmd = [
        'ribotish', 'quality',
        '-b', bam_path,
        '-g', '/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gtf',
        '-o', f'{output_prefix}_qual.txt',
        '-f', f'{output_prefix}_qual.pdf',
        '-r', f'{output_prefix}.para.py',
        '-l', '25,35',
        '-p', '6',
        '-v'
    ]
    
    # Convert the -40,20 range into a proper format
    # We'll pass it as two separate numbers that ribotish can parse
    cmd.extend(['-d', '-40,20'])
    
    print(f"Running command: {' '.join(cmd)}")
    
    try:
        # Run the command and capture output
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Command output:")
        print(result.stdout)
        if result.stderr:
            print("Errors/Warnings:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running ribotish: {e}")
        print("Error output:")
        print(e.stderr)
        sys.exit(1)

def main():
    # Create output directory
    os.makedirs('quality_results', exist_ok=True)
    
    # Define our samples
    samples = [
        ('LZT103-1', 'WT_Rep1'),
        ('LZT103-2', 'WT_Rep2'),
        ('LZT104-1', 'imb2_Rep1'),
        ('LZT104-2', 'imb2_Rep2')
    ]
    
    # Process each sample
    for sample_id, description in samples:
        print(f"\nProcessing {sample_id} ({description})")
        bam_path = f'/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads/{sample_id}_uniq_sort.bam'
        output_prefix = f'quality_results/{sample_id}'
        run_ribotish_quality(bam_path, output_prefix)
    
    print("\nQuality control complete for all samples")

if __name__ == '__main__':
    main()
