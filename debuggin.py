import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
from typing import List, Dict, Set, Tuple

def get_chromosome_from_transcript(transcript_id: str) -> str:
    """
    Extracts chromosome number from Arabidopsis transcript ID.
    
    Args:
        transcript_id: TAIR-style transcript ID (e.g., 'AT1G01020.2')
        
    Returns:
        Chromosome identifier (e.g., '1')
    """
    # Remove 'transcript:' prefix if present
    if transcript_id.startswith('transcript:'):
        transcript_id = transcript_id[len('transcript:'):]
    
    # Extract chromosome number using regex
    match = re.match(r'AT(\d|C|M)G', transcript_id)
    if match:
        chr_id = match.group(1)
        # Convert chromosome identifiers
        chr_map = {
            'C': 'Pt',  # Chloroplast -> Plastid
            'M': 'Mt'   # Mitochondria
        }
        return chr_map.get(chr_id, str(chr_id))
    return None

def get_transcript_positions(db: gffutils.FeatureDB, 
                           transcript_id: str) -> Tuple[str, int, int]:
    """
    Gets chromosome and positions for a transcript.
    
    Args:
        db: GFF database
        transcript_id: Transcript identifier
        
    Returns:
        Tuple of (chromosome, start, end)
    """
    try:
        feature = db[transcript_id]
        chromosome = get_chromosome_from_transcript(transcript_id)
        return chromosome, feature.start, feature.end
    except Exception as e:
        print(f"Error getting positions for {transcript_id}: {str(e)}")
        return None, None, None


def count_upstream_reads(bam_file: str,
                        db: gffutils.FeatureDB,
                        transcripts: Set[str],
                        upstream_distance: int = 50) -> Dict[str, int]:
    """
    Counts reads in regions upstream of mAUGs using chromosome-level coordinates.
    
    Args:
        bam_file: Path to BAM file
        db: GFF database
        transcripts: Set of transcript IDs
        upstream_distance: Distance upstream to analyze
        
    Returns:
        Dictionary of read counts per transcript
    """
    counts = defaultdict(int)
    processed = 0
    errors = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get available chromosomes in BAM file
        valid_chromosomes = set(bam.references)
        print(f"\nValid chromosomes in BAM: {', '.join(valid_chromosomes)}")
        
        for transcript_id in transcripts:
            try:
                # Get chromosome and positions
                chr_id, start, end = get_transcript_positions(db, transcript_id)
                
                if chr_id not in valid_chromosomes:
                    continue
                
                # Get the CDS start position (using transcript feature)
                feature = db[transcript_id]
                cds_features = list(db.children(feature, featuretype='CDS'))
                
                if not cds_features:
                    continue
                    
                # Get first CDS position (mAUG)
                first_cds = min(cds_features, key=lambda x: x.start)
                cds_start = first_cds.start
                
                # Count reads in upstream window
                for read in bam.fetch(chr_id, 
                                    max(0, cds_start - upstream_distance),
                                    cds_start):
                    counts[transcript_id] += 1
                
                processed += 1
                if processed % 100 == 0:
                    print(f"Processed {processed} transcripts...")
                    
            except Exception as e:
                errors += 1
                if errors <= 5:
                    print(f"Error processing {transcript_id}: {str(e)}")
                continue
    
    print(f"\nSuccessfully processed {processed} transcripts")
    print(f"Encountered errors with {errors} transcripts")
    return counts


def examine_bam_references(bam_file: str, num_examples: int = 5) -> Set[str]:
    """
    Examines the reference names in a BAM file and prints examples.
    
    Args:
        bam_file: Path to BAM file
        num_examples: Number of example references to print
        
    Returns:
        Set of all reference names in the BAM file
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        references = set(bam.references)
        print(f"\nFound {len(references)} unique references in BAM file")
        print("Example reference names:")
        for ref in list(references)[:num_examples]:
            print(f"  {ref}")
    return references

def examine_transcript_ids(transcripts: Set[str], num_examples: int = 5) -> None:
    """
    Prints examples of transcript IDs from the GFF database.
    
    Args:
        transcripts: Set of transcript IDs
        num_examples: Number of examples to print
    """
    print("\nExample transcript IDs from GFF:")
    for transcript_id in list(transcripts)[:num_examples]:
        print(f"  {transcript_id}")

def clean_transcript_id(transcript_id: str) -> str:
    """
    Standardizes transcript ID format by removing prefixes and cleaning up the ID.
    
    Args:
        transcript_id: Raw transcript ID
        
    Returns:
        Cleaned transcript ID
    """
    # First, print the original ID for debugging
    print(f"Cleaning transcript ID: {transcript_id}")
    
    # Remove common prefixes
    prefixes_to_remove = ['transcript:', 'gene:', 'mRNA:']
    cleaned_id = transcript_id
    for prefix in prefixes_to_remove:
        if cleaned_id.startswith(prefix):
            cleaned_id = cleaned_id[len(prefix):]
    
    # Print the cleaned ID
    print(f"Cleaned to: {cleaned_id}")
    return cleaned_id

def verify_bam_references(bam_file: str, transcripts: Set[str]) -> Set[str]:
    """
    Verifies which transcripts are present in the BAM file.
    Now includes detailed debugging information.
    
    Args:
        bam_file: Path to BAM file
        transcripts: Set of transcript IDs to verify
        
    Returns:
        Set of valid transcript IDs
    """
    valid_transcripts = set()
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_references = set(bam.references)
        print(f"\nBAM file contains {len(bam_references)} references")
        
        # Print some example BAM references
        print("\nExample BAM references:")
        for ref in list(bam_references)[:5]:
            print(f"  {ref}")
        
        # Check each transcript
        print("\nChecking transcript IDs against BAM references...")
        for transcript_id in transcripts:
            clean_id = clean_transcript_id(transcript_id)
            if clean_id in bam_references:
                valid_transcripts.add(clean_id)
            elif len(valid_transcripts) == 0:  # Only print first few failures
                print(f"Failed to find: {clean_id}")
    
    return valid_transcripts

def main():
    """
    Main function orchestrating the analysis of ribosome footprints upstream of main AUG codons.
    The analysis follows these steps:
    1. Load and validate input files
    2. Process transcript information from GFF database
    3. Identify eligible transcripts based on leader sequences and uORFs
    4. Count and normalize upstream reads for each Ribo-seq file
    5. Calculate statistics and save results
    """
    # Define file paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    UORF_GFF = os.path.join(BASE_DIR, "systemPipeR/uorf.gff")
    
    RIBO_SEQ_FILES = [
        os.path.join(BASE_DIR, "unique_reads/LZT103-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT103-2_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-1_uniq_sort.bam"),
        os.path.join(BASE_DIR, "unique_reads/LZT104-2_uniq_sort.bam")
    ]
    
    # Step 1: Load and validate input files
    print("\nStep 1: Loading input files...")
    try:
        db = gffutils.FeatureDB(GFF3_FILE + '.db')
        print("Successfully loaded GFF database")
    except Exception as e:
        print(f"Error loading GFF database: {str(e)}")
        return
    
    # Step 2: Load transcript information
    print("\nStep 2: Processing transcript information...")
    print("Loading uORF positions...")
    uorf_positions = load_uorf_positions(UORF_GFF)
    print(f"Found uORFs in {len(uorf_positions)} transcripts")
    
    print("Calculating transcript leader lengths...")
    leader_lengths = get_transcript_lengths(db)
    print(f"Calculated leader lengths for {len(leader_lengths)} transcripts")
    
    # Step 3: Identify eligible transcripts
    print("\nStep 3: Identifying eligible transcripts...")
    eligible_transcripts = get_eligible_transcripts(
        leader_lengths, 
        uorf_positions,
        min_leader_length=100  # Minimum leader length requirement
    )
    print(f"Found {len(eligible_transcripts)} eligible transcripts")
    
    # Step 4: Process each Ribo-seq file
    print("\nStep 4: Processing Ribo-seq files...")
    results = []
    
    for bam_file in RIBO_SEQ_FILES:
        sample_name = os.path.basename(bam_file)
        print(f"\nProcessing {sample_name}...")
        
        try:
            # Get total reads for normalization
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                total_reads = bam.count()
                print(f"Total mapped reads: {total_reads:,}")
            
            # Count upstream reads using chromosome-aware function
            counts = count_upstream_reads(
                bam_file=bam_file,
                db=db,
                transcripts=eligible_transcripts,
                upstream_distance=50  # 50nt upstream window
            )
            
            # Skip processing if no counts were found
            if not counts:
                print(f"Warning: No valid counts found in {sample_name}")
                continue
            
            # Calculate transcript abundances (using placeholder values for now)
            # In practice, you would calculate this from your RNA-seq data
            transcript_abundances = {tid: 1.0 for tid in eligible_transcripts}
            
            # Normalize counts
            normalized_counts = normalize_counts(
                counts=counts,
                total_reads=total_reads,
                transcript_abundances=transcript_abundances
            )
            
            # Calculate statistics
            stats = calculate_statistics(normalized_counts)
            
            # Store results for this sample
            sample_result = {
                'sample': sample_name,
                'normalized_counts': normalized_counts,
                'total_reads': total_reads,
                'stats': stats
            }
            results.append(sample_result)
            
            # Print sample statistics
            print(f"\nStatistics for {sample_name}:")
            for stat_name, stat_value in stats.items():
                if stat_value is not None:
                    print(f"{stat_name.capitalize()}: {stat_value:.2f}")
                else:
                    print(f"{stat_name.capitalize()}: N/A")
                    
        except Exception as e:
            print(f"Error processing {sample_name}: {str(e)}")
            continue
    
    # Step 5: Save results
    print("\nStep 5: Saving results...")
    if not results:
        print("Warning: No results to save")
        return
        
    try:
        # Prepare results for output
        output_data = []
        for result in results:
            for transcript_id, count in result['normalized_counts'].items():
                output_data.append({
                    'sample': result['sample'],
                    'transcript_id': transcript_id,
                    'normalized_count': count,
                    'q3': result['stats']['q3'],
                    'median': result['stats']['median'],
                    'mean': result['stats']['mean']
                })
        
        # Save to file
        output_file = os.path.join(BASE_DIR, "systemPipeR/maug_upstream_analysis.tsv")
        pd.DataFrame(output_data).to_csv(output_file, sep='\t', index=False)
        print(f"Results saved to {output_file}")
        
        # Print summary statistics
        print("\nAnalysis Summary:")
        print(f"Total transcripts analyzed: {len(eligible_transcripts)}")
        print(f"Samples processed: {len(results)}")
        print("Done!")
        
    except Exception as e:
        print(f"Error saving results: {str(e)}")

if __name__ == "__main__":
    main()
