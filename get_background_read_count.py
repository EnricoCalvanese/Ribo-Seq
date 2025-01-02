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

def load_uorf_positions(uorf_gff: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Reads and parses the uORF GFF file to identify transcripts containing upstream AUG codons.
    This function creates a mapping between transcript IDs and their uORF positions.
    
    Args:
        uorf_gff: Path to the GFF file containing uORF annotations
        
    Returns:
        Dictionary where keys are transcript IDs and values are lists of (start, end) positions
    """
    uorf_positions = defaultdict(list)
    print(f"Reading uORF file: {uorf_gff}")
    
    try:
        with open(uorf_gff) as f:
            for line in f:
                # Skip comment lines in GFF file
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                # Ensure line has all required GFF fields
                if len(fields) < 9:
                    continue
                    
                transcript_id = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                uorf_positions[transcript_id].append((start, end))
                
        print(f"Found uORFs in {len(uorf_positions)} transcripts")
        return uorf_positions
        
    except FileNotFoundError:
        print(f"Error: uORF GFF file not found at {uorf_gff}")
        return {}
    except Exception as e:
        print(f"Error processing uORF file: {str(e)}")
        return {}

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

def get_transcript_lengths(db: gffutils.FeatureDB) -> Dict[str, int]:
    """
    Calculates the length of 5' leader sequences for all transcripts in the database.
    The leader sequence is defined as the region from transcript start to the first CDS.
    
    Args:
        db: GFF database containing transcript and CDS annotations
        
    Returns:
        Dictionary mapping transcript IDs to their 5' leader lengths
    """
    leader_lengths = {}
    processed = 0
    
    for transcript in db.features_of_type('mRNA'):
        try:
            # Find all CDS features for this transcript
            cds_features = list(db.children(transcript, featuretype='CDS'))
            
            if cds_features:
                # Get the start position of the first CDS
                first_cds = min(cds_features, key=lambda x: x.start)
                # Calculate leader length
                leader_length = first_cds.start - transcript.start
                leader_lengths[transcript.id] = leader_length
            
            processed += 1
            if processed % 1000 == 0:
                print(f"Processed {processed} transcripts...")
                
        except Exception as e:
            print(f"Warning: Could not process transcript {transcript.id}: {str(e)}")
            continue
    
    print(f"Calculated leader lengths for {len(leader_lengths)} transcripts")
    return leader_lengths

def get_eligible_transcripts(leader_lengths: Dict[str, int],
                           uorf_positions: Dict[str, List[Tuple[int, int]]],
                           min_leader_length: int = 100) -> Set[str]:
    """
    Identifies transcripts that meet our criteria:
    1. Have a 5' leader sequence ≥ 100 nucleotides
    2. Do not contain any upstream AUG codons
    
    Args:
        leader_lengths: Dictionary of transcript leader lengths
        uorf_positions: Dictionary of uORF positions
        min_leader_length: Minimum required leader length (default: 100)
        
    Returns:
        Set of transcript IDs that meet all criteria
    """
    eligible = set()
    
    for transcript_id, length in leader_lengths.items():
        # Check if leader is long enough
        if length >= min_leader_length:
            # Check if transcript has no uORFs
            if transcript_id not in uorf_positions:
                eligible.add(transcript_id)
    
    print(f"Found {len(eligible)} eligible transcripts out of {len(leader_lengths)} total")
    return eligible
                               
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
       
def clean_transcript_id(transcript_id: str) -> str:
    """
    Standardizes transcript ID format by removing prefixes and cleaning up the ID.
    
    For example:
    - "transcript:AT5G53480.1" becomes "AT5G53480.1"
    - "AT5G53480.1" remains "AT5G53480.1"
    
    Args:
        transcript_id: Raw transcript ID from GFF or BAM file
        
    Returns:
        Cleaned transcript ID
    """
    # Remove common prefixes
    prefixes_to_remove = ['transcript:', 'gene:', 'mRNA:']
    cleaned_id = transcript_id
    for prefix in prefixes_to_remove:
        if cleaned_id.startswith(prefix):
            cleaned_id = cleaned_id[len(prefix):]
    return cleaned_id

def verify_bam_references(bam_file: str, transcripts: Set[str]) -> Set[str]:
    """
    Verifies which transcripts are actually present in the BAM file and
    returns the set of valid transcripts.
    
    Args:
        bam_file: Path to BAM file
        transcripts: Set of transcript IDs to verify
        
    Returns:
        Set of transcript IDs that exist in the BAM file
    """
    valid_transcripts = set()
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get all reference names from BAM file
        bam_references = set(bam.references)
        
        # Check each transcript
        for transcript_id in transcripts:
            clean_id = clean_transcript_id(transcript_id)
            if clean_id in bam_references:
                valid_transcripts.add(clean_id)
    
    print(f"Found {len(valid_transcripts)} valid transcripts in BAM file")
    return valid_transcripts

def normalize_counts(counts: Dict[str, int],
                    total_reads: int,
                    transcript_abundances: Dict[str, float]) -> Dict[str, float]:
    """
    Normalizes read counts by total sequencing depth and transcript abundance.
    
    Args:
        counts: Raw read counts per transcript
        total_reads: Total number of mapped reads in the sample
        transcript_abundances: RPKM values for each transcript
        
    Returns:
        Dictionary of normalized read counts
    """
    normalized = {}
    
    for transcript_id, count in counts.items():
        if transcript_id in transcript_abundances and transcript_abundances[transcript_id] > 0:
            # Normalize by total reads (RPM) and transcript abundance
            normalized[transcript_id] = (count / total_reads * 1e6) / transcript_abundances[transcript_id]
    
    return normalized

def calculate_statistics(normalized_counts: Dict[str, float]) -> Dict[str, float]:
    """
    Calculates statistics with improved error handling and validation.
    """
    stats = {'q3': None, 'median': None, 'mean': None}
    
    if not normalized_counts:
        print("Warning: No normalized counts available for statistics calculation")
        return stats
    
    counts_list = list(normalized_counts.values())
    
    if not counts_list:
        print("Warning: Empty counts list")
        return stats
    
    try:
        # Calculate statistics only if we have valid data
        if len(counts_list) > 0:
            stats['q3'] = float(np.percentile(counts_list, 75))
            stats['median'] = float(np.median(counts_list))
            stats['mean'] = float(np.mean(counts_list))
            
            # Validate results
            if not all(isinstance(v, (int, float)) for v in [stats['q3'], stats['median'], stats['mean']]):
                print("Warning: Invalid statistics calculated")
                stats = {'q3': None, 'median': None, 'mean': None}
    except Exception as e:
        print(f"Warning: Error calculating statistics: {str(e)}")
        stats = {'q3': None, 'median': None, 'mean': None}
    
    return stats

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
