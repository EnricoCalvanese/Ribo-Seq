import gffutils
import pyfaidx
from Bio import SeqIO
from Bio.Seq import Seq
import os
from collections import defaultdict

def create_db_from_gtf(gtf_file, db_file="arabidopsis.db"):
    """
    Create a gffutils database from GTF file if it doesn't exist.
    This database allows efficient querying of genomic features.
    """
    if not os.path.exists(db_file):
        print(f"Creating database from {gtf_file}...")
        db = gffutils.create_db(gtf_file, db_file, force=True, 
                              merge_strategy='merge',
                              id_spec={'transcript': ['transcript_id']},
                              gtf_transcript_key='transcript_id')
        print("Database creation complete!")
    return gffutils.FeatureDB(db_file)

def get_five_prime_leader(transcript, fasta, db):
    """
    Extract 5' leader sequence for a given transcript.
    Returns None if no valid leader sequence is found.
    """
    try:
        t = db[transcript]
        
        # Find the first CDS for this transcript
        cds_list = list(db.children(transcript, featuretype='CDS', order_by='start'))
        if not cds_list:
            return None
            
        if t.strand == '+':
            leader_start = t.start
            leader_end = cds_list[0].start - 1
        else:
            leader_start = cds_list[-1].end + 1
            leader_end = t.end
        
        # Verify leader sequence coordinates are valid
        if leader_end <= leader_start:
            return None
            
        # Extract sequence
        seq = fasta[t.seqid][leader_start-1:leader_end].seq
        if t.strand == '-':
            seq = str(Seq(seq).reverse_complement())
        
        # Verify sequence is not empty
        if not seq or len(seq) == 0:
            return None
            
        return seq
    except Exception as e:
        print(f"Error processing transcript {transcript}: {str(e)}")
        return None

def get_transcript_ids_for_gene(gene_id, db):
    """
    Get all transcript IDs associated with a gene ID.
    Handles both versioned and unversioned gene IDs.
    """
    # Remove version number if present
    base_gene_id = gene_id.split('.')[0]
    
    # Query the database for all transcripts of this gene
    transcript_ids = []
    try:
        # First try exact match
        for feature in db.children(gene_id, featuretype='transcript'):
            transcript_ids.append(feature.id)
        
        # If no exact match, try base gene ID
        if not transcript_ids:
            pattern = f"{base_gene_id}%"
            for feature in db.features_of_type('transcript', order_by='start'):
                if feature.id.startswith(base_gene_id):
                    transcript_ids.append(feature.id)
    except Exception as e:
        print(f"Warning: Error finding transcripts for gene {gene_id}: {str(e)}")
    
    return transcript_ids

def process_transcript_list(transcript_file, fasta_file, gtf_file, output_file):
    """
    Process a list of transcripts and write their 5' leader sequences to FASTA.
    Handles both file-based input and direct transcript ID input.
    Supports both versioned and unversioned gene IDs.
    """
    db = create_db_from_gtf(gtf_file)
    fasta = pyfaidx.Fasta(fasta_file)
    
    # If transcript_file is None, process all transcripts from the database
    if transcript_file is None:
        print("Processing all transcripts from GTF...")
        transcripts = [f.id for f in db.features_of_type('transcript')]
    else:
        print(f"Reading gene/transcript list from {transcript_file}...")
        with open(transcript_file) as f:
            input_ids = [line.strip() for line in f]
            
        # Expand gene IDs to their corresponding transcript IDs
        transcripts = []
        for input_id in input_ids:
            transcript_ids = get_transcript_ids_for_gene(input_id, db)
            if transcript_ids:
                transcripts.extend(transcript_ids)
            else:
                print(f"Warning: No transcripts found for {input_id}")
    
    print(f"Writing sequences to {output_file}...")
    processed = 0
    with open(output_file, 'w') as out:
        for transcript in transcripts:
            seq = get_five_prime_leader(transcript, fasta, db)
            if seq and len(seq) > 0:
                out.write(f">{transcript}\n{seq}\n")
                processed += 1
            if processed % 1000 == 0 and processed > 0:
                print(f"Processed {processed} transcripts...")
    
    print(f"Successfully processed {processed} transcripts with valid leader sequences")

def main():
    # Base paths
    base_dir = "/global/scratch/users/enricocalvane/riboseq/imb2/unique_reads"
    ref_dir = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference"
    
    # Reference files
    gtf_file = os.path.join(ref_dir, "Arabidopsis_thaliana.TAIR10.60.gtf")
    fasta_file = os.path.join(ref_dir, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.join(base_dir, "leader_sequences")
    os.makedirs(output_dir, exist_ok=True)
    
    # Input files for active uORFs
    uorf_files = {
        "down": "active_uORF_transcripts_TE-down.txt",
        "up": "active_uORF_transcripts_TE-up.txt",
        "nc": "active_uORF_transcripts_TE-nc.txt"
    }
    
    # Input files for translatome
    translatome_files = {
        "down": "TEdown_genes.txt",
        "up": "TEup_genes.txt",
        "nc": "TEnc_genes.txt"
    }
    
    print("Processing active uORF sequences...")
    for category, filename in uorf_files.items():
        input_file = os.path.join(base_dir, filename)
        output_file = os.path.join(output_dir, f"uORF_{category}_leaders.fa")
        process_transcript_list(input_file, fasta_file, gtf_file, output_file)
    
    print("\nProcessing translatome sequences...")
    for category, filename in translatome_files.items():
        input_file = os.path.join(base_dir, filename)
        output_file = os.path.join(output_dir, f"translatome_{category}_leaders.fa")
        process_transcript_list(input_file, fasta_file, gtf_file, output_file)
    
    print("\nProcessing full Arabidopsis transcriptome...")
    full_transcriptome_output = os.path.join(output_dir, "full_transcriptome_leaders.fa")
    process_transcript_list(None, fasta_file, gtf_file, full_transcriptome_output)

if __name__ == "__main__":
    main()
