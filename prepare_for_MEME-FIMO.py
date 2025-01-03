import gffutils
import pyfaidx
from Bio import SeqIO
from Bio.Seq import Seq
import os
import sys
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('sequence_prep.log'),
        logging.StreamHandler()
    ]
)

def get_transcript_ids_for_gene(gene_id, db):
    """Get all transcript IDs associated with a gene ID."""
    base_gene_id = gene_id.split('.')[0]
    transcript_ids = []
    try:
        # First try exact match
        for feature in db.children(gene_id, featuretype='transcript'):
            transcript_ids.append(feature.id)
        
        # If no exact match, try base gene ID
        if not transcript_ids:
            for feature in db.features_of_type('transcript', order_by='start'):
                if feature.id.startswith(base_gene_id):
                    transcript_ids.append(feature.id)
    except Exception as e:
        logging.error(f"Error finding transcripts for gene {gene_id}: {str(e)}")
    
    return transcript_ids

def get_five_prime_leader(transcript, fasta, db):
    """Extract 5' leader sequence for a given transcript."""
    try:
        t = db[transcript]
        cds_list = list(db.children(transcript, featuretype='CDS', order_by='start'))
        if not cds_list:
            logging.warning(f"No CDS found for transcript {transcript}")
            return None
            
        if t.strand == '+':
            leader_start = t.start
            leader_end = cds_list[0].start - 1
        else:
            leader_start = cds_list[-1].end + 1
            leader_end = t.end
        
        if leader_end <= leader_start:
            logging.warning(f"Invalid leader coordinates for transcript {transcript}")
            return None
            
        seq = fasta[t.seqid][leader_start-1:leader_end].seq
        if t.strand == '-':
            seq = str(Seq(seq).reverse_complement())
        
        if not seq or len(seq) == 0:
            logging.warning(f"Empty sequence for transcript {transcript}")
            return None
            
        return seq
    except Exception as e:
        logging.error(f"Error processing transcript {transcript}: {str(e)}")
        return None

def process_single_category(input_file, fasta_file, gtf_file, output_file, db_file="arabidopsis.db"):
    """Process a single category of sequences."""
    logging.info(f"Starting processing of {input_file if input_file else 'full transcriptome'}")
    
    try:
        # Create or load database
        if not os.path.exists(db_file):
            logging.info(f"Creating database from {gtf_file}...")
            db = gffutils.create_db(gtf_file, db_file, force=True,
                                  merge_strategy='merge',
                                  id_spec={'transcript': ['transcript_id']},
                                  gtf_transcript_key='transcript_id')
        db = gffutils.FeatureDB(db_file)
        
        # Load genome
        fasta = pyfaidx.Fasta(fasta_file)
        
        # Get transcript list
        if input_file is None:
            transcripts = [f.id for f in db.features_of_type('transcript')]
        else:
            with open(input_file) as f:
                input_ids = [line.strip() for line in f]
                transcripts = []
                for input_id in input_ids:
                    trans_ids = get_transcript_ids_for_gene(input_id, db)
                    if trans_ids:
                        transcripts.extend(trans_ids)
                    else:
                        logging.warning(f"No transcripts found for {input_id}")
        
        # Process sequences
        processed = 0
        with open(output_file, 'w') as out:
            for transcript in transcripts:
                seq = get_five_prime_leader(transcript, fasta, db)
                if seq and len(seq) > 0:
                    out.write(f">{transcript}\n{seq}\n")
                    out.flush()  # Ensure immediate writing to file
                    processed += 1
                if processed % 100 == 0 and processed > 0:
                    logging.info(f"Processed {processed} transcripts for {output_file}")
        
        logging.info(f"Completed processing {processed} transcripts for {output_file}")
        
    except Exception as e:
        logging.error(f"Error in process_single_category: {str(e)}")
        raise

def main():
    if len(sys.argv) < 2:
        print("Usage: script.py <category>")
        print("Categories: uorf_down, uorf_up, uorf_nc, translatome_down, translatome_up, translatome_nc, full_transcriptome")
        sys.exit(1)
    
    category = sys.argv[1]
    
    # Base paths
    base_dir = "/global/scratch/users/enricocalvane/riboseq/imb2"
    ref_dir = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference"
    
    # Reference files
    gtf_file = os.path.join(ref_dir, "Arabidopsis_thaliana.TAIR10.60.gtf")
    fasta_file = os.path.join(ref_dir, "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
    
    # Create output directory
    output_dir = os.path.join(base_dir, "unique_reads/leader_sequences")
    os.makedirs(output_dir, exist_ok=True)
    
    # Category to filename mapping
    category_map = {
        "uorf_down": (os.path.join("systemPipeR/uORF_counts/AllGenes/active_uORF_transcripts_TE-down.txt"), "uORF_down_leaders.fa"),
        "uorf_up": (os.path.join("systemPipeR/uORF_counts/AllGenes/active_uORF_transcripts_TE-up.txt"), "uORF_up_leaders.fa"),
        "uorf_nc": (os.path.join("systemPipeR/uORF_counts/AllGenes/active_uORF_transcripts_TE-nc.txt"), "uORF_nc_leaders.fa"),
        "translatome_down": ("unique_reads/TEdown_genes.txt", "translatome_down_leaders.fa"),
        "translatome_up": ("unique_reads/TEup_genes.txt", "translatome_up_leaders.fa"),
        "translatome_nc": ("unique_reads/TEnc_genes.txt", "translatome_nc_leaders.fa"),
        "full_transcriptome": (None, "full_transcriptome_leaders.fa")
    }
    
    if category not in category_map:
        print(f"Invalid category: {category}")
        sys.exit(1)
    
    input_file, output_file = category_map[category]
    input_path = os.path.join(base_dir, input_file) if input_file else None
    output_path = os.path.join(output_dir, output_file)
    
    process_single_category(input_path, fasta_file, gtf_file, output_path)

if __name__ == "__main__":
    main()
