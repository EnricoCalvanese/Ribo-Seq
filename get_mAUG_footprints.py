import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
import os
import gffutils
import sqlite3
from typing import List, Dict, Tuple

def validate_gff_file(gff3_file: str) -> bool:
    """
    Validate that the GFF3 file exists and is readable.
    
    Args:
        gff3_file: Path to the GFF3 file
        
    Returns:
        bool: True if file is valid
        
    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file isn't readable
    """
    if not os.path.exists(gff3_file):
        raise FileNotFoundError(f"GFF3 file not found: {gff3_file}")
    
    if not os.access(gff3_file, os.R_OK):
        raise PermissionError(f"Cannot read GFF3 file: {gff3_file}")
    
    return True

def create_gff_db(gff3_file: str, force_rebuild: bool = False) -> gffutils.FeatureDB:
    """
    Create or connect to a GFF database with comprehensive error handling.
    
    Args:
        gff3_file: Path to the GFF3 file
        force_rebuild: If True, rebuild database even if it exists
        
    Returns:
        gffutils database object
        
    Raises:
        Various exceptions with descriptive error messages
    """
    # Validate input file
    validate_gff_file(gff3_file)
    
    db_path = gff3_file + '.db'
    
    # Remove existing database if force rebuild
    if force_rebuild and os.path.exists(db_path):
        os.remove(db_path)
    
    try:
        # First try to connect to existing database
        if os.path.exists(db_path):
            print(f"Connecting to existing GFF database: {db_path}")
            return gffutils.FeatureDB(db_path)
        
        print(f"Creating new GFF database: {db_path}")
        # Create new database with careful error handling
        db = gffutils.create_db(
            gff3_file,
            dbfn=db_path,
            force=True,  # Overwrite any existing partially created database
            merge_strategy='create_unique',
            id_spec={'gene': ['ID', 'Name'],
                    'mRNA': ['ID', 'transcript_id'],
                    'CDS': ['ID', 'Parent']},
            verbose=True,  # Show progress
            force_merge_fields=['ID'],
            keep_order=True,
            sort_attribute_values=True
        )
        
        # Verify database was created successfully
        test_db = gffutils.FeatureDB(db_path)
        # Try a simple query to validate the database
        next(test_db.all_features(), None)
        
        return test_db
        
    except sqlite3.OperationalError as e:
        raise RuntimeError(f"Database creation failed due to SQLite error: {str(e)}\n"
                         f"Try deleting the existing database file: {db_path}")
    
    except Exception as e:
        # If database creation fails, ensure partial database is cleaned up
        if os.path.exists(db_path):
            os.remove(db_path)
        raise RuntimeError(f"Failed to create GFF database: {str(e)}")

def main():
    # File paths
    BASE_DIR = "/global/scratch/users/enricocalvane/riboseq/imb2"
    GFF3_FILE = "/global/scratch/users/enricocalvane/riboseq/Xu2017/tair10_reference/Arabidopsis_thaliana.TAIR10.60.gff3"
    
    try:
        print("Creating GFF database...")
        # Try creating database with force rebuild
        db = create_gff_db(GFF3_FILE, force_rebuild=True)
        print("GFF database created successfully!")
        
        # Verify database functionality
        print("Validating database...")
        feature_count = sum(1 for _ in db.all_features())
        print(f"Database contains {feature_count} features")
        
        # Continue with rest of script...
        
    except Exception as e:
        print(f"Error: {str(e)}")
        print("\nTroubleshooting steps:")
        print("1. Check GFF3 file permissions and format")
        print("2. Delete any existing .db files")
        print("3. Ensure enough disk space")
        print("4. Check for special characters in file path")
        raise

if __name__ == "__main__":
    main()
