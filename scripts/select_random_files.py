# scripts/select_random_files.py
import argparse
import glob
import random
from pathlib import Path
import re
import pandas as pd
import os
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description='Select random files from a directory')
    parser.add_argument('--input-folder', required=True, help='Input file folder')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--sizes', nargs='+', type=int, required=True, 
                      help='List of sample sizes')
    return parser.parse_args()

def parse_filename(filename):
    """
    Parse filename to extract virus name and record ID.
    """
    # Get base filename without extension
    base = Path(filename).stem
    
    # Split by CF- or EF-
    parts = re.split(r'[CE]F-', base)
    
    if len(parts) != 2:
        return None, None
    
    virus_part = parts[0].rstrip('_')  # Remove trailing underscore if present
    remaining = parts[1]
    
    # Convert underscores to spaces in virus name
    virus_name = virus_part.replace('_', ' ')
    
    # Get record ID (everything before the next underscore)
    record_id = remaining.split('_')[0]
    
    return virus_name, record_id

def create_summary_table(files, output_dir, host_organism):
    """
    Create summary table with parsed filename information.
    """
    records = []
    for file in files:
        virus_name, record_id = parse_filename(file)
        if virus_name and record_id:  # Only add if parsing succeeded
            records.append({
                'pdb': os.path.basename(file),  # Just the filename, not full path
                'virus': virus_name,
                'record_id': record_id
            })
    
    # Create DataFrame and save to TSV
    df = pd.DataFrame(records)
    summary_file = output_dir / f"{host_organism}_summary_1000.tsv"
    df.to_csv(summary_file, sep='\t', index=False)
    print(f"Created summary table with {len(df)} entries in {summary_file}")

def select_nested_random_files(files, sizes, seed=42):
    """
    Select random files ensuring larger sets contain smaller sets.
    """
    random.seed(seed)
    
    # Sort sizes to ensure proper nesting
    sizes = sorted(sizes)
    
    # Initialize selections dictionary
    selections = {}
    current_files = []
    remaining_files = files.copy()
    
    for size in sizes:
        # Calculate how many additional files we need
        additional_needed = size - len(current_files)
        
        if additional_needed > 0:
            # Select additional files from remaining pool
            if len(remaining_files) < additional_needed:
                raise ValueError(
                    f"Not enough files to select {additional_needed} more files. "
                    f"Only {len(remaining_files)} files remaining."
                )
            
            new_files = random.sample(remaining_files, additional_needed)
            
            # Remove selected files from remaining pool
            for file in new_files:
                remaining_files.remove(file)
            
            # Add new files to current selection
            current_files.extend(new_files)
        
        # Store current selection for this size
        selections[size] = current_files.copy()
    
    return selections

def main():
    args = parse_args()
    
    # Get list of all input files
    all_files = glob.glob(args.input_folder)
    if not all_files:
        raise ValueError(f"No files found matching pattern: {args.input_folder}")
    
    try:
        # Get output directory and host organism from input folder
        input_path = Path(args.input_folder)
        output_dir = input_path.parent.parent
        host_organism = output_dir.name
        
        # Get nested selections
        selections = select_nested_random_files(
            files=all_files,
            sizes=args.sizes,
            seed=args.seed
        )
        
        # Create directories and copy files for each selection
        for size, selected_files in selections.items():
            # Create directory for this size
            size_dir = output_dir / f"selected_files_{host_organism}_{size}"
            size_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy selected files to the new directory
            for file in selected_files:
                # Copy the file, preserving only the filename
                shutil.copy2(file, size_dir / Path(file).name)
            print(f"Created directory {size_dir} with {len(selected_files)} files")
            
            # Create summary table for 1000-file set
            if size == 1000:
                create_summary_table(selected_files, output_dir, host_organism)
            
            # Verify subset relationship
            if size > min(args.sizes):
                prev_size = max(s for s in args.sizes if s < size)
                prev_files = set(selections[prev_size])
                current_files = set(selected_files)
                assert prev_files.issubset(current_files), \
                    f"Set of size {size} does not contain all files from size {prev_size}"
    
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    main()
