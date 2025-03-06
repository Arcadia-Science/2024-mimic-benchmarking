import os
import sys

import pandas as pd

# Define file paths
tsv_file = sys.argv[1]
pdb_folder = sys.argv[2]

# Load the TSV file
df = pd.read_csv(tsv_file, sep="\t")

# Count unique Record ID entries
unique_record_ids = df["Record ID"].nunique()
total_record_ids = df["Record ID"].count()
print(f"📝 Total 'Record ID' entries in TSV: {total_record_ids}")
print(f"📝 Unique 'Record ID' entries in TSV: {unique_record_ids}")

# Count .pdb files in the folder
pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith(".pdb")]
num_pdb_files = len(pdb_files)
print(f"📁 Number of .pdb files in folder: {num_pdb_files}")

# Compare counts
if unique_record_ids == num_pdb_files:
    print("✅ The total number of unique Record ID entries matches the number of .pdb files.")
else:
    print("⚠️ Mismatch detected!")
    print(f"🔍 Difference: {abs(unique_record_ids - num_pdb_files)}")
