import os
import sys
import json

# Inputs
response_dir = sys.argv[1]
output_dir = sys.argv[2]

# Base URL for PDB download
base_pdb_url = "https://viro3d.cvr.gla.ac.uk/api/pdb/"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop through JSON response files
for response_file in os.listdir(response_dir):
    if response_file.endswith(".json"):
        response_path = os.path.join(response_dir, response_file)
        with open(response_path, 'r') as f:
            try:
                data = json.load(f)
                virus_name = data.get("virus_name", "unknown_virus").replace(" ", "_")
                protein_structures = data.get("protein_structures", [])
                for protein in protein_structures:
                    record_id = protein.get("record_id")
                    if record_id:
                        # Construct PDB URL
                        pdb_url = f"{base_pdb_url}{record_id}_relaxed.pdb"
                        # Construct output file path
                        pdb_output_path = os.path.join(output_dir, f"{virus_name}_{record_id}_relaxed.pdb")
                        # Download the PDB file
                        print(f"Downloading {pdb_url} to {pdb_output_path}...")
                        os.system(f"curl -s -o {pdb_output_path} {pdb_url}")
            except json.JSONDecodeError:
                print(f"Failed to decode JSON from {response_file}")
