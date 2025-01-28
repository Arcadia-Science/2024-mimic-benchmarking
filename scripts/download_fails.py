import json
import os
import re
import sys

import pandas as pd
import requests

# Input and output paths
input_dir = sys.argv[1]  # Directory containing PDB files
summary_file = sys.argv[2]  # Path to the summary file
metadata_dir = sys.argv[3]  # Directory containing JSON metadata files
fails_file = sys.argv[4]  # Path to downloadedviro3d_pdbs_fails.txt

# Base URL for downloading the PDB files
base_pdb_url = "https://viro3d.cvr.gla.ac.uk/api/pdb/"

# Load the summary file as a DataFrame
if os.path.exists(summary_file):
    summary_data = pd.read_csv(summary_file, sep="\t")
    summary_data.columns = summary_data.columns.str.strip()
else:
    # Create an empty DataFrame with required columns if the summary file doesn't exist
    summary_data = pd.DataFrame(
        columns=[
            "Virus Name",
            "Record ID",
            "ESMFold pLDDT",
            "ColabFold pLDDT",
            "Chosen Method",
            "structure_file",
        ]
    )


# Function to extract metadata from a filename
def extract_metadata(filename):
    match = re.match(r"(.*?)_(CF-|EF-)([A-Z0-9._]+)_relaxed\.pdb", filename)
    if match:
        virus_name = match.group(1).replace("_", " ")
        record_id = match.group(3)
        return virus_name, record_id
    return None, None


# Function to search all JSONs for metadata
def find_metadata_in_jsons(record_id, metadata_dir):
    """
    Searches all JSON files in the metadata directory for the given record_id
    and extracts the corresponding virus name, ESMFold pLDDT, and ColabFold pLDDT.
    """
    for json_file in os.listdir(metadata_dir):
        if json_file.endswith(".json"):
            json_path = os.path.join(metadata_dir, json_file)
            try:
                with open(json_path) as f:
                    metadata = json.load(f)

                    # Extract virus name
                    virus_name = metadata.get("virus_name", "").replace(" ", "_")

                    # Check if the record_id exists in the protein_structures list
                    protein_structures = metadata.get("protein_structures", [])
                    for protein in protein_structures:
                        if protein.get("record_id") == record_id:
                            esm_pLDDT = float(protein.get("esmfold_log_pLDDT", 0))
                            colab_pLDDT = float(protein.get("colabfold_json_pLDDT", 0))
                            return virus_name, esm_pLDDT, colab_pLDDT
            except json.JSONDecodeError:
                print(f"Error: Failed to decode JSON file {json_file}. Skipping...")
            except Exception as e:
                print(f"Unexpected error while processing {json_file}: {e}")

    # Return None if no match is found
    print(f"Record ID {record_id} not found in any JSON files.")
    return None, None, None


# Function to download a file and update the summary
def download_and_update(api_url, virus_name, record_id, chosen_method):
    filename = f"{virus_name.replace(' ', '_')}_{chosen_method[:2]}-{record_id}_relaxed.pdb"
    output_path = os.path.join(input_dir, filename)

    try:
        response = requests.get(api_url, stream=True)
        if response.status_code == 200:
            with open(output_path, "wb") as pdb_output:
                for chunk in response.iter_content(chunk_size=8192):
                    pdb_output.write(chunk)
            print(f"Successfully downloaded {output_path}")

            # Append to summary file
            if record_id not in summary_data["Record ID"].values:
                esm_pLDDT, colab_pLDDT = None, None
                # If JSON metadata exists, retrieve pLDDT scores
                json_filename = f"{virus_name.replace('_', ' ')}.json"
                json_path = os.path.join(metadata_dir, json_filename)
                if os.path.exists(json_path):
                    _, esm_pLDDT, colab_pLDDT = find_metadata_in_jsons(record_id, metadata_dir)

                # Add a new row to the summary file
                summary_data.loc[len(summary_data)] = {
                    "Virus Name": virus_name,
                    "Record ID": record_id,
                    "ESMFold pLDDT": esm_pLDDT,
                    "ColabFold pLDDT": colab_pLDDT,
                    "Chosen Method": ("ColabFold" if chosen_method == "CF" else "ESMFold"),
                    "structure_file": filename,
                }
                print(f"Added {record_id} to summary.")
            return True
        else:
            print(f"Failed to download {api_url}. HTTP Status: {response.status_code}")
            return False
    except Exception as e:
        print(f"Error downloading {api_url}: {e}")
        return False


# Step 1: Retry URLs in `fails_file`
if os.path.exists(fails_file):
    with open(fails_file) as f:
        failed_urls = [line.strip() for line in f if line.strip()]

    for url in failed_urls:
        filename = os.path.basename(url)
        virus_name, record_id = extract_metadata(filename)
        chosen_method = "CF" if "CF-" in filename else "EF"

        if not virus_name or not record_id:
            print(f"Could not extract metadata from {filename}. Searching JSONs...")
            record_id_match = re.search(r"(EF-|CF-)([A-Z0-9._]+)_relaxed", filename)
            if record_id_match:
                record_id = record_id_match.group(2)
                virus_name, esm_pLDDT, colab_pLDDT = find_metadata_in_jsons(record_id, metadata_dir)
                if not virus_name:
                    print(f"Record ID {record_id} not found in metadata. Skipping.")
                    continue
        print(f"Retrying download for {url}...")
        if not download_and_update(url, virus_name, record_id, chosen_method):
            # Try the opposite prefix if the download fails
            opposite_url = url.replace("CF-", "EF-") if "CF-" in url else url.replace("EF-", "CF-")
            opposite_method = "CF" if chosen_method == "EF" else "EF"
            print(f"Retrying with opposite prefix: {opposite_url}")
            download_and_update(opposite_url, virus_name, record_id, opposite_method)

# Step 2: Process 22-byte files
for pdb_file in os.listdir(input_dir):
    input_path = os.path.join(input_dir, pdb_file)

    if not os.path.isfile(input_path) or os.path.getsize(input_path) != 22:
        continue

    print(f"File {pdb_file} is 22 bytes. Attempting to download the opposite prefix...")

    virus_name, record_id = extract_metadata(pdb_file)
    if not virus_name or not record_id:
        print(f"Could not extract metadata from {pdb_file}. Searching JSONs...")
        record_id_match = re.search(r"(EF-|CF-)([A-Z0-9._]+)_relaxed", pdb_file)
        if record_id_match:
            record_id = record_id_match.group(2)
            virus_name, esm_pLDDT, colab_pLDDT = find_metadata_in_jsons(record_id, metadata_dir)
            if not virus_name:
                print(f"Record ID {record_id} not found in metadata. Skipping.")
                continue

    if "EF-" in pdb_file:
        api_filename = pdb_file.replace("EF-", "CF-")
        chosen_method = "CF"
    elif "CF-" in pdb_file:
        api_filename = pdb_file.replace("CF-", "EF-")
        chosen_method = "EF"
    else:
        print(f"File {pdb_file} does not contain a recognizable prefix. Skipping.")
        continue

    pdb_url = f"{base_pdb_url}{api_filename}"
    download_and_update(pdb_url, virus_name, record_id, chosen_method)

# Save the updated summary file
summary_data.to_csv(summary_file, sep="\t", index=False)
print(f"Updated summary file saved to {summary_file}")
