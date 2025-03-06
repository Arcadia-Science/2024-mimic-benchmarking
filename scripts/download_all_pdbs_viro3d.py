import json
import os
import sys

import pandas as pd

# Input arguments
response_dir = sys.argv[1]  # Directory with JSON response files
output_dir = sys.argv[2]  # Directory to save downloaded PDB files
summary_path = sys.argv[3]  # Path to save the summary table
failed_downloads_path = sys.argv[4]  # Path to save failed downloads (TSV file)

# Base URL for downloading PDB files
base_pdb_url = "https://viro3d.cvr.gla.ac.uk/api/pdb/"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Initialize storage
summary_data = []
failed_downloads = []

# DEBUG: Print the directory contents
print(f"📂 Checking response directory: {response_dir}")
print(os.listdir(response_dir))

# Loop through JSON response files
for response_file in os.listdir(response_dir):
    if response_file.endswith(".json"):
        response_path = os.path.join(response_dir, response_file)
        print(f"📄 Processing: {response_file}")

        with open(response_path) as f:
            try:
                data = json.load(f)

                # Handle cases where JSON is a dict instead of list
                if isinstance(data, dict):
                    protein_structures = data.get("protein_structures", [])
                elif isinstance(data, list):
                    protein_structures = data
                else:
                    print(f"❌ Unexpected JSON format in {response_file}")
                    continue

                # DEBUG: Check if JSON file contains proteins
                print(f"🔍 Found {len(protein_structures)} protein structures in {response_file}")

                for protein in protein_structures:
                    record_id = protein.get("record_id")
                    virus_name = (
                        protein.get("Virus_name_s_", "unknown_virus").lower().replace(" ", "_")
                    )
                    uniprot_id = protein.get("uniprot_id", "N/A")

                    print(
                        f"🧐 Extracted: record_id={record_id}, "
                        f"virus_name={virus_name}, uniprot_id={uniprot_id}"
                    )

                    if not record_id:
                        print(f"⚠️ Skipping entry with missing record_id: {protein}")
                        continue

                    esm_pLDDT = float(protein.get("esmfold_log_pLDDT", 0))
                    colab_pLDDT = float(protein.get("colabfold_json_pLDDT", 0))

                    # Determine prefix based on pLDDT score
                    if colab_pLDDT > esm_pLDDT:
                        prefix = "CF-"
                        chosen_method = "ColabFold"
                    else:
                        prefix = "EF-"
                        chosen_method = "ESMFold"

                    pdb_url = f"{base_pdb_url}{prefix}{record_id}_relaxed.pdb"
                    pdb_output_path = os.path.join(
                        output_dir, f"{virus_name}_{prefix}{record_id}_relaxed.pdb"
                    )

                    print(f"🚀 Attempting to download: {pdb_url}")

                    # Run curl and capture errors
                    exit_code = os.system(
                        f"curl -s -o '{pdb_output_path}' {pdb_url} || echo '❌ Download failed'"
                    )

                    if exit_code == 0 and os.path.exists(pdb_output_path):
                        file_size = os.path.getsize(pdb_output_path)

                        if file_size > 22:
                            print(f"✅ Successfully downloaded {pdb_url}")
                            summary_data.append(
                                {
                                    "Virus Name": virus_name,
                                    "Record ID": record_id,
                                    "UniProt ID": uniprot_id,
                                    "ESMFold pLDDT": esm_pLDDT,
                                    "ColabFold pLDDT": colab_pLDDT,
                                    "Chosen Method": chosen_method,
                                }
                            )
                        else:
                            print(f"⚠️ Corrupted 22-byte file detected: {pdb_output_path}")
                            os.remove(pdb_output_path)

                            # Log the failure
                            failed_downloads.append(
                                [
                                    virus_name,
                                    record_id,
                                    chosen_method,
                                    pdb_url,
                                    pdb_output_path,
                                    "Corrupted (22 bytes)",
                                ]
                            )
                    else:
                        print(f"❌ Failed to download {pdb_url}")
                        failed_downloads.append(
                            [
                                virus_name,
                                record_id,
                                chosen_method,
                                pdb_url,
                                pdb_output_path,
                                "Download failed",
                            ]
                        )

            except json.JSONDecodeError:
                print(f"❌ JSON decoding failed for {response_file}")
            except ValueError:
                print(f"⚠️ Invalid pLDDT value in {response_file}")

# DEBUG: Check if we processed anything
if not summary_data:
    print("⚠️ WARNING: No PDBs were successfully downloaded!")

# Save the summary table
summary_df = pd.DataFrame(summary_data)
summary_df = pd.DataFrame(summary_data).drop_duplicates(keep="first")
summary_df.to_csv(summary_path, sep="\t", index=False)
print(f"✅ Summary table saved to {summary_path}")

# Save failed downloads in a structured TSV file
failed_df = pd.DataFrame(
    failed_downloads,
    columns=[
        "Virus Name",
        "Record ID",
        "Chosen Method",
        "Failed URL",
        "Intended Filename",
        "Failure Reason",
    ],
)
failed_df.to_csv(failed_downloads_path, sep="\t", index=False)
print(f"❌ Failed downloads saved to {failed_downloads_path}")
