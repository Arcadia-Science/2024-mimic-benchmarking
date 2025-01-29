import json
import os
import sys

import pandas as pd

response_dir = sys.argv[1]
output_dir = sys.argv[2]
summary_path = sys.argv[3]  # Path to save the summary table
failed_downloads_path = sys.argv[4]  # Path to save the list of failed downloads

base_pdb_url = "https://viro3d.cvr.gla.ac.uk/api/pdb/"

os.makedirs(output_dir, exist_ok=True)

# Initialize a list to store summary information
summary_data = []

# Initialize a list to store failed downloads
failed_downloads = []

# Loop through JSON response files
for response_file in os.listdir(response_dir):
    if response_file.endswith(".json"):
        response_path = os.path.join(response_dir, response_file)
        with open(response_path) as f:
            try:
                data = json.load(f)
                virus_name = data.get("virus_name", "unknown_virus").replace(" ", "_")
                protein_structures = data.get("protein_structures", [])

                for protein in protein_structures:
                    # Extract record-level details
                    record_id = protein.get("record_id")
                    esm_pLDDT = float(protein.get("esmfold_log_pLDDT", 0))
                    colab_pLDDT = float(protein.get("colabfold_json_pLDDT", 0))

                    # Decide prefix based on the higher pLDDT score
                    if colab_pLDDT > esm_pLDDT:
                        prefix = "CF-"
                        chosen_method = "ColabFold"
                    else:
                        prefix = "EF-"
                        chosen_method = "ESMFold"

                    if record_id:
                        pdb_url = f"{base_pdb_url}{prefix}{record_id}_relaxed.pdb"
                        pdb_output_path = os.path.join(
                            output_dir, f"{virus_name}_{prefix}{record_id}_relaxed.pdb"
                        )

                        print(f"Downloading {pdb_url} to {pdb_output_path}...")
                        # Download the file
                        exit_code = os.system(f"curl -s -o '{pdb_output_path}' {pdb_url}")

                        if (
                            exit_code == 0
                            and os.path.exists(pdb_output_path)
                            and os.path.getsize(pdb_output_path) > 22
                        ):
                            # Log successful download
                            print(f"Successfully downloaded {pdb_url}")
                            # Append details to the summary data
                            summary_data.append(
                                {
                                    "Virus Name": virus_name,
                                    "Record ID": record_id,
                                    "ESMFold pLDDT": esm_pLDDT,
                                    "ColabFold pLDDT": colab_pLDDT,
                                    "Chosen Method": chosen_method,
                                }
                            )
                        else:
                            # Log failed download
                            print(f"Failed to download {pdb_url}")
                            failed_downloads.append(pdb_url)
            except json.JSONDecodeError:
                print(f"Failed to decode JSON from {response_file}")
            except ValueError:
                print(f"Invalid pLDDT value in {response_file}")

# Create a DataFrame from the summary data
summary_df = pd.DataFrame(summary_data)

# Save the summary table to the specified path
summary_df.to_csv(summary_path, sep="\t", index=False)
print(f"Summary table saved to {summary_path}")

# Save the list of failed downloads to a file
with open(failed_downloads_path, "w") as failed_file:
    for url in failed_downloads:
        failed_file.write(f"{url}\n")
print(f"Failed downloads saved to {failed_downloads_path}")
