import json
import os
import re
import sys

import pandas as pd

input_dir = sys.argv[1]  # Directory containing PDB files
summary_file = sys.argv[2]  # Path to the summary file
metadata_dir = sys.argv[3]  # Path to metadata directory (JSON)
fails_file = sys.argv[4]  # Path to structured failed downloads TSV file
new_summary_file = sys.argv[5]  # Path to save the updated summary file

print(f"✅ Input Directory: {input_dir}")
print(f"✅ Summary File: {summary_file}")
print(f"✅ Metadata Directory: {metadata_dir}")
print(f"✅ Fails File: {fails_file}")
print(f"✅ New Summary File: {new_summary_file}")

# Load existing summary data
if os.path.exists(summary_file):
    summary_data = pd.read_csv(summary_file, sep="\t")
else:
    print("❌ Summary file not found. Exiting.")
    sys.exit(1)


# Function to extract metadata (UniProt ID, pLDDT) from JSON files
def extract_metadata(record_id):
    for json_file in os.listdir(metadata_dir):
        if json_file.endswith(".json"):
            json_path = os.path.join(metadata_dir, json_file)
            with open(json_path) as f:
                try:
                    data = json.load(f)
                    protein_structures = data.get("protein_structures", [])
                    for protein in protein_structures:
                        if protein.get("record_id") == record_id:
                            return {
                                "Virus Name": protein.get("Virus_name_s_", "Unknown").replace(
                                    " ", "_"
                                ),
                                "Record ID": record_id,
                                "UniProt ID": protein.get("uniprot_id", "N/A"),
                                "ESMFold pLDDT": float(protein.get("esmfold_log_pLDDT", 0)),
                                "ColabFold pLDDT": float(protein.get("colabfold_json_pLDDT", 0)),
                            }
                except json.JSONDecodeError:
                    print(f"⚠️ Skipping corrupted JSON file: {json_file}")
    return None


# Function to retry failed downloads
def download_and_update(api_url, virus_name, record_id, old_method):
    global summary_data

    # Extract prefix (CF- or EF-)
    match = re.search(r"(CF-|EF-)", api_url)
    original_prefix = match.group(0) if match else ""

    # Determine the opposite prefix
    opposite_prefix = "CF-" if original_prefix == "EF-" else "EF-"

    # Construct opposite URL
    opposite_api_url = api_url.replace(original_prefix, opposite_prefix)

    # Define filenames
    original_filename = f"{virus_name.replace(' ', '_')}_{original_prefix}{record_id}_relaxed.pdb"
    opposite_filename = f"{virus_name.replace(' ', '_')}_{opposite_prefix}{record_id}_relaxed.pdb"

    original_output_path = os.path.join(input_dir, original_filename)
    opposite_output_path = os.path.join(input_dir, opposite_filename)

    # Step 1: Retry the original failed URL
    print(f"🔄 Retrying failed URL: {api_url}")
    curl_command = f'curl -s -L -A "Mozilla/5.0" -o "{original_output_path}" "{api_url}"'
    os.system(curl_command)

    # Check if the original URL was successful
    if os.path.exists(original_output_path) and os.path.getsize(original_output_path) > 22:
        print(f"✅ Successfully downloaded {original_output_path}")
        chosen_method = old_method
        structure_file = original_filename
    else:
        print(f"❌ Original retry failed. Trying opposite prefix: {opposite_api_url}")

        # Step 2: If first retry fails, try opposite prefix
        curl_command = (
            f'curl -s -L -A "Mozilla/5.0" -o "{opposite_output_path}" "{opposite_api_url}"'
        )
        os.system(curl_command)

        # Check if the opposite file was successful
        if os.path.exists(opposite_output_path) and os.path.getsize(opposite_output_path) > 22:
            print(f"✅ Successfully downloaded {opposite_output_path}")
            chosen_method = "ColabFold" if old_method == "ESMFold" else "ESMFold"
            structure_file = opposite_filename
        else:
            print(f"❌ Both retry attempts failed for {record_id}. Removing corrupted files.")
            remove_corrupted_files([original_output_path, opposite_output_path])
            return False

    # If download was successful, update the summary
    metadata = extract_metadata(record_id) or {
        "Virus Name": virus_name.replace(" ", "_"),
        "Record ID": record_id,
        "UniProt ID": "N/A",
        "ESMFold pLDDT": None,
        "ColabFold pLDDT": None,
    }

    if record_id in summary_data["Record ID"].values:
        summary_data.loc[summary_data["Record ID"] == record_id, "Chosen Method"] = chosen_method
        summary_data.loc[summary_data["Record ID"] == record_id, "structure_file"] = structure_file
        print(f"🔄 Updated summary: {record_id} now uses {chosen_method} and {structure_file}")
    else:
        print(f"🆕 Adding missing record {record_id} with {chosen_method}")
        summary_data = pd.concat(
            [
                summary_data,
                pd.DataFrame(
                    [{**metadata, "Chosen Method": chosen_method, "structure_file": structure_file}]
                ),
            ],
            ignore_index=True,
        )

    return True


# Function to remove 22-byte corrupted files
def remove_corrupted_files(files):
    for file in files:
        if os.path.exists(file) and os.path.getsize(file) == 22:
            print(f"🗑️ Removing corrupted 22-byte file: {file}")
            os.remove(file)


# Process all failed downloads
if os.path.exists(fails_file):
    failed_df = pd.read_csv(fails_file, sep="\t")

    for _, row in failed_df.iterrows():
        virus_name = row["Virus Name"]
        record_id = row["Record ID"]
        chosen_method = row["Chosen Method"]
        failed_url = row["Failed URL"]

        download_and_update(failed_url, virus_name, record_id, chosen_method)


# Function to save the updated summary file
def save_updated_summary():
    global summary_data
    summary_data = summary_data.drop_duplicates(subset=["Virus Name", "Record ID"], keep="first")
    summary_data.to_csv(new_summary_file, sep="\t", index=False)
    print(f"✅ Updated summary file saved to {new_summary_file}")


# Final cleanup: Remove any remaining 22-byte files
for pdb_file in os.listdir(input_dir):
    file_path = os.path.join(input_dir, pdb_file)
    if os.path.isfile(file_path) and os.path.getsize(file_path) == 22:
        print(f"🗑️ Removing corrupted 22-byte file: {file_path}")
        os.remove(file_path)

# Run all necessary steps
try:
    save_updated_summary()
except Exception as e:
    print(f"❌ Error while updating summary file: {e}")
    sys.exit(1)

# Ensure the output file exists before the script exits
if not os.path.exists(new_summary_file):
    print(f"❌ Critical error: {new_summary_file} was not created!")
    sys.exit(1)
