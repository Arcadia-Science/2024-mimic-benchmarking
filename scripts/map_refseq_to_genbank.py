import os
import json
import pandas as pd
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Map RefSeq IDs to GenBank IDs using Viro3D metadata.")
parser.add_argument("--metadata", required=True, help="Path to viral_structure_metadata.tsv")
parser.add_argument("--viro3d_metadata", required=True, help="Path to directory containing Viro3D JSON metadata files")
parser.add_argument("--output", required=True, help="Output path for the mapping file")
args = parser.parse_args()

# Load metadata
metadata_tsv = args.metadata
viro3d_metadata_dir = args.viro3d_metadata
output_mapping_file = args.output

metadata = pd.read_csv(metadata_tsv, sep="\t")

# Ensure the 'nomburg_protein_name' column exists
if "nomburg_protein_name" not in metadata.columns:
    raise ValueError("The column 'nomburg_protein_name' is missing from the metadata file.")

nomburg_protein_names = metadata["nomburg_protein_name"]
nomburg_entries = [
    {
        "product_name": entry.split("__")[0],
        "refseq_id": entry.split("__")[1],
        "virus_name": entry.split("__")[2],
        "length": metadata.loc[metadata["nomburg_protein_name"] == entry, "length"].values[0]
    }
    for entry in nomburg_protein_names if isinstance(entry, str) and len(entry.split("__")) == 4
]

# Initialize mapping results
mapping_results = []

# Process each Nomburg entry
for entry in nomburg_entries:
    refseq_id = entry["refseq_id"]
    product_name = entry["product_name"]
    virus_name_safe = entry["virus_name"].replace(" ", "_")
    nomburg_length = entry["length"]

    # Path to the corresponding JSON file
    json_file_path = os.path.join(viro3d_metadata_dir, f"{virus_name_safe}.json")

    if not os.path.exists(json_file_path):
        print(f"[DEBUG]: Metadata file for {virus_name_safe} not found.")
        continue

    # Load JSON metadata
    with open(json_file_path, "r") as f:
        viro3d_data = json.load(f)

    # Extract protein structures
    protein_structures = viro3d_data.get("protein_structures", [])
    matches = []
    match_priority = None  # Track the highest-priority match type

    for protein in protein_structures:
        genbank_id = protein.get("genbank_id", "")
        protlen = protein.get("protlen")

        # Match by 'product' field
        product_field = protein.get("product", "").lower()
        if product_field and product_name.lower() == product_field:
            matches.append((genbank_id, "Product", protlen))
            match_priority = "Product"

        # Match by 'genbank_name_curated' field
        curated_name = protein.get("genbank_name_curated", "").lower()
        if "product:" in curated_name:
            extracted_product = curated_name.split("product:")[1].split(";")[0].strip()
        else:
            extracted_product = ""

        if extracted_product and product_name.lower().replace("_", " ") in extracted_product:
            if match_priority is None or match_priority != "Product":
                matches.append((genbank_id, "GenBank Name Curated (Extracted Product)", protlen))
                match_priority = "GenBank Name Curated (Extracted Product)"

        # Match by length (as a fallback)
        if not matches and protlen and nomburg_length and protlen == nomburg_length:
            matches.append((genbank_id, "Length", protlen))

    # Deduplicate matches and filter by length if multiple matches exist
    filtered_matches = []
    if matches:
        # Keep only matches with the highest-priority type
        filtered_by_priority = [match for match in matches if match[1] == match_priority]

        # If multiple matches remain, filter by length
        if len(filtered_by_priority) > 1:
            filtered_matches = [
                match for match in filtered_by_priority if match[2] == nomburg_length
            ]
        else:
            filtered_matches = filtered_by_priority

    # Save all matches for this RefSeq ID
    if filtered_matches:
        for match in filtered_matches:
            mapping_results.append({
                "RefSeq_ID": refseq_id,
                "GenBank_ID": match[0],
                "Match_Type": match[1]
            })
    else:
        mapping_results.append({
            "RefSeq_ID": refseq_id,
            "GenBank_ID": "No Match Found",
            "Match_Type": "No Match"
        })

# Write results to output file
with open(output_mapping_file, "w") as f:
    f.write("RefSeq_ID\tGenBank_ID\tMatch_Type\n")
    for result in mapping_results:
        f.write(f"{result['RefSeq_ID']}\t{result['GenBank_ID']}\t{result['Match_Type']}\n")

print(f"[INFO]: Mapping complete. Results saved to {output_mapping_file}")
