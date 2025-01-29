import re
import sys

import pandas as pd

# Input and output file paths
vmr_metadata_path = sys.argv[1]
virushostdb_path = sys.argv[2]
output_path = sys.argv[3]


def clean_and_split_ids(column, regex=r"NC_\d{6}"):
    """
    Extracts valid REFSEQ IDs (e.g., NC_123456) from a column with concatenated text.
    """
    column = column.astype(str)  # Ensure column is string
    return column.apply(lambda x: re.findall(regex, x))  # Find all matching IDs


# Load the data
vmr_metadata = pd.read_csv(vmr_metadata_path, sep="\t")
virushostdb = pd.read_csv(virushostdb_path, sep="\t")

# Clean and split REFSEQ IDs in both datasets
vmr_metadata["Virus REFSEQ accession"] = clean_and_split_ids(vmr_metadata["Virus REFSEQ accession"])
virushostdb["refseq id"] = clean_and_split_ids(virushostdb["refseq id"])

# Drop rows with empty or NaN `Virus REFSEQ accession` in vmr_metadata
vmr_metadata = vmr_metadata[vmr_metadata["Virus REFSEQ accession"].str.len() > 0].reset_index(
    drop=True
)

print(f"Filtered VMR Metadata: {len(vmr_metadata)} rows remain after removing blanks")

# Explode both datasets to normalize the REFSEQ IDs (one row per ID)
vmr_metadata_exploded = vmr_metadata.explode("Virus REFSEQ accession")
virushostdb_exploded = virushostdb.explode("refseq id")

print(f"Exploded VMR Metadata: {len(vmr_metadata_exploded)} rows")
print(f"Exploded VirusHostDB: {len(virushostdb_exploded)} rows")

# Perform the merge on cleaned REFSEQ IDs
merged_data = pd.merge(
    vmr_metadata_exploded,
    virushostdb_exploded,
    left_on="Virus REFSEQ accession",
    right_on="refseq id",
    how="left",
)

print(f"Merged Data: {len(merged_data)} rows")

# Debugging: Check for rows with missing `host tax id`
unmatched_rows = merged_data[merged_data["host tax id"].isna()]
print(f"Unmatched Rows: {len(unmatched_rows)} rows (These have no matching host tax id)")

# Filter rows where host tax id is 9606 (humans)
if "host tax id" in merged_data.columns:
    filtered_data = merged_data[merged_data["host tax id"] == 9606].copy()
    print(f"Filtered Data (host tax id == 9606): {len(filtered_data)} rows")
else:
    print("Warning: 'host tax id' column not found in merged data.")
    filtered_data = pd.DataFrame(columns=merged_data.columns)

# Always include these specific virus names, regardless of host tax id
virus_names_to_keep = [
    "Thogoto virus",  # https://doi.org/10.3201/eid2802.211270
    "Piry virus",  # https://doi.org/10.1007/bf01241673
    "Vesicular stomatitis New Jersey virus",  # https://doi.org/10.1056/nejm196711092771901
]

# Filter rows matching these virus names
extra_rows = merged_data[merged_data["Virus name(s)"].isin(virus_names_to_keep)].copy()
print(f"Extra Rows (specific virus names): {len(extra_rows)} rows")

# Combine filtered rows and extra rows
final_data = pd.concat([filtered_data, extra_rows]).drop_duplicates()

# Deduplicate the merged data to keep only one row per unique `Virus name(s)`
deduplicated_data = final_data.drop_duplicates(subset=["Virus name(s)"])
print(f"Deduplicated Data: {len(deduplicated_data)} rows")

# Save the deduplicated data
deduplicated_data.to_csv(output_path, sep="\t", index=False)
print(f"Saved deduplicated data to {output_path}")
