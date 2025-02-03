import argparse
import json
import os

import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge JSON metadata with summary file based on Record ID"
    )
    parser.add_argument(
        "--summary-file",
        required=True,
        help="Input summary file with Record IDs",
    )
    parser.add_argument(
        "--json-dir",
        required=True,
        help="Directory containing JSON metadata files",
    )
    parser.add_argument(
        "--output-file",
        required=True,
        help="Output path for merged metadata file",
    )
    return parser.parse_args()


def find_metadata_in_json(record_id, json_dir):
    """
    Searches JSON files for the given record_id and extracts metadata fields.
    """
    metadata_fields = [
        "protlen",
        "Genome_composition",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Virus_names",
    ]

    for json_file in os.listdir(json_dir):
        if json_file.endswith(".json"):
            json_path = os.path.join(json_dir, json_file)
            try:
                with open(json_path) as f:
                    metadata = json.load(f)

                # Check if this JSON file contains the record_id
                for protein in metadata.get("protein_structures", []):
                    if protein.get("record_id") == record_id:
                        return {field: protein.get(field, None) for field in metadata_fields}
            except json.JSONDecodeError:
                print(f"Error: Could not decode {json_file}. Skipping.")
            except Exception as e:
                print(f"Unexpected error processing {json_file}: {e}")

    return {field: None for field in metadata_fields}


def main():
    args = parse_args()

    # Load summary file
    print(f"Loading summary file: {args.summary_file}")
    summary_data = pd.read_csv(args.summary_file, sep="\t")
    print(f"Loaded {len(summary_data)} rows from summary file.")

    # Initialize new columns
    metadata_fields = [
        "protlen",
        "Genome_composition",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Virus_names",
    ]
    for field in metadata_fields:
        summary_data[field] = None

    # Count JSON files
    json_files = [f for f in os.listdir(args.json_dir) if f.endswith(".json")]
    print(f"Found {len(json_files)} JSON files in metadata directory.")

    # Process records
    records_found = 0
    records_not_found = 0

    for index, row in tqdm(summary_data.iterrows(), total=len(summary_data)):
        record_id = row["Record ID"]
        metadata_values = find_metadata_in_json(record_id, args.json_dir)

        if any(value is not None for value in metadata_values.values()):
            records_found += 1
        else:
            records_not_found += 1
            print(f"\nRecord ID not found: {record_id}")

        # Update metadata columns
        for field in metadata_fields:
            summary_data.at[index, field] = metadata_values[field]

    # Save results
    print(f"\nSaving updated summary to {args.output_file}")
    summary_data.to_csv(args.output_file, sep="\t", index=False)

    print("\nProcessing complete!")
    print(f"Total records processed: {len(summary_data)}")
    print(f"Records with metadata found: {records_found}")
    print(f"Records not found: {records_not_found}")


if __name__ == "__main__":
    main()
