import argparse
import json
import os

import pandas as pd


def preprocess_viro3d_metadata(viro3d_metadata_dir):
    """
    Preprocess the Viro3D metadata to create a mapping of GenBank_ID -> genbank_name_curated.
    """
    genbank_mapping = {}
    for json_file in os.listdir(viro3d_metadata_dir):
        if json_file.endswith(".json"):
            json_path = os.path.join(viro3d_metadata_dir, json_file)
            with open(json_path) as f:
                data = json.load(f)
                for protein in data.get("protein_structures", []):
                    genbank_id = protein.get("genbank_id")
                    genbank_name_curated = protein.get(
                        "genbank_name_curated", "No curated name found"
                    )
                    if genbank_id and genbank_name_curated:
                        genbank_mapping[genbank_id] = genbank_name_curated
    return genbank_mapping


def find_refseq_pdb_file(refseq_id, viral_structures_dir):
    """
    Find the file name in the viral structures directory for a given RefSeq ID.
    """
    pattern = f"*{refseq_id}__*.pdb"
    for _, _, files in os.walk(viral_structures_dir):
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                return file
    return "No PDB Found"


def main(args):
    # Preprocess Viro3D metadata
    print("Preprocessing Viro3D metadata...")
    genbank_mapping = preprocess_viro3d_metadata(args.viro3d_metadata_dir)
    print(f"Loaded {len(genbank_mapping)} GenBank ID mappings.")

    # Load the random1000 RefSeq list
    refseq_ids = []
    with open(args.random1000_refseq_list) as f:
        refseq_ids = [line.strip() for line in f]

    # Load mapping file: RefSeq -> GenBank
    mapping = {}
    with open(args.mapping_file) as map_file:
        for line in map_file:
            cols = line.strip().split("\t")
            if len(cols) >= 2:
                refseq, genbank = cols[0], cols[1]
                if genbank != "No Match Found":
                    mapping[refseq] = genbank

    # Create the output table
    output_data = []
    for refseq_id in refseq_ids:
        genbank_id = mapping.get(refseq_id, "No Match")
        refseq_pdb_filename = find_refseq_pdb_file(refseq_id, args.viral_structures_dir)
        genbank_entry_name = (
            genbank_mapping.get(genbank_id, "No curated name found")
            if genbank_id != "No Match"
            else "No Match"
        )
        output_data.append([refseq_id, genbank_id, refseq_pdb_filename, genbank_entry_name])

    # Save the output table
    output_df = pd.DataFrame(
        output_data,
        columns=[
            "RefSeq_ID",
            "GenBank_ID",
            "RefSeq_PDB_Filename",
            "GenBank_Entry_Name",
        ],
    )
    output_df.to_csv(args.output_table, sep="\t", index=False)
    print(f"Random1000 table saved to {args.output_table}")


if __name__ == "__main__":
    import fnmatch

    parser = argparse.ArgumentParser(
        description="Create a Random1000 table with additional details."
    )
    parser.add_argument(
        "--random1000_refseq_list",
        required=True,
        help="Path to the random1000 refseq list.",
    )
    parser.add_argument(
        "--mapping_file", required=True, help="Path to RefSeq-to-GenBank mapping file."
    )
    parser.add_argument(
        "--viro3d_metadata_dir",
        required=True,
        help="Path to the Viro3D metadata directory.",
    )
    parser.add_argument(
        "--viral_structures_dir",
        required=True,
        help="Path to the viral structures directory.",
    )
    parser.add_argument("--output_table", required=True, help="Path to the output table.")
    args = parser.parse_args()
    main(args)
