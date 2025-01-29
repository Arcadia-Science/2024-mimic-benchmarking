import argparse
import os
import random
import shutil


def main(args):
    mapping = {}
    # Load mapping file: RefSeq -> GenBank
    with open(args.mapping_file) as map_file:
        for line in map_file:
            cols = line.strip().split("\t")
            if len(cols) >= 2:
                refseq, genbank = cols[0], cols[1]
                if genbank != "No Match Found" and genbank not in mapping.values():
                    mapping[refseq] = genbank

    # Reverse mapping to group RefSeq IDs by GenBank ID
    genbank_to_refseq = {}
    for refseq, genbank in mapping.items():
        if genbank not in genbank_to_refseq:
            genbank_to_refseq[genbank] = []
        genbank_to_refseq[genbank].append(refseq)

    # Initialize cumulative sets for inclusivity
    cumulative_genbank_ids = set()
    cumulative_refseq_ids = set()

    # Process each subset
    for i, refseq_output in enumerate(args.refseq_subsets):
        genbank_list = args.genbank_lists[i]
        genbank_subset_dir = args.genbank_subset_dirs[i]
        refseq_subset_dir = args.refseq_subset_dirs[i]
        random_size = args.random_sizes[i]

        os.makedirs(genbank_subset_dir, exist_ok=True)
        os.makedirs(refseq_subset_dir, exist_ok=True)

        # Filter available GenBank IDs to those with corresponding PDBs
        available_genbank_ids = []
        for genbank_id in genbank_to_refseq.keys():
            for _, _, files in os.walk(args.viro3d_dir):
                if any(f"{genbank_id}_" in file for file in files):
                    available_genbank_ids.append(genbank_id)
                    break

        # Exclude already selected GenBank IDs and select additional ones
        remaining_genbank_ids = list(set(available_genbank_ids) - cumulative_genbank_ids)
        additional_genbank_ids = random.sample(
            remaining_genbank_ids, max(0, random_size - len(cumulative_genbank_ids))
        )
        cumulative_genbank_ids.update(additional_genbank_ids)

        with open(genbank_list, "w") as genbank_file:
            genbank_file.write("\n".join(cumulative_genbank_ids))

        # Copy unique PDB files for each GenBank ID
        processed_genbank_ids = set()
        for genbank_id in additional_genbank_ids:
            if genbank_id in processed_genbank_ids:
                continue  # Skip duplicates

            found = False
            for root, _, files in os.walk(args.viro3d_dir):
                for file in files:
                    if f"{genbank_id}_" in file:
                        shutil.copy(
                            os.path.join(root, file),
                            os.path.join(genbank_subset_dir, file),
                        )
                        processed_genbank_ids.add(genbank_id)
                        found = True
                        break
                if found:
                    break

        # Ensure the GenBank subset directory has the correct size
        while len(os.listdir(genbank_subset_dir)) < random_size:
            remaining_genbank_ids = list(set(available_genbank_ids) - cumulative_genbank_ids)
            if not remaining_genbank_ids:
                break  # No more GenBank IDs to add
            additional_genbank_id = random.choice(remaining_genbank_ids)
            cumulative_genbank_ids.add(additional_genbank_id)

            for root, _, files in os.walk(args.viro3d_dir):
                for file in files:
                    if f"{additional_genbank_id}_" in file:
                        shutil.copy(
                            os.path.join(root, file),
                            os.path.join(genbank_subset_dir, file),
                        )
                        break

        # Select one RefSeq ID per GenBank ID and write to RefSeq output
        additional_refseq_ids = []
        for genbank_id in additional_genbank_ids:
            if genbank_id in genbank_to_refseq:
                # Select the first RefSeq ID for this GenBank ID
                additional_refseq_ids.append(genbank_to_refseq[genbank_id][0])

        cumulative_refseq_ids.update(additional_refseq_ids)

        with open(refseq_output, "w") as refseq_file:
            refseq_file.write("\n".join(cumulative_refseq_ids))

        # Copy one matching PDB file for each RefSeq ID
        for refseq_id in additional_refseq_ids:
            found = False
            for root, _, files in os.walk(args.viral_structures_dir):
                for file in files:
                    if f"__{refseq_id}__" in file:
                        shutil.copy(
                            os.path.join(root, file),
                            os.path.join(refseq_subset_dir, file),
                        )
                        found = True
                        break
                if found:
                    break

        # Ensure the RefSeq subset directory has the correct size
        while len(os.listdir(refseq_subset_dir)) < random_size:
            remaining_refseq_ids = list(set(mapping.keys()) - cumulative_refseq_ids)
            if not remaining_refseq_ids:
                break  # No more RefSeq IDs to add
            additional_refseq_id = random.choice(remaining_refseq_ids)
            cumulative_refseq_ids.add(additional_refseq_id)

            for root, _, files in os.walk(args.viral_structures_dir):
                for file in files:
                    if f"__{additional_refseq_id}__" in file:
                        shutil.copy(
                            os.path.join(root, file),
                            os.path.join(refseq_subset_dir, file),
                        )
                        break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset GenBank IDs and corresponding PDBs.")
    parser.add_argument(
        "--genbank_lists", nargs="+", required=True, help="Paths to GenBank ID lists."
    )
    parser.add_argument(
        "--genbank_subset_dirs",
        nargs="+",
        required=True,
        help="Paths to output GenBank PDB subset directories.",
    )
    parser.add_argument(
        "--refseq_subset_dirs",
        nargs="+",
        required=True,
        help="Paths to output RefSeq PDB subset directories.",
    )
    parser.add_argument(
        "--refseq_subsets",
        nargs="+",
        required=True,
        help="Paths to output RefSeq ID lists.",
    )
    parser.add_argument(
        "--random_sizes",
        nargs="+",
        type=int,
        required=True,
        help="Subset sizes for random GenBank IDs.",
    )
    parser.add_argument(
        "--mapping_file", required=True, help="Path to RefSeq-to-GenBank mapping file."
    )
    parser.add_argument("--viro3d_dir", required=True, help="Directory containing Viro3D PDBs.")
    parser.add_argument(
        "--viral_structures_dir",
        required=True,
        help="Directory containing viral structures.",
    )
    args = parser.parse_args()
    main(args)
