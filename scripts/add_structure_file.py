import sys

import pandas as pd

# Input and output paths
summary_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# Load the summary file
summary_data = pd.read_csv(summary_file_path, sep="\t")

# Function to construct the structure_file
def construct_structure_file(row):
    # Convert Chosen Method to string to avoid TypeError
    chosen_method = str(row["Chosen Method"])  # Ensure it's a string

    # Determine the prefix based on the Chosen Method
    if "x" in chosen_method:
        prefix = "CF-"  # ColabFoldx
    elif "ColabFold" in chosen_method:
        prefix = "CF-"  # ColabFold
    elif "ESMFold" in chosen_method:
        prefix = "EF-"  # ESMFold
    else:
        prefix = ""

    # Construct the structure file name
    return f"{row['Virus Name']}_{prefix}{row['Record ID']}_relaxed.pdb"

# Apply the function to construct the structure_file column
summary_data["structure_file"] = summary_data.apply(construct_structure_file, axis=1)

# Save the updated summary file
summary_data.to_csv(output_file_path, sep="\t", index=False)
