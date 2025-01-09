import pandas as pd
import sys

# Input and output paths
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the TSV file
df = pd.read_csv(input_file, sep="\t")

# Extract unique virus names (assuming the column name is 'virus')
unique_viruses = df['virus'].unique()

# Save to a file
with open(output_file, 'w') as f:
    for virus in unique_viruses:
        f.write(f"{virus}\n")
