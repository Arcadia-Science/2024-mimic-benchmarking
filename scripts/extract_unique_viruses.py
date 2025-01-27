import sys

import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t")

unique_viruses = df["Virus name(s)"].unique()

with open(output_file, "w") as f:
    for virus in unique_viruses:
        f.write(f"{virus}\n")
