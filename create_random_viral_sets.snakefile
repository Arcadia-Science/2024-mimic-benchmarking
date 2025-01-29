"""
This snakefile creates standardized datasets of random viral protein structures from human-infecting viruses from Viro3D.
"""

from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("benchmarking_data/")
INPUT_DIRPATH = Path("inputs")

# Read in host metadata
### we are only retrieving human data in this implementation but I left the structure to retrieve for multiple organisms
host_metadata = pd.read_csv("inputs/viral/host-information.csv", header=0).set_index(
    "organism", drop=False
)
HOST_ORGANISMS = host_metadata["organism"].unique().tolist()

###########################################################
## Rules
###########################################################


rule all:
    input:
        expand(
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / "{host_organism}_summary_1000.tsv",
            host_organism=HOST_ORGANISMS,
        ),


rule select_random_files:
    input:
        files = OUTPUT_DIRPATH / "random_protein_sets" / "viral" / "{host_organism}" / "viro3d_{host_organism}_pdbs"
    output:
        dir_50 = directory(OUTPUT_DIRPATH/"random_protein_sets"/"viral"/"{host_organism}"/"selected_files_{host_organism}_50"),
        dir_100 = directory(OUTPUT_DIRPATH/"random_protein_sets"/"viral"/"{host_organism}"/"selected_files_{host_organism}_100"),
        dir_500 = directory(OUTPUT_DIRPATH/"random_protein_sets"/"viral"/"{host_organism}"/"selected_files_{host_organism}_500"),
        dir_1000 = directory(OUTPUT_DIRPATH/"random_protein_sets"/"viral"/"{host_organism}"/"selected_files_{host_organism}_1000"),
        summary = OUTPUT_DIRPATH/"random_protein_sets"/"viral"/"{host_organism}"/"{host_organism}_summary_1000.tsv"
    params:
        seed = 42,
        sizes = [50, 100, 500, 1000],
        script = "scripts/select_random_files.py"
    shell:
        "python {params.script} "
        "--input-folder '{input.files}/*' "
        "--seed {params.seed} "
        "--sizes {params.sizes} "
