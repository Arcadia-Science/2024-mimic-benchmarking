from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("outputs")
INPUT_DIRPATH = Path("inputs")

# Read in host metadata, which we'll use to link the organism name to identifiers like taxid and
# uniprot proteome id. Setting the index allows us to use the organism name as a look up to find
# the correct value for other metadata.
host_metadata = pd.read_csv("inputs/viral/host-information.csv", header=0).set_index(
    ["organism"], drop=False
)
host_metadata = host_metadata[host_metadata["organism"] == "human"]
HOST_ORGANISMS = host_metadata["organism"].unique().tolist()

#####################################################################
## Benchmark foldseek
#####################################################################

ALIGNMENT_TYPE = ["1", "2"]  # 1: TMAlign, 2: 3di + AA
TMALIGN_FAST = ["0"]  # 0: off
EXACT_TMSCORE = ["1"]  # 1: on
TMSCORE_THRESHOLD = ["0", "0.25", "0.5"]
EXHAUSTIVE_SEARCH = ["0"]  # 0: off


rule benchmark_foldseek_against_viral_proteins:
    input:
        pdbs="benchmarking_data/human_query/",
        # downloaded by get_all_viral_structures.snakefile 
        protein_structures_dir="benchmarking_data/random_protein_sets/viral/{host_organism}/viro3d_{host_organism}_pdbs"
    output:
        tsv=OUTPUT_DIRPATH
        / "{host_organism}"
        / "foldseek"
        / "human_query"
        / "raw"
        / "foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}_exhaustivesearch{exhaustive_search}.tsv",
    conda:
        "envs/foldseek.yml"
    benchmark:
        "benchmarks/{host_organism}/foldseek/human_query/foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}_exhaustivesearch{exhaustive_search}.tsv"
    threads: 7
    shell:
        """
        foldseek easy-search \
            {input.pdbs} \
            {input.protein_structures_dir} \
            {output.tsv} \
            tmp_foldseek \
            -e inf \
            --max-seqs 21000 \
            --alignment-type {wildcards.alignment_type} \
            --tmalign-fast {wildcards.tmalign_fast} \
            --exact-tmscore {wildcards.exact_tmscore} \
            --exhaustive-search {wildcards.exhaustive_search} \
            --tmscore-threshold {wildcards.tmscore_threshold} \
            --format-output query,target,qlen,tlen,alnlen,alntmscore,qtmscore,ttmscore,lddt,prob,qcov,tcov,pident,bits,evalue,qstart,qend,tstart,tend \
            --format-mode 4 \
            --threads {threads}
        """


rule combine_foldseek_results_with_metadata:
    input:
        foldseek_tsv=rules.benchmark_foldseek_against_viral_proteins.output.tsv,
        human_metadata_csv=INPUT_DIRPATH / "human_metadata_combined.csv.gz",
        host_lddt_tsv=INPUT_DIRPATH / "human_proteome_pdb_structure_quality.tsv",
        viral_metadata_tsv="benchmarking_data/merged_viral_metadata.tsv",
    output:
        tsv=OUTPUT_DIRPATH
        / "{host_organism}"
        / "foldseek"
        / "human_query"
        / "processed"
        / "foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}_exhaustivesearch{exhaustive_search}.tsv",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_results_with_metadata.R \
            --input_results {input.foldseek_tsv} \
            --target viral \
            --input_human_metadata {input.human_metadata_csv} \
            --input_host_lddt {input.host_lddt_tsv} \
            --input_viral_metadata {input.viral_metadata_tsv} \
            --output {output.tsv}
        """


rule all:
    default_target: True
    input:
        expand(rules.combine_foldseek_results_with_metadata.output.tsv,
               host_organism = HOST_ORGANISMS,
               alignment_type=ALIGNMENT_TYPE,
               tmalign_fast=TMALIGN_FAST,
               exact_tmscore=EXACT_TMSCORE,
               tmscore_threshold=TMSCORE_THRESHOLD,
               exhaustive_search=EXHAUSTIVE_SEARCH,
           )
