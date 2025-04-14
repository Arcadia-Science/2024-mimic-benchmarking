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


CONTROLS = [
    "bcl2",
    "c4bp",
    "cd47",
    "chemokine",
    "eif2a",
    "ifngr",
    "il10",
    "il18bp",
    "kinase",
    "lfg4",
    "nsp16",
    "nsp5",
    "ccr1",
    "ccr2",
    "helicase",
    "c1l",
    "c1lpt1",
    "c1lpt2",
]

###########################################################
## Download ProteinCartography scripts
###########################################################


rule download_proteincartography_scripts:
    """
    ProteinCartography (https://github.com/Arcadia-Science/ProteinCartography) contains many scripts
    for interacting with UniProt and AlphaFold. Because PC is not pip-installable, we download the
    scripts and environments we need in this workflow. We take this approach instead of making a
    copy of each script inside this repo. Technically, this should be broken out into many rules
    (one per file), but that seemed unnecessarily verbose so I confined it to one rule. 
    
    Note the envs that these scripts require (envs/plotting.yml, envs/web_apis.yml) need to already
    be present in the repo so they are duplicated.

    An alternative to this approach would be to use
    [Git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules). We may update to this
    approach in the future.
    """
    output:
        # Create an empty file to use a pointer for this rule runnign successfully.
        # This will allow us to refer to the ProteinCartography scripts by name/filepath instead of
        # by snakemake output syntax.
        txt=touch("scripts/ProteinCartography_scripts_downloaded.txt"),
        api_utils="ProteinCartography/api_utils.py",
        artifact_generation_utils="ProteinCartography/tests/artifact_generation_utils.py",
        assess_pdbs="ProteinCartography/assess_pdbs.py",
        color_utils="ProteinCartography/color_utils.py",
        constants="ProteinCartography/constants.py",
        download_pdbs="ProteinCartography/download_pdbs.py",
        fetch_accession="ProteinCartography/fetch_accession.py",
        fetch_uniprot_metadata="ProteinCartography/fetch_uniprot_metadata.py",
        file_utils="ProteinCartography/file_utils.py",
        map_refseq_ids="ProteinCartography/map_refseq_ids.py",
        mocks="ProteinCartography/tests/mocks.py",
    params:
        commit="88160fcf098347a29124488f445ed1d9ad72bc12",
    shell:
        """
        curl -JLo {output.api_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/api_utils.py
        curl -JLo {output.artifact_generation_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/artifact_generation_utils.py
        curl -JLo {output.assess_pdbs} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/assess_pdbs.py
        curl -JLo {output.color_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/color_utils.py
        curl -JLo {output.constants} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/constants.py
        curl -JLo {output.download_pdbs} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/download_pdbs.py
        curl -JLo {output.fetch_accession} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/fetch_accession.py
        curl -JLo {output.fetch_uniprot_metadata} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/fetch_uniprot_metadata.py
        curl -JLo {output.file_utils} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/file_utils.py
        curl -JLo {output.map_refseq_ids} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/map_refseq_ids.py
        curl -JLo {output.mocks} https://raw.githubusercontent.com/Arcadia-Science/ProteinCartography/{params.commit}/ProteinCartography/tests/mocks.py
        """


#####################################################################
## Download host proteome structures & metadata
## ----------------------------------------------
## This section downloads all 20,656 protein structurs associated
## with the canonical sequences of human UniProt proteins.
## For most of our benchmarking, we compare against all 20k human
## proteins. This helps us assess completeness of results.
#####################################################################


rule download_uniprot_proteome_canonical_sequence_ids:
    """
    The following rule uses the UniProt FASTA file that records one protein per gene for each of our
    reference proteomes as the source of UniProt Protein IDs. This allows us to use only the
    canonical host protein for comparisons, as these are more metadata-complete and facilitate a
    cleaner analysis.
    """
    output:
        txt=OUTPUT_DIRPATH / "{host_organism}" / "host_proteome_canonical_protein_ids.txt",
    conda:
        "envs/seqkit.yml"
    params:
        uniprot_proteome_id=lambda wildcards: host_metadata.loc[
            wildcards.host_organism, "uniprot_proteome_id"
        ],
        taxon_id=lambda wildcards: host_metadata.loc[wildcards.host_organism, "taxon_id"],
    shell:
        """
        curl -JL https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{params.uniprot_proteome_id}/{params.uniprot_proteome_id}_{params.taxon_id}.fasta.gz | \
            seqkit seq --only-id --name | cut -d'|' -f2 > {output.txt}
        """


rule download_host_pdbs:
    """
    Download all PDB files from AlphaFold.
    While this outputs many PDBs, we don't have any operations in Snakemake where the snakefile
    needs to be aware of all of the PDB accessions. Therefore, instead of treating this as a
    checkpoint and complicating the DAG, we'll only designate the directory as output.
    """
    input:
        py=rules.download_proteincartography_scripts.output.txt,
        txt=rules.download_uniprot_proteome_canonical_sequence_ids.output.txt,
    output:
        protein_structures_dir=directory(
            OUTPUT_DIRPATH / "{host_organism}" / "host_proteome_pdb_structures"
        ),
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/download_pdbs.py \
            --input {input.txt} \
            --output {output.protein_structures_dir} \
            --max-structures 100000
        """


#####################################################################
## Benchmark foldseek
#####################################################################

ALIGNMENT_TYPE = ["1", "2"]  # 1: TMAlign, 2: 3di + AA
TMALIGN_FAST = ["0"]  # 0: off, 1: on
EXACT_TMSCORE = ["1"]  # 0: off, 1: on
TMSCORE_THRESHOLD = ["0", "0.5"]


rule benchmark_foldseek_against_human_proteome:
    input:
        pdbs="benchmarking_data/controls/{control}/",
        protein_structures_dir=rules.download_host_pdbs.output.protein_structures_dir,
    output:
        tsv=OUTPUT_DIRPATH
        / "{host_organism}"
        / "foldseek"
        / "{control}"
        / "raw"
        / "foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}.tsv",
    conda:
        "envs/foldseek.yml"
    benchmark:
        "benchmarks/{host_organism}/foldseek/{control}/foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}.tsv"
    threads: 7
    params:
        exhaustive_search=lambda wildcards: "0" if wildcards.tmscore_threshold == "0.5" else "1",
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
            --exhaustive-search {params.exhaustive_search} \
            --tmscore-threshold {wildcards.tmscore_threshold} \
            --format-output query,target,qlen,tlen,alnlen,alntmscore,qtmscore,ttmscore,lddt,prob,qcov,tcov,pident,bits,evalue,qstart,qend,tstart,tend \
            --format-mode 4 \
            --threads {threads}
        """


rule combine_foldseek_results_with_metadata:
    input:
        foldseek_tsv=rules.benchmark_foldseek_against_human_proteome.output.tsv,
        human_metadata_csv=INPUT_DIRPATH / "human_metadata_combined.csv.gz",
        host_lddt_tsv=INPUT_DIRPATH / "human_proteome_pdb_structure_quality.tsv",
        query_metadata_tsv="benchmarking_data/merged_viral_metadata_human.tsv",
        query_lddt_tsv="benchmarking_data/controls/control_plddt.tsv",
    output:
        tsv=OUTPUT_DIRPATH
        / "{host_organism}"
        / "foldseek"
        / "{control}"
        / "processed"
        / "foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}.tsv",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/combine_results_with_metadata.R \
            --input_results {input.foldseek_tsv} \
            --input_human_metadata {input.human_metadata_csv} \
            --input_host_lddt {input.host_lddt_tsv} \
            --input_query_metadata {input.query_metadata_tsv} \
            --input_query_lddt {input.query_lddt_tsv} \
            --output {output.tsv}
        """


rule download_viro3d_cluster_info:
    output:
        csv=INPUT_DIRPATH / "viral" / "foldseekCluster90_1E-5_meta.csv",
    shell:
        """
        curl -JLo {output} https://raw.githubusercontent.com/ulad-litvin/viro3d-analysis/refs/heads/main/4_structural_clustering/foldseekCluster90_1E-5/foldseekCluster90_1E-5_meta.csv
        """


rule run_gmmselection:
    input:
        csv=rules.download_viro3d_cluster_info.output.csv,
        tsvs=expand(
            OUTPUT_DIRPATH
            / "{{host_organism}}"
            / "foldseek"
            / "{control}"
            / "processed"
            / "foldseek_alignmenttype{alignment_type}_tmalignfast{tmalign_fast}_exacttmscore{exact_tmscore}_tmscorethreshold{tmscore_threshold}.tsv",
            control=CONTROLS,
            alignment_type=ALIGNMENT_TYPE,
            tmalign_fast=TMALIGN_FAST,
            exact_tmscore=EXACT_TMSCORE,
            tmscore_threshold=TMSCORE_THRESHOLD,
        ),
    output:
        csv=OUTPUT_DIRPATH
        / "{host_organism}"
        / "selected_mimics"
        / "gmmviro3d_benchmarking041425.csv",
        csv_detailed=OUTPUT_DIRPATH
        / "{host_organism}"
        / "selected_mimics"
        / "gmmviro3d_benchmarking041425_detailed.csv",
    conda:
        "envs/scikitlearn.yml"
    shell:
        """
        python scripts/gmmselection.py \
            --root-dir outputs/human/foldseek \
            --clusters-csv {input.csv} \
            --output {output.csv} \
            --detailed-output {output.csv_detailed}
        """


rule all:
    default_target: True
    input:
        expand(
            rules.run_gmmselection.output.csv,
            host_organism=HOST_ORGANISMS,
        ),
