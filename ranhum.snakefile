from pathlib import Path
import pandas as pd

OUTPUT_DIRPATH = Path("benchmarking_data/")
INPUT_DIRPATH = Path("inputs")

# Read in host metadata
### we are only retrieving human data in this implementation but I left the structure to retrieve for multiple organisms
host_metadata = pd.read_csv("inputs/viral/host-information.csv", header=0).set_index("organism", drop=False)
HOST_ORGANISMS = host_metadata["organism"].unique().tolist()
RANDOM_SUBSET = ["50", "100", "500", "1000"]

###########################################################
## Rules
###########################################################

rule all:
    """
    Final rule that ensures all outputs are generated for all host organisms.
    """
    input:
        expand(
            OUTPUT_DIRPATH / "random_protein_sets" / "{host_organism}" / "random{random_subset}_pdb_quality.tsv",
            host_organism=HOST_ORGANISMS,
            random_subset=RANDOM_SUBSET
        )

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
        # Create an empty file to use a pointer for this rule running successfully.
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

rule download_uniprot_proteome_canonical_sequence_ids:
    """
    The following rule uses the UniProt FASTA file that records one protein per gene for each of our
    reference proteomes as the source of UniProt Protein IDs. This allows us to use only the
    canonical host protein for comparisons, as these are more metadata-complete and facilitate a
    cleaner analysis.
    """
    output:
        txt=OUTPUT_DIRPATH / "random_protein_sets" / "{host_organism}" / "host_proteome_canonical_protein_ids.txt",
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

rule create_random_host_list:
    """
    Create a list of random host proteins that will be used to create random subsets.
    """
    input:
        txt=rules.download_uniprot_proteome_canonical_sequence_ids.output.txt
    output:
        txt=temp(OUTPUT_DIRPATH / "random_protein_sets" / "{host_organism}" / "temp_random.txt")
    shell:
        """
        awk 'BEGIN {{srand()}} {{print rand() "\\t" $0}}' {input.txt} | sort -k1,1n | cut -f2- > {output.txt}
        """

rule select_random_host:
    """
    Generate progressively larger subsets (50, 100, 500, 1000) of protein IDs.
    Sets are inclusive, so all 50 proteins in the 50 subset occur in the 100 subset.
    """
    input:
        txt=rules.create_random_host_list.output.txt
    output:
        txt=OUTPUT_DIRPATH / "random_protein_sets" / "{host_organism}" / "random{random_subset}_proteins.txt",
    shell:
        """
        head -n {wildcards.random_subset} {input.txt} > {output.txt}
        """

rule download_host_pdbs:
    """
    Download PDB files for each subset of protein IDs into separate subfolders for each host organism.
    """
    input:
        txt=rules.select_random_host.output.txt
    output:
        randir=directory(OUTPUT_DIRPATH / "random_protein_sets" / "{host_organism}" / "random{random_subset}_pdbs"),
    conda:
        "envs/web_apis.yml"
    shell:
        """
        python ProteinCartography/download_pdbs.py \
            --input {input.txt} \
            --output {output.randir} \
            --max-structures {wildcards.random_subset}
        """

rule assess_pdbs:
    """
    Assess the quality of PDB files for each subset and each host organism.
    """
    input:
        randir=rules.download_host_pdbs.output.randir,
    output:
        tsv=OUTPUT_DIRPATH / "random_protein_sets" / "{host_organism}" / "random{random_subset}_pdb_quality.tsv",
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py \
            --input {input.randir} \
            --output {output.tsv}
        """
