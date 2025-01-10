"""
This snakefile creates standardized datasets of viral protein structures from Viro3D and Nomburg et al.
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
            / "random1000_table.tsv",
            host_organism=HOST_ORGANISMS,
        ),


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


rule download_kegg_virushostdb:
    """
    Downloaded 2024-12-09 version
    MD5sum 19ddaac3ffafc290793a69d367739b88
    """
    output:
        tsv=INPUT_DIRPATH / "viral" / "virushostdb.tsv",
    shell:
        """
        curl -JLo {output.tsv} https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
        """


rule download_ncbi_taxdump_information:
    """
    Tax dump contains files to go from a NCBI taxonomy ID to full lineage or to name.
    """
    output:
        tar=INPUT_DIRPATH / "taxdump" / "taxdump.tar.gz",
        dmp=INPUT_DIRPATH / "taxdump" / "nodes.dmp",
    params:
        taxdump_dirpath=INPUT_DIRPATH / "taxdump",
    shell:
        """
        curl -JLo {output.tar} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && \
            tar xf {output.tar} -C {params.taxdump_dirpath}
        """


rule retrieve_ncbi_taxonomy_lineages_for_viral_taxids:
    """
    The NCBI taxonomy IDs don't match between KEGG virushostdb and Nomburg et al.
    Often, virushostdb uses a strain-level taxid while Nomburg uses a species.
    """
    input:
        tsv=rules.download_kegg_virushostdb.output.tsv,
        dmp=rules.download_ncbi_taxdump_information.output.dmp,
    output:
        tsv=OUTPUT_DIRPATH / "random_protein_sets" / "viral" / "viral_lineages.tsv",
    conda:
        "envs/taxonkit.yml"
    params:
        taxdump_dirpath=INPUT_DIRPATH / "taxdump",
    shell:
        """
        taxonkit lineage \
            --data-dir {params.taxdump_dirpath} \
            --show-lineage-taxids \
            --taxid-field 1 \
            --out-file {output.tsv} \
            {input.tsv}
        """


rule download_nomburg_supplementary_table:
    output:
        xlsx=INPUT_DIRPATH / "41586_2024_7809_MOESM4_ESM.xlsx",
    conda:
        "envs/web_apis.yml"
    shell:
        """
        curl -JLo {output.xlsx} https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07809-y/MediaObjects/41586_2024_7809_MOESM4_ESM.xlsx
        """


rule filter_nomburg_viruses_by_host:
    input:
        csv=INPUT_DIRPATH / "viral" / "host-information.csv",
        xlsx=rules.download_nomburg_supplementary_table.output.xlsx,
        tsv=rules.retrieve_ncbi_taxonomy_lineages_for_viral_taxids.output.tsv,
    output:
        tsv=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_metadata0.tsv",
        txt=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_paths.txt",
    conda:
        "envs/tidyverse.yml"
    shell:
        """
        Rscript scripts/filter_nomburg_viruses_by_host.R \
            --host_organism {wildcards.host_organism} \
            --host_metadata {input.csv} \
            --nomburg_metadata {input.xlsx} \
            --virushostdb {input.tsv} \
            --output_tsv {output.tsv} \
            --output_txt {output.txt}
        """


rule bind_length_info:
    input:
        viral_metadata=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_metadata0.tsv",
        length_info=INPUT_DIRPATH / "viral" / "human-viral-lengths.tsv",
    output:
        enriched_metadata=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_metadata.tsv",
    shell:
        """
        query_col=$(head -1 {input.length_info} | awk -F'\t' '{{for (i=1; i<=NF; i++) if ($i == "query") print i}}')
        nomburg_col=$(head -1 {input.viral_metadata} | awk -F'\t' '{{for (i=1; i<=NF; i++) if ($i == "nomburg_protein_name") print i}}')
        awk -v query_col="$query_col" -v nomburg_col="$nomburg_col" 'BEGIN {{FS=OFS="\t"}} 
             NR==FNR {{lengths[$(query_col)]=$2; next}}
             FNR==1 {{print $0, "length"; next}}
             {{print $0, ($(nomburg_col) in lengths ? lengths[$(nomburg_col)] : "NA")}}' \
             {input.length_info} {input.viral_metadata} > {output.enriched_metadata}
        """


rule download_nomburg_eukaryotic_virus_structures:
    output:
        zipf=INPUT_DIRPATH / "viral" / "Nomburg_2023_structures.zip",
    shell:
        """
        curl -JLo {output} https://zenodo.org/records/10291581/files/Nomburg_2023_structures.zip?download=1
        """


rule decompress_viral_structures:
    """
    Only decompress Nomburg structures that we want to compare against a given host.
    """
    input:
        zipf=INPUT_DIRPATH / "viral" / "Nomburg_2023_structures.zip",
        txt=rules.filter_nomburg_viruses_by_host.output.txt,
    output:
        dest_dir=directory(
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / "viral_structures"
        ),
    shell:
        """
        python scripts/decompress_viral_structures.py \
            --input_txt {input.txt} \
            --zip_file {input.zipf} \
            --dest_dir {output.dest_dir}
        """


rule assess_pdbs_viral:
    """
    Calculates the quality of all Nomburg PDBs.
    """
    input:
        rules.download_proteincartography_scripts.output.txt,
        protein_structures_dir=rules.decompress_viral_structures.output.dest_dir,
    output:
        tsv=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_quality.tsv",
    conda:
        "envs/plotting.yml"
    shell:
        """
        python ProteinCartography/assess_pdbs.py \
            --input  {input.protein_structures_dir} \
            --output {output.tsv}
        """


rule extract_unique_virus_names:
    """
    Extract unique virus names from viral_structure_metadata.tsv.
    """
    input:
        metadata_tsv=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_metadata.tsv",
    output:
        unique_viruses_txt=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "unique_virus_names.txt",
    shell:
        """
        python scripts/extract_unique_viruses.py {input.metadata_tsv} {output.unique_viruses_txt}
        """


rule fetch_viro3d_structures_metadata:
    """
    Fetch all viral structure metadata for filtered Nomburg unique viruses from the Viro3D API.
    """
    input:
        unique_viruses=rules.extract_unique_virus_names.output.unique_viruses_txt,
    output:
        metadata_dir=directory(
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / "viro3d_metadata"
        ),
    shell:
        """
        mkdir -p {output.metadata_dir}
        while read virus; do
            # Replace underscores with spaces for the API query
            formatted_virus_name=$(echo "$virus" | sed 's/_/ /g')
            # Encode the virus name for the URL
            encoded_virus=$(python3 -c "import urllib.parse; print(urllib.parse.quote('''$formatted_virus_name'''))")
            # Replace underscores for safe file naming
            safe_virus_name=$(echo "$virus")

            # Debug: Print the formatted and safe virus names
            echo "Formatted virus name: $formatted_virus_name" >> {output.metadata_dir}/debug.log
            echo "Encoded virus name for URL: $encoded_virus" >> {output.metadata_dir}/debug.log
            echo "Safe virus name: $safe_virus_name" >> {output.metadata_dir}/debug.log

            response_file={output.metadata_dir}/"${{safe_virus_name}}.json"
            curl -s -X 'GET' \
                "https://viro3d.cvr.gla.ac.uk/api/proteins/virus_name/?qualifier=$encoded_virus" \
                -H 'accept: application/json' > "$response_file"

            # Debug: Log the response status
            if [[ -s "$response_file" ]]; then
                echo "Successfully fetched data for $formatted_virus_name" >> {output.metadata_dir}/debug.log
            else
                echo "No data fetched for $formatted_virus_name" >> {output.metadata_dir}/debug.log
            fi
        done < {input.unique_viruses}
        """


rule download_all_pdbs:
    """
    Download all PDB files for viruses listed in the Viro3D metadata.
    """
    input:
        metadata_dir=rules.fetch_viro3d_structures_metadata.output.metadata_dir,
    output:
        pdb_dir=directory(
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / "viro3d_all_pdbs"
        ),
    shell:
        """
        python scripts/download_all_pdbs_viro3d.py {input.metadata_dir} {output.pdb_dir}
        """


rule extract_nomburg_protein_names:
    """
    Extract nomburg_protein_name column from the viral_structure_metadata.tsv file.
    """
    input:
        metadata_tsv=rules.filter_nomburg_viruses_by_host.output.tsv,
    output:
        nomburg_protein_names_txt=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "nomburg_protein_names.txt",
    shell:
        """
        awk -F'\t' 'NR>1 {{print $1}}' {input.metadata_tsv} > {output.nomburg_protein_names_txt}
        """


rule map_refseq_to_genbank:
    """
    Map RefSeq IDs from 'nomburg_protein_name' column to GenBank IDs using Viro3D metadata.
    """
    input:
        metadata_tsv=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viral_structure_metadata.tsv",
        viro3d_metadata_dir=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viro3d_metadata",
    output:
        mapping_file=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "refseq_to_genbank_mapping.tsv",
    shell:
        """
        python scripts/map_refseq_to_genbank.py \
            --metadata {input.metadata_tsv} \
            --viro3d_metadata {input.viro3d_metadata_dir} \
            --output {output.mapping_file}
        """


rule create_genbank_and_refseq_pdb_folders:
    """
    Create random sets of 50,100,500,1000 viral proteins from Viro3D and then pull out the corresponding Nomburg sets.
    """
    input:
        mapping_file=rules.map_refseq_to_genbank.output.mapping_file,
        viro3d_dir=rules.download_all_pdbs.output.pdb_dir,
        viral_structures_dir=rules.decompress_viral_structures.output.dest_dir,
    output:
        genbank_lists=[
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / f"random{size}_genbank_ids.txt"
            for size in [50, 100, 500, 1000]
        ],
        genbank_subset_dirs=[
            directory(
                OUTPUT_DIRPATH
                / "random_protein_sets"
                / "viral"
                / "{host_organism}"
                / f"random{size}_viro3d_pdbs"
            )
            for size in [50, 100, 500, 1000]
        ],
        refseq_subset_dirs=[
            directory(
                OUTPUT_DIRPATH
                / "random_protein_sets"
                / "viral"
                / "{host_organism}"
                / f"random{size}_refseq_pdbs"
            )
            for size in [50, 100, 500, 1000]
        ],
        refseq_subsets=[
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / f"random{size}_refseq_ids.txt"
            for size in [50, 100, 500, 1000]
        ],
    params:
        script="scripts/create_genbank_and_refseq_pdb_folders.py",
        random_sizes=[50, 100, 500, 1000],
    shell:
        """
        python {params.script} \
            --mapping_file {input.mapping_file} \
            --viro3d_dir {input.viro3d_dir} \
            --viral_structures_dir {input.viral_structures_dir} \
            --genbank_lists {output.genbank_lists} \
            --genbank_subset_dirs {output.genbank_subset_dirs} \
            --refseq_subset_dirs {output.refseq_subset_dirs} \
            --refseq_subsets {output.refseq_subsets} \
            --random_sizes {params.random_sizes}
        """


rule create_random1000_table:
    """
    Create a random1000 table with RefSeq_ID, GenBank_ID, RefSeq_PDB_Filename, and GenBank_Entry_Name. Since this is the largest superset, it should also have the info for the smaller sets.
    """
    input:
        random1000_refseq_list=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "random1000_refseq_ids.txt",
        mapping_file=rules.map_refseq_to_genbank.output.mapping_file,
        viro3d_metadata_dir=rules.fetch_viro3d_structures_metadata.output.metadata_dir,
        viral_structures_dir=rules.decompress_viral_structures.output.dest_dir,
    output:
        output_table=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "random1000_table.tsv",
    params:
        script="scripts/create_random1000_table.py",
    shell:
        """
        python {params.script} \
            --random1000_refseq_list {input.random1000_refseq_list} \
            --mapping_file {input.mapping_file} \
            --viro3d_metadata_dir {input.viro3d_metadata_dir} \
            --viral_structures_dir {input.viral_structures_dir} \
            --output_table {output.output_table}
        """
