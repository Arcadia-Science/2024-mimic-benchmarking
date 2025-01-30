"""
This snakefile creates standardized datasets of viral protein structures from Viro3D.
It identifies viral structures from viruses that infect humans and downloads them.
Then, it selects random subsets of those viruses. Each larger subset is inclusive of the smaller
subset.
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
        # Updated VMR with VirusHostDB merged data
        INPUT_DIRPATH / "viral" / "vmr_metadata_with_virushostdb.tsv",
        # Metadata directory fetched from Viro3D API
        expand(
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / "viro3d_metadata",
            host_organism=HOST_ORGANISMS,
        ),
        # Final merged metadata
        OUTPUT_DIRPATH / "merged_viral_metadata.tsv",
        # Final logs for checking and replacing failed downloads
        expand(
            OUTPUT_DIRPATH
            / "random_protein_sets"
            / "viral"
            / "{host_organism}"
            / "viro3d_all_pdbs_logs"
            / "check_and_replace_done.txt",
            host_organism=HOST_ORGANISMS,
        ),


rule download_kegg_virushostdb:
    """
    Downloaded 2024-12-09 version
    MD5sum 19ddaac3ffafc290793a69d367739b88
    This file includes information about which viruses (species, lineage, taxid) infect humans.
    We use the information in this file to determine which subsets of viro3d structures are
    from viruses that infect humans.
    """
    output:
        tsv=INPUT_DIRPATH / "viral" / "virushostdb.tsv",
    shell:
        """
        curl -JLo {output.tsv} https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv
        """


rule download_vmr_metadata:
    """
    VMR_MSL38_v2 version downloaded. This is what Viro3D works with. 
    This is metadata from the International Committee on Taxonomy of Viruses. Entries include the virus name, isolate designation, 
    suggested abbreviation, GenBank accession number, genome composition, and host source.
    """
    output:
        xlsx=INPUT_DIRPATH / "viral" / "vmr_metadata.xlsx",
        tsv=INPUT_DIRPATH / "viral" / "vmr_metadata.tsv",
    conda:
        "envs/xlsx2csv.yml"
    shell:
        """
        curl -JLo {output.xlsx} https://ictv.global/sites/default/files/VMR/VMR_MSL38_v2.xlsx
        xlsx2csv {output.xlsx} --delimiter '\t' > {output.tsv}
        """


rule merge_virushostdb_into_vmr:
    """
    Merge all columns from virushostdb.tsv into vmr_metadata.tsv based on matching REFSEQ IDs.
    """
    input:
        virushostdb=INPUT_DIRPATH / "viral" / "virushostdb.tsv",
        vmr_metadata=INPUT_DIRPATH / "viral" / "vmr_metadata.tsv",
    output:
        merged_vmr=INPUT_DIRPATH / "viral" / "vmr_metadata_with_virushostdb.tsv",
    shell:
        """
        python scripts/merge_virushostdb_into_vmr.py {input.vmr_metadata} {input.virushostdb} {output.merged_vmr}
        """


rule extract_unique_virus_names:
    """
    Extract unique virus names from vmr_metadata_with_virushostdb.tsv.
    """
    input:
        metadata_tsv=INPUT_DIRPATH / "viral" / "vmr_metadata_with_virushostdb.tsv",
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
    Fetch all viral structure metadata for filtered viruses from the Viro3D API.
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
        # Read each virus name directly from the input file
        while read virus; do
            # Define the response file using the raw virus name
            response_file={output.metadata_dir}/"${{virus}}.json"
        
            # Encode the virus name for the API call
            encoded_virus=$(python3 -c "import urllib.parse; print(urllib.parse.quote('''$virus'''))")
        
            # Make the API request and save the response
            curl -s -X 'GET' \
                "https://viro3d.cvr.gla.ac.uk/api/proteins/virus_name/?qualifier=$encoded_virus" \
                -H 'accept: application/json' > "$response_file"
        
            # Check if the response file is not empty and log the result
            if [[ -s "$response_file" ]]; then
                echo "Successfully fetched data for $virus" >> {output.metadata_dir}/debug.log
            else
                echo "No data fetched for $virus" >> {output.metadata_dir}/debug.log
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
            / "viro3d_{host_organism}_pdbs"
        ),
        summary=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "downloadedviro3d_pdbs.txt",
        fails=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "downloadedviro3d_pdbs_fails.txt",
    shell:
        """
        python scripts/download_all_pdbs_viro3d.py {input.metadata_dir} {output.pdb_dir} {output.summary} {output.fails}
        """


rule download_fails:
    """
    Checks the failed downloads list and viro3d_all_pdbs folder for files with 22-byte sizes and redownloads them with the opposite prefix (CF or EF).
    Updates the summary file to add the new downloads.
    """
    input:
        dir=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viro3d_{host_organism}_pdbs",
        summary=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "downloadedviro3d_pdbs.txt",
        fails=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "downloadedviro3d_pdbs_fails.txt",
        metadata_dir=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viro3d_metadata",
    output:
        done_file=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "viro3d_all_pdbs_logs"
        / "check_and_replace_done.txt",
    shell:
        """
        python scripts/download_fails.py {input.dir} {input.summary} {input.metadata_dir} {input.fails}
        find {input.dir} -type f -size 22c -exec rm -f {{}} +
        echo "Check and replace completed on $(date)" > {output.done_file}
        """


rule add_structure_file_column:
    """
    Adds a structure_file column to the summary file based on Virus Name, Chosen Method, and Record ID.
    """
    input:
        summary_file=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "downloadedviro3d_pdbs.txt",
    output:
        updated_summary_file=OUTPUT_DIRPATH
        / "random_protein_sets"
        / "viral"
        / "{host_organism}"
        / "downloadedviro3d_pdbs_updated.txt",
    shell:
        """
        python scripts/add_structure_file.py {input.summary_file} {output.updated_summary_file}
        """


rule merge_metadata:
    """
    Merges summary file with metadata from JSON files.
    """
    input:
        summary_file = OUTPUT_DIRPATH / "random_protein_sets" / "viral" / "{host_organism}" / "downloadedviro3d_pdbs_updated.txt",
        json_dir = rules.fetch_viro3d_structures_metadata.output.metadata_dir,
    output:
        merged_metadata = OUTPUT_DIRPATH / "merged_viral_metadata.tsv"
    shell:
        """
        python scripts/merge_metadata.py \
            --summary-file {input.summary_file} \
            --json-dir {input.json_dir} \
            --output-file {output.merged_metadata}
        """
