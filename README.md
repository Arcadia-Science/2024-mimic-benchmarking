# 2024-mimic-benchmarking

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

This repository contains a control dataset and a workflow for indentifying viral protein structural mimicry.
The code and data in this repo are associated with the pub, "[How do we define and find viral protein structural mimicry?](https://doi.org/10.57844/arcadia-1eu9-gcsx)."

We define this as a viral protein with one or more domains structurally resembling a host protein(s) and interacting with some of the same binding partners or substrates.
This definition allows us to use structural similarity as a starting point for identifying viral proteins that may modulate host biology.
In this repository, we sought to develop and assess a method to detect viral protein structural mimicry.


## Installation and Setup

This repository uses Snakemake to run pipelines and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
mamba env create -n benchmark --file envs/dev.yml
conda activate benchmark
```

Snakemake manages rule-specific environments via the `conda` directive and using environment files in the [envs/](./envs/) directory. Snakemake itself is installed in the main development conda environment as specified in the [dev.yml](./envs/dev.yml) file.

To start the pipeline, run:

```{bash}
snakemake --software-deployment-method conda -j 8
```

## Benchmarking data

We curated a set of viral proteins to assess the method we developed to detect viral protein structural mimicry.
These proteins either:
1. Have experimental evidence for mimicry of a host protein. These proteins have a known correct match or matches to proteins in a host proteome.
2. Share common folds with host proteins. These proteins do not have a known correct match to a protein in the host proteome.

We downloaded viral protein structures and metadata for relevant proteins from [Viro3D](https://viro3d.cvr.gla.ac.uk).
We provide this data in [benchmarking_data](./benchmarking_data).
We also provide host target structures for mimicked proteins.
We downloaded these from [Alphafold](https://alphafold.ebi.ac.uk/).

## Benchmarking tests

We designed the [control_benchmarking.snakefile](./control_benchmarking.snakefile) to test which structural comparison tools maximized true positive hits while minimizing false positive hits.
We used structures of viral proteins and compared each structure against all human proteins (UniProt proteome with single structures in AlphaFold database).
See a list of all viral proteins we used [here](./benchmarking_data/controls/control_metadata.tsv).

We used Foldseek in 3Di+AA mode and TM-align mode to compare viral structures to human structures.
We used the following parameters:
* `--alignment-type`
    * `1`: TM-align mode. TM-align mode uses a structural superposition approach based on backbone geometry to measure structure similarity.
    * `2`: 3Di + AA mode.  3Di+AA uses a hybrid alignment approach that encodes 3D geometry and amino acid identity to measure structure similarity.
* `--tmalign-fast 0`. This parameter only impacts TM-align mode. Turning fast mode OFF increases the accuracy of TM-score estimates.
* `--exact-tmscore 1`: This parameter increases the accuracy of TM-score estimates.
* `--tmscore-threshold`: Alignment TM-score at which to threshold.
    * `0`: TM-score threshold is set to 0. When combined with `--exhaustive-search 1`, this returns all combinations of query vs. database. 
    * `0.5`: TM-score threshold is set to 0.5. When combined with `--exhaustive-search 0`, only hits with an alignment TM-score above this threshold are returned. 

We compare the performance of Foldseek in 3Di+AA and TM-align modes combined with downstream statistical modeling to determine the most accurate approach for detecing viral protein structural mimicry.
For more information, see the pub, "[How do we define and find viral protein structural mimicry?](https://doi.org/10.57844/arcadia-1eu9-gcsx)." 

## Benchmarking outputs

The [benchmarking snakefile](./control_benchmarking.snakefile) outputs files in the following directory structure.
All outputs are contained in an outputs/human directory.
This directory then contains the following subdirectories:
 
* `foldseek/`: Contains Foldseek results from running the commands described in the previous section.
    * `{control}`: one of bcl2, c1l, c1lpt1, c1lpt2, c4bp, ccr1, ccr2, cd47, chemokine, eif2a, helicase, ifngr, il10, il18bp, kinase, lfg4, nsp16 or nsp5. We named each folder after the human protein that the virus mimics when possible, and when not, we use protein family name or the viral protein name. We ran each Foldseek comparison separately on each set of controls (ex. we queried with all Bcl-2 viral mimics).
        * `raw/`: Raw Foldseek results in TSV format.
        * `processed/`: Processed Foldseek results in TSV format. These files include corrected alignment TM-scores and additional metadata for both viral (query) and human (target) proteins.
* `selected_mimics/`: contains results from running Gaussian mixture modeling on the Foldseek results to predict viral structural mimicry.
    * `gmmviro3d_benchmarking041725.csv`: Summary of results from statistical modeling. Columns include:
        * `source_key`: The input processed Foldseek results file that Gaussian mixture modeling is based on.
        * `original_cluster_id`: Viro3D database structure/sequence cluster identifier.
        * `alignment_type`: Foldseek alignment type. Either 1 (TM-align mode) or 2 (3Di+AA mode).
        * `feature_set`: Set of features used to perform Gaussian mixture modeling.
        * `best_by`: The feature used to select the "best" cluster from Gaussian mixture modeling.
        * `best_cluster`: The "best" cluster identifier from Gaussian mixture modeling, representing viral proteins that potentially mimic host proteins.
        * `best_cluster_size`: The number of viral proteins in the best cluster.
        * `best_qtmscore`: The highest query TM-score for a viral-host protein match in the best cluster.
        * `best_neg_log_evalue`: The highest negative log E-value (corresponding to the lowest E-value) for a viral-host protein match in the best cluster.
        * `qtmscore_difference_vs_nextbest`: The difference between the mean query TM-score in the best cluster and the mean query TM-score in the next best cluster. This column highlights how separation there is between the best cluster and other clusters.
        * `neg_log_evalue_difference_vs_nextbest`: The difference between the mean negative log E-value in the best cluster and the mean negative log E-value in the next best cluster. This column highlights how separation there is between the best cluster and other clusters.
        * `cluster_members_host_genes`: List of distinct host genes (proteins) that are matched by viral proteins in the best cluster.
        * `cluster_members_host_functions`: List of distinct host gene (protein) functions for the proteins that are matched by the viral proteins in the best cluster.
        * `cluster_member_queries`: List of viral protein Viro3D accessions for the viral proteins in the best cluster.
        * `genbank_names`: List of GenBank names for the viral proteins in the best cluster.
        * `total_unique_host_genes`: The total number of distinct host genes matched by any viral protein in the Viro3D cluster used for Gaussian mixture modeling.
        * `best_cluster_unique_host_genes`: The total number of distinct host genes matched by viral proteins in the best cluster.
        * `evalue_min`: The minimum E-value of a viral-host protein match in the best cluster.
        * `evalue_max`: The maximum E-value of a viral-host protein match in the best cluster.
        * `evalue_mean`: The mean E-value of all viral-host protein matches in the best cluster.
        * `probability_min`: The minimum probability of a viral-host protein match in the best cluster.
        * `probability_max`: The maximium probability of a viral-host protein match in the best cluster.
        * `qtmscore_min`: The minimum query TM-score of a viral-host protein match in the best cluster.
        * `qtmscore_max`: The maximum query TM-score of a viral-host protein match in the best cluster.
        * `silhouette_score`: Measures how similar an object is to its own cluster compared to other clusters. It ranges from –1 to +1, where +1 indicates that points in the best cluster are well-matched to other points in their own cluster and far from those not in the cluster (good quality cluster). -1 indicates that points in the best cluster are better matched to points in another cluster (bad quality cluster).
    * `gmmviro3d_benchmarking041725_detailed.csv`: Foldseek results for each viral protein selected by Gaussian mixture modeling, as well as Gaussian mixture modeling statistics and additional metadata.


## Repository folder structure

* [benchmarking_data/](./benchmarking_data): Structure files (PDB format) and metadata for the structures used as control files for the benchmarking workflow.
* [envs/](./envs): This repository uses conda to manage software installations and versions. All software packages required for running the snakefiles and notebooks in this repo are recorded in this folder.
* [figures/](./figures/): This directory contains notebooks that created some figures in the pub as well as figure files generated by code in the [scripts/](./scripts) folder.
* [inputs/](./inputs/): Input files required to run the Snakefiles in this repository.
* [scripts/](./scripts): Python and R scripts used by the Snakefiles or to generate figures in this repository.
* [`LICENSE`](./LICENSE): License specifying the re-use terms for the code in this repository.
* [`README.md`](./README.md): File outlining the contents of this repository.
* [`control_benchmarking.snakefile`](./control_benchmarking.snakefile): Snakefile that executes the benchmarking workflow to assess how to detect viral protein structural mimicry.
* [.github/](./.github), [.vscode/](./.vscode), [.gitignore](./.gitignore), [.pre-commit-config.yaml](./.pre-commit-config.yaml), [Makefile](./Makefile), [pyproject.toml](./pyproject.toml): Files that control the developer behavior of the repository.

## Compute specifications

We ran the [control_benchmarking.snakefile](control_benchmarking.snakefile) on an AWS EC2 instance type m5a.4xlarge running Ubuntu (20.04) 20240122 (AMI ID ami-07eb000b3340966b0).
The pipeline runs in about a 24 hours.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
