import argparse
import os
import re

import numpy as np
import pandas as pd
from sklearn.mixture import BayesianGaussianMixture
from sklearn.preprocessing import StandardScaler


def find_and_process_tsv_files(root_directory, verbose=True):
    """
    Find and process TSV files in the specified directory that match given name patterns.

    Args:
        root_directory (str): The root directory to search for TSV files.
        verbose (bool): Whether to print progress messages.

    Returns:
        dict: Dictionary of DataFrames with keys as 'grandparent_filename'.
    """
    tsv_files_dict = {}
    files_checked = 0
    files_kept = 0

    required_substrings = ["threshold0.5", "tmalignfast0", "exacttmscore1"]
    excluded_substrings = ["exhaustivesearch1"]

    if verbose:
        print(f"Starting search in: {root_directory}")

    for dirpath, _, filenames in os.walk(root_directory):
        if "processed" not in dirpath:
            continue
            
        for filename in filenames:
            files_checked += 1

            if (
                filename.endswith(".tsv")
                and all(sub in filename for sub in required_substrings)
                and not any(sub in filename for sub in excluded_substrings)
            ):
                full_path = os.path.join(dirpath, filename)

                try:
                    df = pd.read_csv(full_path, sep="\t", low_memory=False)
                    grandparent = os.path.basename(os.path.dirname(os.path.dirname(full_path)))
                    key = f"{grandparent}_{filename}"
                    tsv_files_dict[key] = df
                    files_kept += 1
                except Exception as e:
                    if verbose:
                        print(f"Failed to read {full_path}: {e}")

    if verbose:
        print(f"Checked {files_checked} files, kept {files_kept}.")

    return tsv_files_dict


def extract_prefix(query):
    """
    Extract EF- or CF- style ID with '_relaxed' suffix from a query string.

    Args:
        query (str): The query string.

    Returns:
        str or None: The matched ID or None if no match is found.
    """
    match = re.search(r"(EF-|CF-)[\w.-]+_relaxed", str(query))
    return match.group(0) if match else None


def process_all_dfs(result_dict, viro3dclusters):
    """
    Process each DataFrame in the result_dict by:
    1. Extracting a specific ID from the 'query' column.
    2. Merging with the viro3dclusters DataFrame on that ID.

    Args:
        result_dict (dict): Dictionary of DataFrames keyed by filename.
        viro3dclusters (DataFrame): DataFrame with 'cluster_member' and 'cluster_id' columns.

    Returns:
        dict: Dictionary with processed DataFrames.
    """
    processed_dict = {}

    for key, df in result_dict.items():
        df_copy = df.copy()

        df_copy["extracted_id"] = df_copy["query"].astype(str).apply(extract_prefix)

        merged_df = df_copy.merge(
            viro3dclusters[["cluster_member", "cluster_id", "genbank_name"]],
            left_on="extracted_id",
            right_on="cluster_member",
            how="left",
        )

        merged_df.drop(columns=["cluster_member"], inplace=True)
        processed_dict[key] = merged_df

    return processed_dict


def process_and_merge_df_pairs(processed_results):
    """
    Processes type1/type2 DataFrame pairs in a dictionary:
    - For type2 entries: transform e-values to -log10 and filter.
    - For type1 entries: merge e-value info from corresponding type2 DataFrames.

    Args:
        processed_results (dict): Dictionary of DataFrames keyed by filename.

    Returns:
        dict: Dictionary with updated DataFrames.
    """
    result_dict = {k: v.copy() for k, v in processed_results.items()}

    # Step 1: Transform and filter type2 entries
    for key, df in result_dict.items():
        if "alignmenttype2" in key and "evalue" in df.columns:
            df["neg_log_evalue"] = -np.log10(df["evalue"].replace(0, 1e-300))
            result_dict[key] = df[(df["alnlen"] > 20) & (df["qtmscore"] > 0.15)]

    # Step 2: Match type1 with type2 by replacing alignmenttype in key
    for type1_key in [k for k in processed_results if "alignmenttype1" in k]:
        type2_key = type1_key.replace("alignmenttype1", "alignmenttype2")
        if type2_key not in processed_results:
            continue  # Skip if no matching type2

        df1 = result_dict[type1_key]
        df2 = result_dict[type2_key]

        for df in (df1, df2):
            if "query_target_pair" not in df.columns:
                df["query_target_pair"] = df["query"] + "|||" + df["target"]

        evalue_info = df2[["query_target_pair", "evalue", "neg_log_evalue"]].copy()
        merged_df = df1.merge(evalue_info, on="query_target_pair", how="inner")
        result_dict[type1_key] = merged_df

    return result_dict


def analyze_cluster(cluster_df, cluster_num, selected_features):
    """
    Analyze a specific cluster and return its statistics.

    Args:
        cluster_df (DataFrame): The full DataFrame with cluster annotations.
        cluster_num (int): Cluster number to analyze.
        selected_features (list): List of feature names.

    Returns:
        dict: Dictionary of summary statistics for the cluster.
    """
    cluster_subset = cluster_df[cluster_df["cluster"] == cluster_num]

    def safe_stat(col, fn):
        return float(fn(cluster_subset[col])) if col in cluster_subset.columns else None

    # Core stats
    stats = {
        "size": len(cluster_subset),
        "qtmscore_mean": safe_stat("qtmscore", np.mean),
        "qtmscore_min": safe_stat("qtmscore", np.min),
        "qtmscore_max": safe_stat("qtmscore", np.max),
        "probability_min": safe_stat("cluster_probability", np.min),
        "probability_max": safe_stat("cluster_probability", np.max),
        "neg_log_evalue_mean": safe_stat("neg_log_evalue", np.mean),
        "neg_log_evalue_min": safe_stat("neg_log_evalue", np.min),
        "neg_log_evalue_max": safe_stat("neg_log_evalue", np.max),
    }

    # Convert neg_log_evalue to e-value scale
    if "neg_log_evalue" in cluster_subset.columns:
        neg_vals = cluster_subset["neg_log_evalue"].dropna()
        if not neg_vals.empty:
            stats.update(
                {
                    "evalue_min": float(10 ** -neg_vals.max()),
                    "evalue_max": float(10 ** -neg_vals.min()),
                    "evalue_mean": float(10 ** -neg_vals.mean()),
                }
            )
        else:
            stats.update({"evalue_min": None, "evalue_max": None, "evalue_mean": None})

    # Categorical lists
    stats["host_genes"] = (
        cluster_subset["host_gene_names_primary"].dropna().unique().tolist()
        if "host_gene_names_primary" in cluster_subset.columns
        else []
    )

    stats["host_functions"] = (
        cluster_subset["host_protein_names"].dropna().unique().tolist()
        if "host_protein_names" in cluster_subset.columns
        else []
    )

    stats["queries"] = (
        cluster_subset["query"].dropna().unique().tolist()
        if "query" in cluster_subset.columns
        else []
    )

    stats["genbank_names"] = (
        cluster_subset["genbank_name"].dropna().unique().tolist()
        if "genbank_name" in cluster_subset.columns
        else []
    )

    return stats


def process_all_dataframes_with_gmm(
    processed_data, viro3dclusters, min_cluster_size=2, feature_combinations=None
):
    """
    Run GMM clustering on all processed DataFrames and extract summary stats.

    Args:
        processed_data (dict): Dictionary of processed DataFrames keyed by filename.
        viro3dclusters (DataFrame): DataFrame containing virus cluster metadata.
        min_cluster_size (int): Minimum size required to consider a cluster.
        feature_combinations (list): List of feature combinations to cluster on.

    Returns:
        (DataFrame, dict): Combined summary DataFrame and full GMM results.
    """

    all_summary_data = []
    all_gmm_results = {}

    for key, df in processed_data.items():
        print(f"\nProcessing: {key}")
        alignment_match = re.search(r"alignmenttype(\d)", key)
        alignment_type = alignment_match.group(1) if alignment_match else "unknown"

        if "cluster_id" not in df.columns:
            print(f"Warning: No cluster_id column found in {key}, skipping.")
            continue

        for cluster_id in df["cluster_id"].dropna().unique():
            cluster_df = df[df["cluster_id"] == cluster_id].copy()
            cluster_df.reset_index(drop=True, inplace=True)

            if cluster_df.empty:
                continue

            if "neg_log_evalue" not in cluster_df.columns and "evalue" in cluster_df.columns:
                cluster_df["neg_log_evalue"] = -np.log10(cluster_df["evalue"].replace(0, 1e-300))
            elif "neg_log_evalue" not in cluster_df.columns:
                print(f"Warning: No evalue column found for {key}, cluster {cluster_id}")
                cluster_df["neg_log_evalue"] = float("nan")

            unique_host_genes = (
                cluster_df["host_gene_names_primary"].dropna().unique()
                if "host_gene_names_primary" in cluster_df.columns
                else []
            )

            gmm_results = {}

            if alignment_type == "1":
                feature_combinations = [
                    ["qtmscore", "alnlen"],
                    ["qtmscore", "neg_log_evalue", "alnlen"],
                ]
            else:
                feature_combinations = [["qtmscore", "neg_log_evalue", "alnlen"]]

            for feature_combination in feature_combinations:
                selected_features = [f for f in feature_combination if f in cluster_df.columns]
                if not selected_features:
                    continue

                feature_df = cluster_df.copy()
                X = feature_df[selected_features].dropna()

                if X.shape[0] < 2:
                    feature_df["cluster"] = 0
                    feature_df["cluster_probability"] = 1.0
                else:
                    try:
                        X_scaled = StandardScaler().fit_transform(X)

                        gmm = BayesianGaussianMixture(
                            weight_concentration_prior_type="dirichlet_process",
                            weight_concentration_prior=0.1,
                            n_components=min(X.shape[0], 40),
                            covariance_type="tied",
                            max_iter=3000,
                            tol=1e-5,
                            random_state=42,
                        )
                        gmm.fit(X_scaled)

                        cluster_labels = gmm.predict(X_scaled)
                        cluster_probs = gmm.predict_proba(X_scaled)

                        feature_df.loc[X.index, "cluster"] = cluster_labels
                        feature_df.loc[X.index, "cluster_probability"] = [
                            probs[label] for probs, label in zip(cluster_probs, cluster_labels)
                        ]
                    except Exception as e:
                        print(f"Error clustering {key}, cluster {cluster_id}: {e}")
                        feature_df["cluster"] = 0
                        feature_df["cluster_probability"] = 1.0

                model_id = f"{key}_Cluster_{cluster_id}_{'_'.join(feature_combination)}_cov-tied"

                cluster_stats = {
                    c: analyze_cluster(feature_df, c, selected_features)
                    for c in feature_df["cluster"].unique()
                }

                for cnum, stats in cluster_stats.items():
                    if stats["size"] == 0:
                        print(f"WARNING: Cluster {cnum} in {model_id} has 0 members!")

                gmm_results[model_id] = {
                    "merged_df": feature_df,
                    "cluster_stats": cluster_stats,
                    "original_cluster_id": cluster_id,
                    "feature_set": feature_combination,
                    "best_by": (
                        "qtmscore"
                        if alignment_type == "1"
                        and set(feature_combination) == {"qtmscore", "alnlen"}
                        else "neg_log_evalue"
                    ),
                    "total_unique_host_genes": len(unique_host_genes),
                    "source_key": key,
                    "alignment_type": alignment_type,
                }

            for model_identifier, result in gmm_results.items():
                cluster_df = result["merged_df"]
                cluster_stats = result["cluster_stats"]
                best_by = result["best_by"]

                metric_key = "qtmscore_mean" if best_by == "qtmscore" else "neg_log_evalue_mean"
                valid_clusters = {
                    cnum: stats.get(metric_key)
                    for cnum, stats in cluster_stats.items()
                    if stats.get(metric_key) is not None and stats.get("size", 0) > 0
                }

                if not valid_clusters:
                    print(f"No valid clusters for {model_identifier}")
                    continue

                sorted_clusters = sorted(valid_clusters.items(), key=lambda x: x[1], reverse=True)
                best_cluster = sorted_clusters[0][0]
                best_stats = cluster_stats[best_cluster]

                next_best_stats = (
                    cluster_stats.get(sorted_clusters[1][0]) if len(sorted_clusters) > 1 else None
                )
                qtmscore_diff = 0.0
                neg_log_evalue_diff = 0.0

                if next_best_stats:
                    if (
                        best_stats.get("qtmscore_mean") is not None
                        and next_best_stats.get("qtmscore_mean") is not None
                    ):
                        qtmscore_diff = (
                            best_stats["qtmscore_mean"] - next_best_stats["qtmscore_mean"]
                        )
                    if (
                        best_stats.get("neg_log_evalue_mean") is not None
                        and next_best_stats.get("neg_log_evalue_mean") is not None
                    ):
                        neg_log_evalue_diff = (
                            best_stats["neg_log_evalue_mean"]
                            - next_best_stats["neg_log_evalue_mean"]
                        )

                summary_entry = {
                    "source_key": result["source_key"],
                    "original_cluster_id": result["original_cluster_id"],
                    "alignment_type": result["alignment_type"],
                    "feature_set": "_".join(result["feature_set"]),
                    "best_by": best_by,
                    "best_cluster": best_cluster,
                    "best_cluster_size": best_stats.get("size"),
                    "best_qtmscore": best_stats.get("qtmscore_mean"),
                    "best_neg_log_evalue": best_stats.get("neg_log_evalue_mean"),
                    "qtmscore_difference_vs_nextbest": qtmscore_diff,
                    "neg_log_evalue_difference_vs_nextbest": neg_log_evalue_diff,
                    "cluster_members_host_genes": ", ".join(best_stats.get("host_genes", [])),
                    "cluster_members_host_functions": ", ".join(
                        best_stats.get("host_functions", [])
                    ),
                    "cluster_member_queries": ", ".join(best_stats.get("queries", [])),
                    "genbank_names": ", ".join(best_stats.get("genbank_names", [])),
                    "total_unique_host_genes": result["total_unique_host_genes"],
                    "best_cluster_unique_host_genes": len(best_stats.get("host_genes", [])),
                    "evalue_min": best_stats.get("evalue_min"),
                    "evalue_max": best_stats.get("evalue_max"),
                    "evalue_mean": best_stats.get("evalue_mean"),
                    "probability_min": best_stats.get("probability_min"),
                    "probability_max": best_stats.get("probability_max"),
                    "qtmscore_min": best_stats.get("qtmscore_min"),
                    "qtmscore_max": best_stats.get("qtmscore_max"),
                }

                all_summary_data.append(summary_entry)

            all_gmm_results.update(gmm_results)

    combined_summary_df = pd.DataFrame(all_summary_data)
    return combined_summary_df, all_gmm_results


def generate_detailed_csv(all_gmm_results, output_path="cluster_analysis_detailed.csv"):
    """
    Generate a detailed CSV containing individual data point information for the best clusters.

    Args:
        all_gmm_results (dict): Dictionary of all GMM results.
        output_path (str): Path to save the detailed CSV file.

    Returns:
        pd.DataFrame: DataFrame containing the detailed information.
    """

    def extract_protein_id(query):
        match = re.search(r"(EF-|CF-)([A-Za-z0-9]+\.[0-9]+)", str(query))
        if match:
            return match.group(2)
        return None

    detailed_data = []

    for model_identifier, result_dict in all_gmm_results.items():  # noqa: B007
        merged_df = result_dict["merged_df"]
        cluster_stats = result_dict["cluster_stats"]
        original_cluster_id = result_dict["original_cluster_id"]
        source_key = result_dict["source_key"]
        feature_set = "_".join(result_dict["feature_set"])

        best_by = result_dict.get("best_by")
        metric_key = "qtmscore_mean" if best_by == "qtmscore" else "neg_log_evalue_mean"

        cluster_means = {
            cnum: stats.get(metric_key)
            for cnum, stats in cluster_stats.items()
            if stats.get(metric_key) is not None
        }

        if not cluster_means:
            continue

        best_cluster = max(cluster_means, key=cluster_means.get)
        best_cluster_df = merged_df[merged_df["cluster"].astype(int) == int(best_cluster)].copy()
        best_cluster_analysis = analyze_cluster(merged_df, best_cluster, result_dict["feature_set"])

        grandparent_folder = source_key.split("_")[0] if "_" in source_key else ""

        # Ensure protein_id is present or compute it
        if "protein_id" not in best_cluster_df.columns and "query" in best_cluster_df.columns:
            best_cluster_df["protein_id"] = best_cluster_df["query"].apply(extract_protein_id)

        for _, row in best_cluster_df.iterrows():
            entry = {
                "source_key": source_key,
                "grandparent_folder": grandparent_folder,
                "original_cluster_id": original_cluster_id,
                "feature_set": feature_set,
                "best_cluster": best_cluster,
                "query": row.get("query", "N/A"),
                "target": row.get("target", "N/A"),
                "qtmscore": row.get("qtmscore", "N/A"),
                "neg_log_evalue": row.get("neg_log_evalue", "N/A"),
                "evalue": row.get("evalue", "N/A"),
                "alnlen": row.get("alnlen", "N/A"),
                "cluster_probability": row.get("cluster_probability", "N/A"),
                "protein_id": row.get("protein_id", "N/A"),
                "host_gene_names_primary": row.get("host_gene_names_primary", "N/A"),
                "host_protein_names": row.get("host_protein_names", "N/A"),
                "genbank_name": row.get("genbank_name", "N/A"),
                "best_cluster_qtmscore_min": best_cluster_analysis.get("qtmscore_min", "N/A"),
                "best_cluster_qtmscore_max": best_cluster_analysis.get("qtmscore_max", "N/A"),
                "best_cluster_neg_log_evalue_min": best_cluster_analysis.get(
                    "neg_log_evalue_min", "N/A"
                ),
                "best_cluster_neg_log_evalue_max": best_cluster_analysis.get(
                    "neg_log_evalue_max", "N/A"
                ),
            }
            detailed_data.append(entry)

    detailed_df = pd.DataFrame(detailed_data)

    for col in detailed_df.columns:
        detailed_df[col] = detailed_df[col].apply(lambda x: "N/A" if x is None else x)

    detailed_df.to_csv(output_path, index=False)
    return detailed_df


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="gmm selection workflow")

    parser.add_argument("--root-dir", required=True, help="Root directory to search for TSV files")

    parser.add_argument("--clusters-csv", required=True, help="Path to the viro3dclusters CSV file")

    parser.add_argument(
        "--output",
        default="combined_clusters_gmm_results.csv",
        help="Output file name for the combined results",
    )

    parser.add_argument(
        "--detailed-output",
        default="cluster_analysis_detailed.csv",
        help="Output file name for the detailed per-row results",
    )

    return parser.parse_args()


def main():
    """Main function to run the workflow"""
    # Parse command line arguments
    args = parse_arguments()

    # Step 1: Load viro3dclusters dataset
    print(f"Loading viro3dclusters from {args.clusters_csv}")
    viro3dclusters = pd.read_csv(args.clusters_csv)

    # Step 2: Find and process TSV files
    foldseek_dict = find_and_process_tsv_files(args.root_dir)
    print(f"Found {len(foldseek_dict)} matching TSV files")

    # Step 3: Process all DataFrames to extract IDs and merge with clusters
    processed_results = process_all_dfs(foldseek_dict, viro3dclusters)

    # Step 4: Process and merge DataFrame pairs
    processed_data = process_and_merge_df_pairs(processed_results)

    # Step 5: Apply GMM clustering and analyze results
    combined_summary, all_results = process_all_dataframes_with_gmm(processed_data, viro3dclusters)

    # Step 6: Generate and save detailed per-row data for best clusters
    detailed_df = generate_detailed_csv(all_results, args.detailed_output)

    # Save the results to the specified output file
    combined_summary.to_csv(args.output, index=False)
    print(f"Combined summary saved to {args.output}")

    return combined_summary, detailed_df


if __name__ == "__main__":
    main()
