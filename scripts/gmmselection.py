import argparse
import os
import re

import arcadia_pycolor as apc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from sklearn.metrics import silhouette_score
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
    match = re.search(r"(EF-|CF-)([\w.-]+_relaxed)", str(query))
    return match.group(2) if match else None


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
        viro3dclusters["match_id"] = (
            viro3dclusters["cluster_member"].astype(str).apply(extract_prefix)
        )

        merged_df = df_copy.merge(
            viro3dclusters[["match_id", "cluster_id", "genbank_name"]],
            left_on="extracted_id",
            right_on="match_id",
            how="left",
        )

        merged_df.drop(columns=["match_id"], inplace=True)
        processed_dict[key] = merged_df

        before = len(df_copy)
        after = len(merged_df)

        if before != after:
            print(f"⚠️ Row count changed during merge: before = {before}, after = {after}")
        else:
            print(f"✅ Row count unchanged: {before} rows")

    return processed_dict


def process_and_merge_df_pairs(processed_results):
    """
    Processes type1/type2 DataFrame pairs in a dictionary:
    - For type2 entries (3di + AA mode): transform e-values to -log10 and filter.
    The filter thresholds were chosen empirically based on our observation that
    3di + AA mode returns many very short/low quality alignments that we do not
    want to follow up with.
    - E-value of zero (which in our data represents very low e-value) is replaced
    with 1e-300 to avoid log calculation errors. 1e-300 is very small relative
    to any other returned e-values in our data.
    - For type1 entries (TMAlign mode): merge e-value info from corresponding type2 DataFrames.
    - Every type1 entry will have a corresponding type2 entry.

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
            raise ValueError(f"No matching type2 file found for: {type1_key}")

        df1 = result_dict[type1_key]
        df2 = result_dict[type2_key]

        for df in (df1, df2):
            if "query_target_pair" not in df.columns:
                df["query_target_pair"] = df["query"] + "|||" + df["target"]

        evalue_info = df2[["query_target_pair", "evalue", "neg_log_evalue"]].copy()
        merged_df = df1.merge(evalue_info, on="query_target_pair", how="inner")
        result_dict[type1_key] = merged_df

    return result_dict


def analyze_cluster(cluster_df, cluster_num):
    """
    Analyze a specific cluster and return its statistics.

    Args:
        cluster_df (DataFrame): The full DataFrame with cluster annotations.
        cluster_num (int): Cluster number to analyze.

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


def process_all_dataframes_with_gmm(processed_data):
    """
    Run GMM clustering on all processed DataFrames and extract summary stats.
    Our main goal is to use this gmm to help us select a distinct winning cluster of 'best' hits.
    The maximum number of components was set to 40; we do not expect our known
    mimic search results to have more than 40 meaningful clusters.
    Weight_concentration_prior controls the model preference for many vs. few clusters.
    We use 0.1, which will prioritize fewer clusters. These settings were in part
    determined empirically in our initial testing of using the gmm framework and we
    don't expect slightly changing them to dramatically impact the outcomes here.

    Args:
        processed_data (dict): Dictionary of processed DataFrames keyed by filename.

    Returns:
        (DataFrame, dict): Combined summary DataFrame and full GMM results.
    """

    all_summary_data = []
    all_gmm_results = {}

    UNCLUSTERED_CLUSTER_ID = 0

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

            if "neg_log_evalue" not in cluster_df.columns:
                print(f"Warning: No evalue column found for {key}, cluster {cluster_id}")

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
                if not all(f in cluster_df.columns for f in feature_combination):
                    print(
                        f"Skipping clustering for {key}, cluster {cluster_id} — "
                        f"missing one or more features: {feature_combination}"
                    )
                    continue

                selected_features = feature_combination

                feature_df = cluster_df.copy()
                X = feature_df[selected_features].dropna()

                clustering_successful = False

                if X.shape[0] < 10:
                    feature_df["cluster"] = UNCLUSTERED_CLUSTER_ID
                    feature_df["cluster_probability"] = np.nan
                    clustering_successful = True
                    sil_score = np.nan
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

                        clustering_successful = True

                        cluster_labels = gmm.predict(X_scaled)

                        if len(set(cluster_labels)) > 1:
                            sil_score = silhouette_score(X_scaled, cluster_labels)
                        else:
                            sil_score = np.nan

                        cluster_probs = gmm.predict_proba(X_scaled)

                        feature_df.loc[X.index, "cluster"] = cluster_labels
                        feature_df.loc[X.index, "cluster_probability"] = [
                            float(probs[label].item())
                            if hasattr(probs[label], "item")
                            else float(probs[label])
                            for probs, label in zip(cluster_probs, cluster_labels)
                        ]

                    except Exception as e:
                        print(f"Error clustering {key}, cluster {cluster_id}: {e}")
                        feature_df["cluster"] = UNCLUSTERED_CLUSTER_ID
                        feature_df["cluster_probability"] = np.nan
                        clustering_successful = False

                if clustering_successful:
                    model_id = (
                        f"{key}_Cluster_{cluster_id}_" f"{'_'.join(feature_combination)}_cov-tied"
                    )

                    cluster_stats = {
                        c: analyze_cluster(feature_df, c) for c in feature_df["cluster"].unique()
                    }

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
                        "silhouette_score": sil_score,
                    }
                else:
                    print(
                        f"Skipping {key}, cluster {cluster_id}, features {feature_combination}"
                        f"due to error."
                    )

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
                qtmscore_diff = None
                neg_log_evalue_diff = None

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
                    "silhouette_score": result["silhouette_score"],
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
        best_cluster_analysis = analyze_cluster(merged_df, best_cluster)

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


def plot_3d_cluster_visualization(
    all_gmm_results,
    output_path="gmm_clusters_3d",
    foldseek_dict=None,
):
    """
    Create a 3D visualization of GMM clusters using Plotly.
    Generates one plot per GMM result (maintaining the same combinations used for clustering).

    Args:
        all_gmm_results (dict): Dictionary containing all GMM results
        output_path (str): Path for output files
        foldseek_dict (dict, optional): Original TSV data for host gene lookup

    Returns:
        dict: Dictionary of figures keyed by model identifier
    """
    apc.plotly.setup()

    # Create target-to-host-gene lookup from original TSV files
    target_to_host_gene = {}
    if foldseek_dict:
        print("Building host gene lookup from original TSV files...")
        for df in foldseek_dict.items():
            if "host_gene_names_primary" in df.columns and "target" in df.columns:
                for _, row in df.iterrows():
                    target = row.get("target")
                    host_gene = row.get("host_gene_names_primary")
                    if pd.notna(target) and pd.notna(host_gene):
                        target_to_host_gene[target] = host_gene

        print(f"Created lookup with {len(target_to_host_gene)} host gene entries")
        if target_to_host_gene:
            sample_entries = list(target_to_host_gene.items())[:3]
            print(f"Sample host gene entries: {sample_entries}")

    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    # Combine primary and secondary palettes
    my_palette = [
        apc.aster,
        apc.lime,
        apc.dragon,
        apc.wish,
    ] + list(apc.palettes.secondary)

    # Always use amber for the best cluster
    BEST_CLUSTER_COLOR = apc.amber

    # Get other colors, excluding aegean
    OTHER_COLORS = [color for color in my_palette if color != BEST_CLUSTER_COLOR]

    # Define a consistent marker size
    MARKER_SIZE = 8  # size 6 better for smaller plot versions

    # Dictionary to store figures
    figures = {}

    # Process each GMM result separately (keeping the same combinations as GMM)
    print(f"Processing {len(all_gmm_results)} GMM results...")

    for model_id, result_data in all_gmm_results.items():
        # Extract key information
        merged_df = result_data["merged_df"]
        cluster_stats = result_data["cluster_stats"]

        # Create figure
        fig = go.Figure()

        # Identify the best cluster (highest neg_log_evalue or qtmscore based on the best_by field)
        best_by = result_data.get("best_by", "neg_log_evalue")
        metric_key = "qtmscore_mean" if best_by == "qtmscore" else "neg_log_evalue_mean"

        cluster_means = {}
        for cluster_num, stats in cluster_stats.items():
            if metric_key in stats and stats[metric_key] is not None:
                cluster_means[cluster_num] = stats[metric_key]

        best_cluster = None
        if cluster_means:
            best_cluster = max(cluster_means.items(), key=lambda x: x[1])[0]

        # Get unique clusters
        unique_clusters = merged_df["cluster"].unique()

        # Create a color map - best cluster gets aegean, others get from palette
        color_idx = 0
        cluster_colors = {}

        # First assign the best cluster color
        if best_cluster is not None and best_cluster in unique_clusters:
            cluster_colors[best_cluster] = BEST_CLUSTER_COLOR

        # Then assign other colors
        for cluster_num in unique_clusters:
            if cluster_num not in cluster_colors:
                cluster_colors[cluster_num] = OTHER_COLORS[color_idx % len(OTHER_COLORS)]
                color_idx += 1

        # Plot each cluster
        for cluster_num in unique_clusters:
            # Get data for this cluster
            cluster_df = merged_df[merged_df["cluster"] == cluster_num].copy()

            # Skip if empty
            if len(cluster_df) == 0:
                continue

            # Get the color for this cluster
            color = cluster_colors[cluster_num]

            # Create hover text
            hover_texts = []

            for _, row in cluster_df.iterrows():
                # Get query and target strings
                query_str = str(row.get("query", "None"))
                viral_gene = str(row.get("genbank_name", "None"))

                # Get host gene directly from the row if available
                host_gene_str = "Unknown"  # Default value when no host gene is found

                if "host_gene_names_primary" in row and pd.notna(row["host_gene_names_primary"]):
                    host_gene_str = str(row["host_gene_names_primary"])

                # Check if this is a best cluster
                is_best = cluster_num == best_cluster
                special_note = (
                    f"<br><b>Best {best_by.replace('_', ' ').title()}" f" cluster</b>"
                    if is_best
                    else ""
                )

                # Construct hover text with the host gene info for this specific row
                hover_text = (
                    f"<b>Query TM-score:</b> {row['qtmscore']:.3f}<br>"
                    f"<b>Neg log E-value:</b> {row['neg_log_evalue']:.3f}<br>"
                    f"<b>Alignment length:</b> {row['alnlen']:.1f}<br>"
                    f"<b>Viral query accession:</b> {query_str}<br>"
                    f"<b>Viral query gene:</b> {viral_gene}<br>"
                    f"<b>Human target gene:</b> {host_gene_str}{special_note}"
                )

                hover_texts.append(hover_text)

            # Create the cluster label - mark the best cluster
            is_best = cluster_num == best_cluster
            best_metric = best_by.replace("_", " ").title()
            cluster_label = f"Cluster {cluster_num}" + (f" (Best {best_metric})" if is_best else "")

            # Add trace for this cluster
            fig.add_trace(
                go.Scatter3d(
                    x=cluster_df["alnlen"],
                    y=cluster_df["neg_log_evalue"],
                    z=cluster_df["qtmscore"],
                    mode="markers",
                    marker=dict(
                        size=MARKER_SIZE,
                        color=color,
                        opacity=0.7,
                        symbol="circle",
                        line=dict(width=1, color="DarkSlateGrey"),
                    ),
                    text=hover_texts,
                    hoverinfo="text",
                    name=cluster_label,
                )
            )

        # Update layout with complete information in title
        fig.update_layout(
            scene=dict(
                xaxis=dict(
                    title="Alignment length",
                    font=dict(size=15),
                    range=[0, 650],
                    tickfont=dict(size=13),
                ),
                yaxis=dict(
                    title="Neg log E-value",
                    font=dict(size=15),
                    range=[0, 25],
                    tickfont=dict(size=13),
                ),
                zaxis=dict(
                    title="Query TM-score", font=dict(size=15), range=[0, 1], tickfont=dict(size=13)
                ),
                aspectmode="manual",
                aspectratio=dict(x=0.7, y=0.7, z=0.7),
            ),
            autosize=True,
            showlegend=False,
            hoverlabel=dict(font=dict(size=12, color="black")),
        )
        apc.plotly.style_plot(fig, monospaced_axes="all")

        # Store figure
        figures[model_id] = fig

        # Save as interactive HTML
        # Create a sanitized filename from the model_id
        sanitized_model_id = model_id.replace("/", "_").replace(" ", "_")
        output_file = os.path.join(output_path, f"{sanitized_model_id}.html")
        fig.write_html(output_file, full_html=True, include_plotlyjs="cdn")

        # Save as static PNG
        png_file = os.path.join(output_path, f"{sanitized_model_id}.png")
        fig.write_image(png_file, width=600, height=600, scale=4)

    print(f"Created {len(figures)} visualizations")
    # Return the figures dictionary
    return figures


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

    parser.add_argument(
        "--plot-dir",
        default="figures/3d_gmm_plots/",
        help="Directory to save the 3D plot HTML files",
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
    combined_summary, all_results = process_all_dataframes_with_gmm(processed_data)

    # Step 6: Generate and save detailed per-row data for best clusters
    detailed_df = generate_detailed_csv(all_results, args.detailed_output)

    # Step 7: Generate 3D visualizations
    os.makedirs(args.plot_dir, exist_ok=True)
    print(f"Generating 3D visualizations in {args.plot_dir}...")
    figures = plot_3d_cluster_visualization(all_results, output_path=args.plot_dir)

    # Step 8: Save the results to the specified output file
    combined_summary.to_csv(args.output, index=False)
    print(f"Combined summary saved to {args.output}")

    return combined_summary, detailed_df, figures


if __name__ == "__main__":
    main()
