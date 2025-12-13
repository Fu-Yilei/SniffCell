import logging
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
import pysam
from functools import partial

import os
import logging
import pandas as pd
from src.anno.vcf_to_df import read_vcf_to_df
from src.get_base_modification_dictionary import get_base_modification_dictionary_basic_supporting_reads_df
from sklearn.metrics import pairwise_distances

logging.basicConfig(level=logging.INFO)

def find_blocks_with_readwise_coverage(df, min_fraction=0.9, min_block_size=3):
    coverage = df.notna().astype(int)
    coverage_array = coverage.values
    num_reads, num_cpgs = coverage_array.shape
    columns = list(df.columns)

    blocks = []

    for block_size in range(min_block_size, num_cpgs + 1):
        for i in range(num_cpgs - block_size + 1):
            block_cov = coverage_array[:, i:i+block_size]

            # Reads with at least one CpG covered in this block
            reads_with_any = (block_cov.sum(axis=1) > 0)
            # Among those, how many have full coverage (i.e., all CpGs covered)
            reads_with_all = (block_cov.sum(axis=1) == block_size)

            total_valid_reads = reads_with_any.sum()
            if total_valid_reads == 0:
                continue  # no informative reads here

            fully_covered_reads = (reads_with_all & reads_with_any).sum()
            fraction = fully_covered_reads / total_valid_reads

            if fraction >= min_fraction:
                blocks.append((columns[i], columns[i+block_size-1], block_size, fraction, total_valid_reads))

    return pd.DataFrame(blocks, columns=['start_CpG', 'end_CpG', 'block_size', 'fraction_covered', 'num_valid_reads'])
    
def annotate_block_read_names(df, test_df):
    coverage = test_df.notna().astype(int)
    coverage_array = coverage.values
    read_names = test_df.index

    blocks_with_reads = []

    columns = list(test_df.columns)

    for _, row in df.iterrows():
        start, end = row['start_CpG'], row['end_CpG']
        block_cols = columns[columns.index(start): columns.index(end) + 1]
        block_cov = coverage[block_cols].values

        # Filter reads with at least one covered CpG in the block
        valid_mask = (block_cov.sum(axis=1) > 0)
        fully_covered_mask = (block_cov.sum(axis=1) == len(block_cols))

        # Get read names that fully cover this block
        covered_reads = read_names[valid_mask & fully_covered_mask]
        blocks_with_reads.append(list(covered_reads))

    # Add read_names column to the original block dataframe
    df = df.copy()
    df['read_names'] = blocks_with_reads
    return df
    
def find_blocks_for_read_list(block_df, read_list, min_fraction=0.9):
    # Use only the read_name part (ignore haplotype)
    read_set = set(read_list)

    matched_blocks = []

    for _, row in block_df.iterrows():
        # Extract read names (remove haplotype level)
        block_reads = set(r if not isinstance(r, tuple) else r[0] for r in row['read_names'])

        if not block_reads:
            continue

        overlap = read_set & block_reads
        fraction = len(overlap) / len(read_set)

        if fraction >= min_fraction:
            matched_blocks.append({
                'start_CpG': row['start_CpG'],
                'end_CpG': row['end_CpG'],
                'block_size': row['block_size'],
                'fraction_of_input_reads': fraction,
                'matching_read_count': len(overlap),
                'read_names': list(overlap)
            })

    return pd.DataFrame(matched_blocks)


def plot_block_with_marked_reads(test_df, start_cpg, end_cpg, marked_read_names, ignore_haplotype=True, title=None):
    cols = list(test_df.columns)
    start_idx = cols.index(start_cpg)
    end_idx = cols.index(end_cpg) + 1
    block_cpgs = cols[start_idx:end_idx]

    sub_df = test_df[block_cpgs]
    coverage_mask = sub_df.notna().any(axis=1)
    covered_reads_df = sub_df[coverage_mask]

    # Determine marked reads
    if ignore_haplotype:
        read_names = covered_reads_df.index.get_level_values(0)
        is_marked_array = pd.Series(read_names, index=covered_reads_df.index).isin(marked_read_names)
    else:
        is_marked_array = pd.Index(covered_reads_df.index).isin(marked_read_names)

    # Assign row colors
    row_colors = pd.Series(
        ['red' if mark else 'gray' for mark in is_marked_array],
        index=covered_reads_df.index
    )

    # Plot
    sns.clustermap(
        covered_reads_df,
        row_colors=row_colors,
        col_cluster=False,
        row_cluster=False,
        cmap='viridis',
        vmin=0,
        vmax=1,
        xticklabels=True,
        yticklabels=False,
        figsize=(min(0.4 * len(block_cpgs), 16), min(0.25 * len(covered_reads_df), 16))
    )
    plt.suptitle(title or f'CpG Block: {start_cpg}–{end_cpg}', y=1.02)
    plt.show()

    
import pysam
from sklearn.cluster import KMeans
from tqdm import tqdm

def kmeans_per_haplotype(test_df, start_cpg, end_cpg, n_clusters=2, random_state=0):
    cols = list(test_df.columns)
    start_idx = cols.index(start_cpg)
    end_idx = cols.index(end_cpg) + 1
    block_cpgs = cols[start_idx:end_idx]

    # Subset block CpGs and keep only rows with at least one non-NaN
    block_df = test_df[block_cpgs]
    block_df = block_df[block_df.notna().any(axis=1)]

    # Group by haplotype
    results = []
    for hap, hap_df in block_df.groupby(level='haplotype'):
        # Drop rows with too many NaNs (optional)
        if hap == -1:   # Skip unphased haplotype
            continue
        hap_df_clean = hap_df.dropna(axis=0, thresh=1)

        if hap_df_clean.shape[0] < n_clusters:
            # print(f"Skipping haplotype {hap}: fewer rows than clusters")
            continue

        # Fill NaNs with column means or other strategy
        filled_df = hap_df_clean.T.fillna(hap_df_clean.T.mean()).T

        # Run KMeans
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
        labels = kmeans.fit_predict(filled_df.values)

        hap_result = pd.DataFrame({
            'read_name': [idx[0] for idx in hap_df_clean.index],
            'haplotype': hap,
            'cluster': labels
        })
        results.append(hap_result)

    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    
def cluster_per_haplotype_hdbscan(
    test_df: pd.DataFrame,
    start_cpg: str,
    end_cpg: str,
    min_cluster_size: int = 2,
    min_samples: int | None = None,
    metric: str = "euclidean",
    cluster_selection_epsilon: float = 0.0,
    allow_single_cluster: bool = True,
    scale: bool = False,
):
    """
    Density-based clustering of reads per haplotype across CpGs [start_cpg..end_cpg]
    using HDBSCAN (no need to specify number of clusters).

    Parameters
    ----------
    test_df : pd.DataFrame
        MultiIndex rows with levels that include 'haplotype' and typically 'read_name'.
        Columns are CpG IDs/positions; must include start_cpg..end_cpg (inclusive).
        Values should be comparable on the same scale (e.g., methylation 0–1 or 0–100).
    start_cpg, end_cpg : str
        Inclusive CpG column labels defining the block to cluster on.
    min_cluster_size : int
        Smallest cluster size HDBSCAN will consider (increase -> fewer, larger clusters).
    min_samples : int or None
        If None, defaults to min_cluster_size. Larger -> more conservative (more noise).
    metric : str
        Distance metric for HDBSCAN (e.g., 'euclidean', 'manhattan', 'cosine').
    cluster_selection_epsilon : float
        Epsilon for cluster selection (use >0 to split dense regions).
    allow_single_cluster : bool
        If True, HDBSCAN can return a single cluster when that’s the best fit.
    scale : bool
        If True, standardize features within each haplotype block (zero mean, unit var).

    Returns
    -------
    pd.DataFrame
        Columns: ['read_name','haplotype','cluster','is_noise','probability','n_clusters']
        - cluster: -1 indicates noise (unclustered).
        - probability: soft cluster membership (0–1) from HDBSCAN.
        - n_clusters: number of non-noise clusters discovered for that haplotype.
    """
    try:
        import hdbscan
    except ImportError as e:
        raise ImportError("Please install hdbscan (pip install hdbscan)") from e

    if not isinstance(test_df.index, pd.MultiIndex):
        raise ValueError("test_df.index must be a MultiIndex with a 'haplotype' level.")
    if "haplotype" not in test_df.index.names:
        raise ValueError("MultiIndex must include a 'haplotype' level.")

    cols = list(test_df.columns)
    if start_cpg not in cols or end_cpg not in cols:
        raise ValueError("start_cpg/end_cpg must be existing columns in test_df.")
    start_idx = cols.index(start_cpg)
    end_idx = cols.index(end_cpg) + 1
    block_cpgs = cols[start_idx:end_idx]

    # Keep rows with at least one observed CpG in the block
    block_df = test_df[block_cpgs]
    block_df = block_df[block_df.notna().any(axis=1)]
    if block_df.empty:
        return pd.DataFrame(columns=["read_name","haplotype","cluster","is_noise","probability","n_clusters"])

    # Helper to extract read_name robustly
    def _get_read_names(idx):
        if "read_name" in idx.names:
            pos = idx.names.index("read_name")
            return [t[pos] for t in idx]
        # fallback: assume first level is read id
        return [t[0] for t in idx]

    # Optional scaler
    if scale:
        from sklearn.preprocessing import StandardScaler

    results = []
    for hap, hap_df in block_df.groupby(level="haplotype"):
        # Drop rows that are entirely NaN in the block (already filtered), and
        # fill remaining NaNs with per-column means (within this haplotype)
        hap_df_clean = hap_df.copy()
        hap_df_clean = hap_df_clean.T.fillna(hap_df_clean.T.mean()).T

        # If still any columns all-NaN (pathological), drop them
        all_nan_cols = hap_df_clean.columns[hap_df_clean.isna().all(axis=0)]
        if len(all_nan_cols) > 0:
            hap_df_clean = hap_df_clean.drop(columns=all_nan_cols)

        # If we ended up with no features, assign single cluster
        if hap_df_clean.shape[1] == 0:
            read_names = _get_read_names(hap_df_clean.index)
            out = pd.DataFrame({
                "read_name": read_names,
                "haplotype": hap,
                "cluster": np.zeros(len(read_names), dtype=int),
                "is_noise": np.zeros(len(read_names), dtype=bool),
                "probability": np.ones(len(read_names), dtype=float),
                "n_clusters": 1
            })
            results.append(out)
            continue

        X = hap_df_clean.values

        # If no variance, assign single cluster
        if np.allclose(np.nanstd(X, axis=0), 0):
            read_names = _get_read_names(hap_df_clean.index)
            out = pd.DataFrame({
                "read_name": read_names,
                "haplotype": hap,
                "cluster": np.zeros(len(read_names), dtype=int),
                "is_noise": np.zeros(len(read_names), dtype=bool),
                "probability": np.ones(len(read_names), dtype=float),
                "n_clusters": 1
            })
            results.append(out)
            continue

        # Optional scaling within haplotype
        if scale:
            X = StandardScaler(with_mean=True, with_std=True).fit_transform(X)

        # HDBSCAN fit
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples if min_samples is not None else min_cluster_size,
            metric=metric,
            cluster_selection_epsilon=cluster_selection_epsilon,
            allow_single_cluster=allow_single_cluster,
            prediction_data=True
        ).fit(X)

        labels = clusterer.labels_
        probs = clusterer.probabilities_
        n_clusters_eff = int(np.unique(labels[labels >= 0]).size)

        # Fallback: if everything is noise, label as single cluster 0 (pragmatic)
        if n_clusters_eff == 0:
            labels = np.zeros_like(labels)
            probs = np.ones_like(probs)
            n_clusters_eff = 1
            is_noise = np.zeros_like(labels, dtype=bool)
        else:
            is_noise = labels == -1

        read_names = _get_read_names(hap_df_clean.index)
        hap_out = pd.DataFrame({
            "read_name": read_names,
            "haplotype": hap,
            "cluster": labels,
            "is_noise": is_noise,
            "probability": probs,
            "n_clusters": n_clusters_eff
        })
        results.append(hap_out)

    return pd.concat(results, ignore_index=True)

# --- Top-level function required for multiprocessing ---
def process_sv_wrapper(
    sv_dict, input_bam_path, reference_genome_path, cell_type_deconv_range,
    min_supporting_read_num, min_cpg_number, min_read_aligned_fraction,
    min_read_supporting_fraction, min_cluster_proportion, max_cluster_proportion,
    min_cluster_distance, output_file_path=None
):
    sv = pd.Series(sv_dict)
    if len(sv['supporting_reads']) < min_supporting_read_num:
        logging.info(f"Skipping SV {sv['id']} due to insufficient supporting reads")
        return None

    chromosome = sv['chr']
    start = sv['ref_start']
    end = sv['ref_end']

    input_bam = pysam.AlignmentFile(input_bam_path, "rb")
    reference_genome = pysam.FastaFile(reference_genome_path)

    region_length = end - start
    if region_length > cell_type_deconv_range * 2:
        input_bam.close()
        reference_genome.close()
        return None

    region_center = (start + end) // 2
    region_start = max(region_center - cell_type_deconv_range, 0)
    region_end = region_center + cell_type_deconv_range

    read_cpg_df = get_base_modification_dictionary_basic_supporting_reads_df(
        input_bam, reference_genome, chromosome,
        (region_start, region_end),
        min_read_length=cell_type_deconv_range/2
    )
    input_bam.close()
    reference_genome.close()

    if read_cpg_df.empty:
        return None

    haplotype_series = read_cpg_df.index.get_level_values('haplotype')
    num_phased = haplotype_series.isin([1, 2]).sum()
    num_total_reads = read_cpg_df.index.get_level_values('read_name').nunique()
    if num_phased / num_total_reads < 0.5:
        return None

    result_df = find_blocks_with_readwise_coverage(read_cpg_df, min_block_size=min_cpg_number, min_fraction=min_read_aligned_fraction)
    annotated_blocks_df = annotate_block_read_names(result_df, read_cpg_df)
    matching_blocks = find_blocks_for_read_list(
        annotated_blocks_df, list(sv['supporting_reads']),
        min_fraction=min_read_supporting_fraction
    )
    if matching_blocks is None or matching_blocks.empty:
        return None

    block_distances = []
    for _, block in matching_blocks.iterrows():
        cluster_df = kmeans_per_haplotype(read_cpg_df, block['start_CpG'], block['end_CpG'], n_clusters=2).set_index(['read_name', 'haplotype'])
        if cluster_df.empty:
            logging.info(f"No clusters formed for block {block['start_CpG']}-{block['end_CpG']}")
            continue
        read_cpg_df['cluster'] = cluster_df['cluster']
        cluster_distances = {}
        block_cpgs = read_cpg_df.columns[read_cpg_df.columns.get_loc(block['start_CpG']):read_cpg_df.columns.get_loc(block['end_CpG']) + 1]
        for hap, hap_df in read_cpg_df.loc[:, list(block_cpgs) + ['cluster']].groupby(level='haplotype'):
            if hap == -1:
                continue

            features = hap_df.drop(columns='cluster')
            clusters = hap_df['cluster']
            cluster_counts = clusters.value_counts(normalize=True)
            centroids = features.groupby(clusters).mean().dropna(axis=1)
            # print(hap_df, features, clusters, centroids)
            if centroids.shape[0] < 2 or centroids.shape[1] < min_cpg_number/2:
                continue
            # print(centroids.shape)
            dist_matrix = pairwise_distances(centroids, metric='euclidean')
            cluster_distances[hap] = pd.DataFrame(
                dist_matrix,
                index=[f'cluster_{i}' for i in centroids.index],
                columns=[f'cluster_{i}' for i in centroids.index]
            )

        block_info = {
            'chr': sv['chr'],
            'start_CpG': block['start_CpG'],
            'end_CpG': block['end_CpG'],
            'block_size': block['block_size'],
            'fraction_of_input_reads': block['fraction_of_input_reads'],
            'matching_read_count': block['matching_read_count'],
        }
        for hap in [1, 2]:
            if hap in cluster_distances:
                try:
                    block_info[f'cluster_distance_hp{hap}'] = cluster_distances[hap].loc['cluster_1.0', 'cluster_0.0']
                except KeyError:
                    block_info[f'cluster_distance_hp{hap}'] = None
            else:
                block_info[f'cluster_distance_hp{hap}'] = None
        block_distances.append(block_info)

    filtered = [b for b in block_distances if (b['cluster_distance_hp1'] or 0) >= min_cluster_distance or (b['cluster_distance_hp2'] or 0) >= min_cluster_distance]
    sorted_blocks = sorted(filtered, key=lambda b: max(b['cluster_distance_hp1'] or 0, b['cluster_distance_hp2'] or 0), reverse=True)

    used_ranges = []
    selected_blocks = []
    for b in sorted_blocks:
        brange = (b['start_CpG'], b['end_CpG'])
        if not any(not (brange[1] < r[0] or brange[0] > r[1]) for r in used_ranges):
            selected_blocks.append(b)
            used_ranges.append(brange)
        if len(selected_blocks) == 3:
            break

    sv_block_df = pd.DataFrame()
    for best_block in selected_blocks:
        cluster_df = kmeans_per_haplotype(read_cpg_df, best_block['start_CpG'], best_block['end_CpG'], n_clusters=2)
        supporting_reads = set(sv['supporting_reads'])
        cluster_assignments = cluster_df[cluster_df['read_name'].isin(supporting_reads)]
        # Build a dictionary: {read_name: (haplotype, cluster)}
        read_cluster_dict = cluster_assignments.set_index('read_name')[['haplotype', 'cluster']].to_dict(orient='index')

        for (hap, cluster_num), _ in cluster_assignments.groupby(['haplotype', 'cluster']):
            total_reads = cluster_df[cluster_df['haplotype'] == hap].shape[0]
            reads_in_cluster = cluster_df[(cluster_df['haplotype'] == hap) & (cluster_df['cluster'] == cluster_num)].shape[0]
            proportion = reads_in_cluster / total_reads if total_reads > 0 else 0
            sv_block_df = pd.concat([sv_block_df, pd.DataFrame([{
                'chr': best_block['chr'],
                'ref_start': sv['ref_start'],
                'ref_end': sv['ref_end'],
                'id': sv['id'],
                'supporting_reads_num': len(supporting_reads),
                'haplotype': hap,
                'cluster_location': f"{best_block['chr']}:{best_block['start_CpG']}-{best_block['end_CpG']}",
                'cluster_distance': best_block[f'cluster_distance_hp{hap}'],
                'cluster': cluster_num,
                'cluster_proportion': proportion,
                'haplotype_read_number': total_reads,
                'cluster_read_number': reads_in_cluster,
                'supporting_reads_assignment': read_cluster_dict,
                'supporting_reads': list(supporting_reads)
            }])], ignore_index=True)
    if output_file_path:
        with open(output_file_path, 'a') as f:
            sv_block_df.to_csv(f, index=False, header=False,sep="\t")  # header=False to avoid repeating column names

    return sv_block_df if not sv_block_df.empty else None

# --- Main wrapper ---
def assign_sv_to_single_celltype_per_region(
    sv_vcf_path, input_bam_path, reference_genome_path, output_bam_folder,
    celltype_proportion_range, min_cpg_number=5, min_read_aligned_fraction=0.7,
    min_cluster_distance=0.5, output_bam=None, min_supporting_read_num=3,
    min_read_supporting_fraction=0.5, cell_type_deconv_range=5000, threads=1, output_file_path=None
):
    input_vcf = pysam.VariantFile(sv_vcf_path)
    mosaic_sv_df = read_vcf_to_df(input_vcf)
    os.makedirs(os.path.join(output_bam_folder, "methylation_sv"), exist_ok=True)
    header_list = [
        "chr", "ref_start", "ref_end", "id", "supporting_reads_num", "haplotype",
        "cluster_location", "cluster_distance", "cluster", "cluster_proportion",
        "haplotype_read_number", "cluster_read_number", "supporting_reads_assignment",
        "supporting_reads"
    ]
    if output_file_path:
        with open(output_file_path, 'a') as f:
            f.write("\t".join(header_list) + "\n")
    process_func = partial(
        process_sv_wrapper,
        input_bam_path=input_bam_path,
        reference_genome_path=reference_genome_path,
        cell_type_deconv_range=cell_type_deconv_range,
        min_supporting_read_num=min_supporting_read_num,
        min_cpg_number=min_cpg_number,
        min_read_aligned_fraction=min_read_aligned_fraction,
        min_read_supporting_fraction=min_read_supporting_fraction,
        min_cluster_proportion=celltype_proportion_range[0],
        max_cluster_proportion=celltype_proportion_range[1],
        min_cluster_distance=min_cluster_distance,
        output_file_path=output_file_path
    )

    sv_dicts = [sv.to_dict() for _, sv in mosaic_sv_df.iterrows()]
    if threads == 1:
        results = list(map(process_func, tqdm(sv_dicts, desc="Processing SVs")))
    else:
        with Pool(threads) as pool:
            results = list(tqdm(pool.imap(process_func, sv_dicts), total=len(sv_dicts), desc="Processing SVs"))

    return pd.concat([r for r in results if r is not None], ignore_index=True)
