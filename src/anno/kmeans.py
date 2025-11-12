from typing import Optional, Mapping, Union, Sequence
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def kmeans_cluster_cells(
    df: pd.DataFrame,
    n_clusters: int = 2,
    random_state: Optional[int] = 42,
    scale: bool = True,
    dmr_row: Optional[Union[pd.Series, Mapping]] = None,
    # NEW: expected mean methylation per cluster (on per-read mean scale, 0..1)
    # e.g., {'best': 0.20, 'rest': 0.65} or [0.20, 0.65]
    expected_cluster_means: Optional[Union[Mapping[str, float], Sequence[float]]] = None,
) -> pd.DataFrame:
    """
    KMeans on numeric CpG columns with optional semi-supervision.

    Semi-supervision: if `expected_cluster_means` is provided, we seed KMeans
    with rows whose per-read mean methylation are closest to those expected means.
    This biases the solution without hard constraints.

    If `dmr_row` is provided and n_clusters==2, we map the cluster that matches
    the DMR direction (hyper/hypo) or proximity to reported means to `best_group`
    and the other to 'Other'.
    """
    data = df.copy()

    # Extract numeric CpG columns
    cpg_cols = [c for c in data.columns if np.issubdtype(data[c].dtype, np.number)]
    X = data[cpg_cols].astype(float)

    # Impute NaNs with column means; drop all-NaN columns
    X_imputed = X.fillna(X.mean()).dropna(axis=1, how="all")
    if X_imputed.shape[1] == 0:
        return data.assign(cluster=np.nan, celltype_or_other="Unknown")

    # Optional scaling
    if scale:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_imputed)
    else:
        X_scaled = X_imputed.values

    # ---- Semi-supervised centroid initialization (guided init) ----
    init_centers = None
    if expected_cluster_means is not None:
        # Normalize to a simple ordered list of targets
        if isinstance(expected_cluster_means, Mapping):
            # keep insertion order if it's an OrderedMapping; otherwise sort by value
            items = list(expected_cluster_means.items())
            # stable order: by expected mean ascending
            items.sort(key=lambda kv: kv[1])
            target_means = [float(v) for _, v in items]
        else:
            target_means = [float(v) for v in expected_cluster_means]

        if len(target_means) == n_clusters:
            # Compute per-read mean methylation (on the imputed data, before scaling)
            read_mean = X_imputed.mean(axis=1).values  # shape (n_samples,)
            # For each target mean, pick the row closest in read_mean
            chosen_idx = []
            used = set()
            for tm in target_means:
                # distances to target mean
                d = np.abs(read_mean - tm)
                # avoid reusing the same row for different centers
                # pick the nearest unused index
                order = np.argsort(d)
                pick = next(i for i in order if i not in used)
                used.add(pick)
                chosen_idx.append(pick)

            # Build init centers from the (scaled) feature vectors of those rows
            init_centers = X_scaled[chosen_idx, :]
        # else: if mismatch in counts, silently fall back to default init

    # ---- Run KMeans (respect guided init if available) ----
    if init_centers is not None:
        km = KMeans(n_clusters=n_clusters, init=init_centers, n_init=1, random_state=random_state)
    else:
        km = KMeans(n_clusters=n_clusters, n_init='auto', random_state=random_state)

    labels = km.fit_predict(X_scaled)
    data["cluster"] = labels  # numeric labels

    # ---- Optional DMR-aware mapping for n_clusters==2 ----
    if dmr_row is not None and n_clusters == 2:
        # per-read mean (same as above)
        read_mean = X_imputed.mean(axis=1)
        cluster_means = read_mean.groupby(labels).mean()  # index {0,1}

        best_group = str(dmr_row.get("best_group", "Unknown"))
        best_dir = dmr_row.get("best_dir", None)
        mb = dmr_row.get("mean_best_value", np.nan)
        mr = dmr_row.get("mean_rest_value", np.nan)
        mb = None if pd.isna(mb) else float(mb)
        mr = None if pd.isna(mr) else float(mr)

        # Decide best cluster: prefer explicit direction; else proximity to mb/mr; else higher mean
        if isinstance(best_dir, str) and best_dir.lower() in ("hyper", "hypo"):
            best_cluster = int(cluster_means.idxmax() if best_dir.lower() == "hyper"
                               else cluster_means.idxmin())
        elif mb is not None and mr is not None:
            d_best = (cluster_means - mb).abs()
            d_rest = (cluster_means - mr).abs()
            tentative = {int(c): ("best" if d_best[c] < d_rest[c] else "rest") for c in cluster_means.index}
            if list(tentative.values()).count("best") == 1:
                best_cluster = [c for c, role in tentative.items() if role == "best"][0]
            else:
                best_cluster = int(d_best.idxmin())
        else:
            best_cluster = int(cluster_means.idxmax())

        data["celltype_or_other"] = np.where(labels == best_cluster, best_group, "Other")
    else:
        data["celltype_or_other"] = data["cluster"].astype(str)

    return data
