from typing import Optional, Mapping, Union
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

def kmeans_cluster_cells(
    df: pd.DataFrame,
    n_clusters: int = 2,
    random_state: Optional[int] = 42,
    scale: bool = True,
    dmr_row: Optional[Union[pd.Series, Mapping]] = None,  # <-- NEW: one row from your DMR table
) -> pd.DataFrame:
    """
    If dmr_row is provided (expects keys/columns: best_group, best_dir, mean_best_value, mean_rest_value),
    replace numeric cluster labels with {best_group, 'Other'} based on direction or proximity.
    """
    data = df.copy()

    # Extract CpG columns (all numeric)
    cpg_cols = [c for c in data.columns if np.issubdtype(data[c].dtype, np.number)]
    X = data[cpg_cols].astype(float)

    # Handle NaN by imputing with column means
    X_imputed = X.fillna(X.mean()).dropna(axis=1, how="all")
    if X_imputed.shape[1] == 0:
        return data.assign(cluster=np.nan, celltype_or_other="Unknown")
    # Optionally scale
    if scale:
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_imputed)
    else:
        X_scaled = X_imputed.values

    # Run KMeans
    km = KMeans(n_clusters=n_clusters, n_init=20, random_state=random_state)
    labels = km.fit_predict(X_scaled)

    data["cluster"] = labels  # keep numeric for debugging

    # --- Optional mapping to cell type vs Others (only when n_clusters==2 and dmr_row is given) ---
    if dmr_row is not None and n_clusters == 2:
        # Per-read mean methylation, then cluster means
        read_mean = X_imputed.mean(axis=1)
        cluster_means = read_mean.groupby(labels).mean()  # index: 0,1

        # Pull fields (tolerant to Series or dict-like)
        best_group = str(dmr_row.get("best_group", "Unknown"))
        best_dir = dmr_row.get("best_dir", None)
        mb = dmr_row.get("mean_best_value", np.nan)
        mr = dmr_row.get("mean_rest_value", np.nan)
        mb = None if pd.isna(mb) else float(mb)
        mr = None if pd.isna(mr) else float(mr)

        # Decide which cluster is "best"
        best_cluster = None
        if isinstance(best_dir, str) and best_dir.lower() in ("hyper", "hypo"):
            if best_dir.lower() == "hyper":
                best_cluster = int(cluster_means.idxmax())
            else:  # hypo
                best_cluster = int(cluster_means.idxmin())
        elif mb is not None and mr is not None:
            # fallback: closer to reported means
            d_best = (cluster_means - mb).abs()
            d_rest = (cluster_means - mr).abs()
            # ensure exactly one best
            tentative = {int(c): ("best" if d_best[c] < d_rest[c] else "rest") for c in cluster_means.index}
            if list(tentative.values()).count("best") == 1:
                best_cluster = [c for c, role in tentative.items() if role == "best"][0]
            else:
                best_cluster = int(d_best.idxmin())
        else:
            # last resort: higher mean = best
            best_cluster = int(cluster_means.idxmax())

        # Map labels
        assignment = np.where(labels == best_cluster, best_group, "Other")
        data["celltype_or_other"] = assignment
    else:
        # No mapping requested → keep numeric labels only
        data["celltype_or_other"] = data["cluster"].astype(str)

    return data
