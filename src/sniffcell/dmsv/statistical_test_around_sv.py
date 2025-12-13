import numpy as np
import pandas as pd
from typing import Dict, Tuple, Iterable, Optional
from math import sqrt
from scipy.stats import mannwhitneyu

def _hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    """Hedges' g with small-sample correction (Welch)."""
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    mx, my = np.mean(x), np.mean(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    # Welch pooled SD
    sp = sqrt(((nx-1)*vx + (ny-1)*vy) / (nx + ny - 2)) if (nx+ny-2) > 0 else np.nan
    if sp == 0 or np.isnan(sp):
        return np.nan
    g = (mx - my) / sp
    # small-sample correction
    J = 1 - (3 / (4*(nx+ny) - 9)) if (nx+ny) > 2 else 1.0
    return g * J

def _cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Cliff's delta: P(X>Y) - P(X<Y), robust to non-normality."""
    # Efficient approximation via ranks: use Mann–Whitney U
    # U / (nx*ny) relates to AUC; delta = 2*AUC - 1
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan
    u, _ = mannwhitneyu(x, y, alternative="two-sided")
    auc = u / (nx * ny)
    return 2*auc - 1

def _bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR over a 1D array; NaNs stay NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    mask = ~np.isnan(p)
    m = mask.sum()
    if m == 0:
        return q
    order = np.argsort(p[mask], kind="mergesort")
    ranked = p[mask][order]
    adj = ranked * m / (np.arange(m) + 1)
    # enforce monotonicity
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    q_vals = np.minimum(adj, 1.0)
    q[mask] = q_vals[np.argsort(order)]
    return q

def _agreement_score(vals: np.ndarray, high_thr=0.7, low_thr=0.3) -> float:
    """
    Fraction of reads that agree with the dominant state.
    Uses 'high' (>=high_thr) vs 'low' (<=low_thr); excludes mid-zone.
    """
    v = vals[~np.isnan(vals)]
    if v.size == 0:
        return np.nan
    high = (v >= high_thr).sum()
    low  = (v <= low_thr).sum()
    total_used = high + low
    if total_used == 0:
        return 0.0
    return max(high, low) / total_used

def statistical_screen_cpgs(
    methyl_matrix1: pd.DataFrame,
    methyl_matrix2: pd.DataFrame,
    *,
    min_reads_per_group: int = 3,
    agreement_min: float = 0.8,     # require strong within-supporting agreement
    sd_max: Optional[float] = 0.15, # optional: supporting SD must be small
    abs_mean_diff_min: float = 0.25,
    cliffs_delta_min: float = 0.3,  # medium effect ~0.33; tweak to taste
    high_thr: float = 0.7,
    low_thr: float = 0.3
) -> pd.DataFrame:
    """
    Return a per-CpG table with robust stats and a 'keep' boolean applying consistency + effect filters.

    Columns:
      n_sup, n_oth, mean_sup, mean_oth, sd_sup, sd_oth,
      delta_mean, cliffs_delta, hedges_g, mwu_p, mwu_q, agreement_sup, keep
    """
    common = methyl_matrix1.columns.intersection(methyl_matrix2.columns)
    out = []

    for cpg in common:
        x = methyl_matrix1[cpg].dropna().values.astype(float)
        y = methyl_matrix2[cpg].dropna().values.astype(float)
        n_sup, n_oth = len(x), len(y)

        if n_sup < min_reads_per_group or n_oth < min_reads_per_group:
            out.append((cpg, n_sup, n_oth, np.nan, np.nan, np.nan, np.nan,
                        np.nan, np.nan, np.nan, np.nan, np.nan, False))
            continue

        mean_sup, mean_oth = float(np.mean(x)), float(np.mean(y))
        sd_sup, sd_oth = float(np.std(x, ddof=1)), float(np.std(y, ddof=1))
        delta_mean = mean_sup - mean_oth

        # Tests / effect sizes
        try:
            mwu_stat, mwu_p = mannwhitneyu(x, y, alternative="two-sided")
        except Exception:
            mwu_p = np.nan
        cd = _cliffs_delta(x, y)
        g  = _hedges_g(x, y)

        # Consistency in supporting reads
        agree = _agreement_score(x, high_thr=high_thr, low_thr=low_thr)

        # Hard gating
        pass_agree = (not np.isnan(agree)) and (agree >= agreement_min)
        pass_sd    = True if sd_max is None else (not np.isnan(sd_sup) and sd_sup <= sd_max)
        pass_diff  = (not np.isnan(delta_mean)) and (abs(delta_mean) >= abs_mean_diff_min)
        pass_cd    = (not np.isnan(cd)) and (abs(cd) >= cliffs_delta_min)

        keep = pass_agree and pass_sd and pass_diff and pass_cd

        out.append((cpg, n_sup, n_oth, mean_sup, mean_oth, sd_sup, sd_oth,
                    delta_mean, cd, g, mwu_p, agree, keep))

    df = pd.DataFrame(out, columns=[
        "cpg", "n_sup", "n_oth", "mean_sup", "mean_oth", "sd_sup", "sd_oth",
        "delta_mean", "cliffs_delta", "hedges_g", "mwu_p", "agreement_sup", "keep"
    ])
    # Per-SV BH over the CpGs available
    df["mwu_q"] = _bh_fdr(df["mwu_p"].values)
    return df.set_index("cpg").sort_values(["keep", "mwu_q", "cliffs_delta", "delta_mean"], ascending=[False, True, False, False])


def get_statistical_tests_around_sv(
    sv_methylation_data,
    haplotype_majority_threshold: float = 0.7,
    test_type: str = "mannwhitneyu",   # kept for compatibility but we’ll compute MWU anyway
    **screen_kwargs
):
    """
    Returns: dict[sv_id] -> DataFrame from statistical_screen_cpgs(...)
    """
    import pandas as pd
    results = {}

    for sv_id, methyl_data in sv_methylation_data.items():
        mm_sup = methyl_data["methyl_matrix_supporting"]
        mm_oth = methyl_data["methyl_matrix_other"]

        # Your existing haplotype-majority restriction logic
        sup_use, oth_use = mm_sup, mm_oth
        can_hap = (
            isinstance(mm_sup.index, pd.MultiIndex)
            and "haplotype" in (mm_sup.index.names or [])
            and isinstance(mm_oth.index, pd.MultiIndex)
            and "haplotype" in (mm_oth.index.names or [])
        )
        if can_hap and len(mm_sup) > 0:
            hap_sup = mm_sup.index.get_level_values("haplotype")
            try:
                hap_sup = hap_sup.astype(int)
            except Exception:
                hap_sup = hap_sup.map(lambda x: {"1": 1, "2": 2, "-1": -1}.get(str(x), -1))
            n_total = len(hap_sup)
            n1 = int((hap_sup == 1).sum()); n2 = int((hap_sup == 2).sum())
            major = 1 if n1 >= n2 else 2
            frac = (n1 if major == 1 else n2) / max(n_total, 1)
            if frac >= haplotype_majority_threshold:
                sup_use = mm_sup.xs(major, level="haplotype", drop_level=True)
                try:
                    oth_use = mm_oth.xs(major, level="haplotype", drop_level=True)
                except KeyError:
                    oth_use = pd.DataFrame(columns=mm_oth.columns)

        # -> Consistency-aware screen
        df_stats = statistical_screen_cpgs(
            sup_use, oth_use,
            # tweak defaults here if you like
            **screen_kwargs
        )
        results[sv_id] = df_stats

    return results
