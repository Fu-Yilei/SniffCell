import re
import numpy as np
import pandas as pd

# --- from your code ---
def _logit(p, eps=1e-6):
    p = np.clip(p, eps, 1.0 - eps)
    return np.log(p/(1.0 - p))

def _inv_logit(z):
    return 1.0/(1.0 + np.exp(-z))

def _harmony_mean_prob(M_region_cols: pd.DataFrame, eps=1e-6) -> pd.Series:
    """Stable combiner: mean on logit scale over AVAILABLE subtype columns."""
    if M_region_cols.shape[1] == 0:
        return pd.Series(np.nan, index=M_region_cols.index)
    Z = _logit(M_region_cols.to_numpy(float), eps=eps)  # B×S
    z_bar = np.nanmean(Z, axis=1)                       # B
    out = _inv_logit(z_bar)
    return pd.Series(out, index=M_region_cols.index)


def build_major_reference_in_region(
    M_df: pd.DataFrame,
    cpg_index: pd.DataFrame,           # must have ['chr','start','end']
    major_to_subtypes: dict,
    region: tuple,                     # (chrom, start, end)
    *,
    overlap_mode: str = "overlap",
    eps: float = 1e-6
):
    """
    Slice M_df to region using cpg_index coordinates,
    then harmonize detailed subtypes -> MAJOR.
    """
    region_chr, region_start, region_end = region

    # 1) Region mask
    same_chr = (cpg_index["chr"].astype(str).values == str(region_chr))
    if overlap_mode == "within":
        mask = same_chr & (cpg_index["start"] >= region_start) & (cpg_index["end"] <= region_end)
    elif overlap_mode == "overlap":
        mask = same_chr & (cpg_index["start"] < region_end) & (cpg_index["end"] > region_start)
    else:
        raise ValueError("overlap_mode must be 'overlap' or 'within'.")

    if not mask.any():
        raise ValueError("No blocks in the specified region.")

    # 2) Keep only subtype cols present
    present = {maj: [c for c in cols if c in M_df.columns] for maj, cols in major_to_subtypes.items()}
    subtype_union = sum(present.values(), [])
    M_region_det = M_df.loc[mask, subtype_union]
    nnz = M_region_det.notna().sum(axis=0)

    present = {maj: [c for c in cols if nnz.get(c, 0) > 0] for maj, cols in present.items()}
    present = {maj: cols for maj, cols in present.items() if len(cols) > 0}
    if not present:
        raise ValueError("After filtering, no informative majors remain in this region.")

    # 3) Harmonize -> major
    majors = {maj: _harmony_mean_prob(M_region_det[cols], eps=eps) for maj, cols in present.items()}
    P_major_df = pd.DataFrame(majors)
    good_mask = P_major_df.notna().any(axis=1)

    return P_major_df.loc[good_mask], M_region_det.loc[good_mask], present


import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Any

def find_consecutive_dmrs_length_biased(
    df: pd.DataFrame,
    cellA_col: str = "Neuron",
    cellB_col: str = "Oligodendrocyte",
    chrom_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
    *,
    # significance knobs
    min_delta: float = 0.30,         # legacy per-block effect threshold (used if min_score=None)
    alpha: float = 0.25,             # concavity for length bias; 0=no bias, 0.25 is mild, 0.5 is stronger
    min_score: float | None = None,  # score threshold; default = min_delta (so 1bp matches legacy)
    # stitching knobs
    min_run_blocks: int = 2,
    max_gap_bp: int = 25,
    min_run_score: float | None = None,  # optional: require a minimum total run score
    keep_block_lists: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Identify DMRs as runs of consecutive, same-direction, significant blocks,
    with a concave length bias: score = |Δ| * (length_bp ** alpha).

    Returns:
      dmrs_df: stitched DMRs with summary stats (including run_score)
      blocks_ann: input df + ['delta','abs_delta','length_bp','score','direction','is_sig']
    """
    if min_score is None:
        # So a 1-bp block needs |Δ| >= min_delta to pass, matching old behavior
        min_score = float(min_delta)

    # --- per-block stats ---
    blocks = df.copy()
    blocks["delta"] = blocks[cellA_col] - blocks[cellB_col]
    blocks["abs_delta"] = blocks["delta"].abs()
    blocks["length_bp"] = (blocks[end_col].astype(int) - blocks[start_col].astype(int)).clip(lower=1)

    # concave length bias
    # NOTE: alpha in (0,1) -> concave; alpha=0.25 is a mild preference for longer blocks
    blocks["score"] = blocks["abs_delta"] * (blocks["length_bp"] ** float(alpha))

    # direction label
    dir_AgtB = f"{cellA_col}>{cellB_col}"
    dir_BgtA = f"{cellB_col}>{cellA_col}"
    blocks["direction"] = np.where(blocks["delta"] > 0, dir_AgtB,
                            np.where(blocks["delta"] < 0, dir_BgtA, "tie"))

    # significance using score
    blocks["is_sig"] = (blocks["score"] >= float(min_score)) & (blocks["direction"] != "tie")

    # --- stitch by chromosome, preserving order ---
    dmrs: List[Dict[str, Any]] = []

    for chrom, sub in blocks.groupby(chrom_col, sort=False):
        sub = sub.sort_values(start_col).reset_index()  # keep original index in 'index'

        run_idx: List[int] = []
        run_dir: str = None
        prev_end = None

        def flush_run():
            if len(run_idx) >= min_run_blocks:
                run = sub.loc[run_idx]
                run_score = float(run["score"].sum())
                if (min_run_score is None) or (run_score >= float(min_run_score)):
                    dmrs.append({
                        chrom_col: chrom,
                        start_col: int(run[start_col].min()),
                        end_col: int(run[end_col].max()),
                        "n_blocks": int(len(run_idx)),
                        "direction": run["direction"].iloc[0],
                        "run_score": run_score,
                        # summary stats
                        "mean_abs_delta": float(run["abs_delta"].mean()),
                        "median_abs_delta": float(run["abs_delta"].median()),
                        "max_abs_delta": float(run["abs_delta"].max()),
                        "mean_score": float(run["score"].mean()),
                        "sum_length_bp": int(run["length_bp"].sum()),
                        "mean_A": float(run[cellA_col].mean()),
                        "mean_B": float(run[cellB_col].mean()),
                        "block_index_list": run["index"].tolist() if keep_block_lists else None,
                    })

        for i, row in sub.iterrows():
            if not row["is_sig"]:
                flush_run()
                run_idx, run_dir, prev_end = [], None, row[end_col]
                continue

            if run_idx:
                close_enough = (row[start_col] <= prev_end) or (row[start_col] - prev_end <= max_gap_bp)
                if (row["direction"] == run_dir) and close_enough:
                    run_idx.append(i)
                else:
                    flush_run()
                    run_idx = [i]
                    run_dir = row["direction"]
            else:
                run_idx = [i]
                run_dir = row["direction"]

            prev_end = row[end_col]

        flush_run()

    dmrs_df = pd.DataFrame(dmrs)
    if not dmrs_df.empty and not keep_block_lists:
        dmrs_df = dmrs_df.drop(columns=["block_index_list"], errors="ignore")

    # sort for readability
    if not dmrs_df.empty:
        dmrs_df = dmrs_df.sort_values([chrom_col, start_col, end_col]).reset_index(drop=True)

    return dmrs_df, blocks
