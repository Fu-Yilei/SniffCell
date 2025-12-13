import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm
def means_from_mapping(df: pd.DataFrame, mapping: Dict[str, List[str]]) -> Dict[str, pd.Series]:
    """
    df: rows = loci, cols = sample IDs (exactly as in mapping lists)
    mapping: {"Neuron": [sample1, sample2, ...], "Oligodendrocyte": [...], ...}
    returns: {"Neuron": mean_series, "Oligodendrocyte": mean_series, ...}
    """
    out = {}
    for group, cols in mapping.items():
        keep = [c for c in cols if c in df.columns]
        if keep:  # skip empty groups cleanly
            out[group] = df[keep].mean(axis=1)
    return out

def call_ct_specific_dmrs(
    *,
    idx_df: pd.DataFrame,                       # columns: chr,start,end,startCpG,endCpG
    mean_by_group: Dict[str, pd.Series],        # {"Neuron": sA, "Oligo": sB, ...} same length/order as idx_df
    diff_threshold: float = 0.40,               # margin vs rest
    rest_std_threshold: float = 0.10,           # require "rest" groups to be similar (low std)
    min_rows: int = 2,
    min_cpgs: int = 3,
    min_bp: int = 0,
    direction: str = "both",                    # {'both','hyper','hypo'}
    max_gap_bp: int = 500,
    bed_out: str = "DMR_celltype_specific.bed",
    per_group_bed_prefix: Optional[str] = None  # writes <prefix>.<group>.bed if set
) -> pd.DataFrame:
    """
    Multi-class cell-type–specific DMRs with a harmonized-rest constraint.
    Outputs a BED file with a comment-style header (#...) as the first line,
    including mean_best_value and mean_rest_value columns.
    """
    required = {"chr","start","end","startCpG","endCpG"}
    assert required.issubset(idx_df.columns), f"idx_df missing: {required - set(idx_df.columns)}"
    groups = list(mean_by_group.keys())
    assert len(groups) >= 2, "Need ≥2 groups"
    n = len(idx_df)
    for g, s in mean_by_group.items():
        assert len(s) == n, f"Length mismatch for group {g}"

    # (N,G) matrix of per-bin means (column order = groups)
    M = np.vstack([mean_by_group[g].to_numpy() for g in groups]).T
    G = len(groups)

    # One-vs-rest margins
    hyper_margin = np.empty((n, G), float)
    hypo_margin  = np.empty((n, G), float)
    rest_std     = np.empty((n, G), float)

    for gi in range(G):
        rest = np.delete(M, gi, axis=1)
        rest_max = np.max(rest, axis=1)
        rest_min = np.min(rest, axis=1)
        hyper_margin[:, gi] = M[:, gi] - rest_max
        hypo_margin[:, gi]  = rest_min - M[:, gi]
        rest_std[:, gi]     = np.std(rest, axis=1, ddof=0)

    # Signed effect + pass mask with rest-consistency
    if direction == "hyper":
        signed_margin = hyper_margin
        pass_mask = (signed_margin >= diff_threshold) & (rest_std <= rest_std_threshold)
    elif direction == "hypo":
        signed_margin = -hypo_margin
        pass_mask = (hypo_margin >= diff_threshold) & (rest_std <= rest_std_threshold)
    elif direction == "both":
        signed_margin = np.where(
            np.abs(hyper_margin) >= np.abs(hypo_margin),
            hyper_margin,
            -hypo_margin
        )
        pass_mask = (np.abs(signed_margin) >= diff_threshold) & (rest_std <= rest_std_threshold)
    else:
        raise ValueError("direction must be one of {'both','hyper','hypo'}")

    # Winner per bin
    winner = np.full(n, -1, int)
    winner_sign = np.zeros(n, int)
    for i in range(n):
        passing = np.where(pass_mask[i])[0]
        if passing.size == 0:
            continue
        j = passing[np.argmax(np.abs(signed_margin[i, passing]))]
        winner[i] = j
        winner_sign[i] = 1 if signed_margin[i, j] >= 0 else -1

    chr_arr   = idx_df["chr"].to_numpy()
    start_arr = idx_df["start"].to_numpy()
    end_arr   = idx_df["end"].to_numpy()
    scpg_arr  = idx_df["startCpG"].to_numpy()
    ecpg_arr  = idx_df["endCpG"].to_numpy()

    # Merge adjacent bins
    regions: List[Tuple[int,int]] = []
    s = None
    last_i = None
    def compatible(i, j):
        if chr_arr[i] != chr_arr[j]: return False
        if winner[i]  != winner[j]:  return False
        if winner_sign[i] != winner_sign[j]: return False
        return (start_arr[i] - end_arr[j]) <= max_gap_bp if max_gap_bp >= 0 else True

    for i in range(n):
        if winner[i] >= 0:
            if s is None:
                s = i
            elif not compatible(i, last_i):
                if last_i is not None and (last_i - s + 1) >= min_rows:
                    regions.append((s, last_i))
                s = i
        else:
            if s is not None and (last_i - s + 1) >= min_rows:
                regions.append((s, last_i))
            s = None
        last_i = i
    if s is not None and (last_i - s + 1) >= min_rows:
        regions.append((s, last_i))

    cols = [
        "chr","start","end","name","score","strand",
        "n_rows","n_cpgs","bp_len",
        "best_group","best_dir","mean_margin","second_best_margin",
        "rest_std_mean",
        "mean_best_value","mean_rest_value",
    ] + [f"mean_{g}" for g in groups]

    if not regions:
        empty_df = pd.DataFrame(columns=cols)
        with open(bed_out, "w") as f:
            f.write("#" + "\t".join(cols) + "\n")
        empty_df.to_csv(bed_out, sep="\t", header=False, index=False, mode="a")
        if per_group_bed_prefix:
            for g in groups:
                with open(f"{per_group_bed_prefix}.{g}.bed", "w") as f:
                    f.write("#" + "\t".join(cols) + "\n")
        return empty_df

    # Summarize regions
    out_rows = []
    means_np = {g: mean_by_group[g].to_numpy() for g in groups}

    for s_idx, e_idx in tqdm(regions, desc="Processing regions"):
        r_chr, r_start, r_end = chr_arr[s_idx], int(start_arr[s_idx]), int(end_arr[e_idx])
        n_rows_ = int(e_idx - s_idx + 1)
        n_cpgs_ = int(ecpg_arr[e_idx] - scpg_arr[s_idx])
        bp_len_ = int(r_end - r_start)
        if (n_rows_ < min_rows) or (n_cpgs_ < min_cpgs) or (bp_len_ < min_bp):
            continue

        gi = winner[s_idx]
        sign = winner_sign[s_idx]
        region_scores = signed_margin[s_idx:e_idx+1, gi]
        mean_signed_margin = float(np.mean(region_scores))
        mean_margin_mag = abs(mean_signed_margin)

        per_group_scores = np.array([float(np.mean(signed_margin[s_idx:e_idx+1, k])) for k in range(G)])
        tmp_abs = np.abs(per_group_scores); tmp_abs[gi] = -np.inf
        second_best = float(np.max(tmp_abs))
        rest_std_region = float(np.mean(rest_std[s_idx:e_idx+1, gi]))

        per_group_means = [float(np.mean(means_np[g][s_idx:e_idx+1])) for g in groups]
        mean_best_value = per_group_means[gi]
        rest_values = [per_group_means[k] for k in range(G) if k != gi]
        mean_rest_value = float(np.mean(rest_values)) if rest_values else np.nan

        best_group = groups[gi]
        best_dir = "hyper" if sign > 0 else "hypo"
        name = f"DMR_{best_group}_{best_dir}_vs_rest"
        score = int(np.clip(round(mean_margin_mag * 1000), 0, 1000))

        out_rows.append([
            r_chr, r_start, r_end, name, score, ".",
            n_rows_, n_cpgs_, bp_len_,
            best_group, best_dir, mean_margin_mag, second_best,
            rest_std_region,
            mean_best_value, mean_rest_value,
            *per_group_means
        ])
        
    dmr_df = pd.DataFrame(out_rows, columns=cols)
    dmr_df.sort_values(["chr","start","end","best_group"], inplace=True, ignore_index=True)

    # Optional: save full table with all columns for downstream analysis
    full_tsv = bed_out + ".full.tsv"
    dmr_df.to_csv(full_tsv, sep="\t", index=False)

    # ---- IGV-friendly BED9 output ----
    bed9 = dmr_df.copy()
    bed9["thickStart"] = bed9["start"]
    bed9["thickEnd"] = bed9["end"]
    bed9["itemRgb"] = 0  # or "0,0,0" if you want explicit RGB

    igv_cols = [
        "chr", "start", "end",
        "name", "score", "strand",
        "thickStart", "thickEnd", "itemRgb",
    ]

    with open(bed_out, "w") as f:
        # optional comment header
        f.write("#" + "\t".join(igv_cols) + "\n")
    bed9[igv_cols].to_csv(bed_out, sep="\t", header=False, index=False, mode="a")

    # Per-group IGV BEDs
    if per_group_bed_prefix:
        for g in groups:
            sub = bed9[bed9["best_group"] == g]
            out_path = f"{per_group_bed_prefix}.{g}.bed"
            with open(out_path, "w") as f:
                f.write("#" + "\t".join(igv_cols) + "\n")
            sub[igv_cols].to_csv(out_path, sep="\t", header=False, index=False, mode="a")

    return dmr_df
