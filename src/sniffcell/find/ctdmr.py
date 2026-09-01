from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from tqdm import tqdm


def means_from_mapping(df: pd.DataFrame, mapping: Dict[str, List[str]]) -> Dict[str, pd.Series]:
    """
    df: rows = loci, cols = sample IDs (exactly as in mapping lists)
    mapping: {"Neuron": [sample1, sample2, ...], "Oligodendrocyte": [...], ...}
    returns: {"Neuron": mean_series, "Oligodendrocyte": mean_series, ...}
    """
    out: Dict[str, pd.Series] = {}
    for group, cols in mapping.items():
        keep = [c for c in cols if c in df.columns]
        if keep:
            out[group] = df[keep].mean(axis=1)
    return out


def _groups_to_mask(group_idx: Sequence[int]) -> int:
    mask = 0
    for idx in group_idx:
        mask |= (1 << int(idx))
    return mask


def _mask_to_groups(mask: int, groups: Sequence[str]) -> list[str]:
    return [groups[i] for i in range(len(groups)) if (mask & (1 << i))]


def _group_label(groups: Sequence[str]) -> str:
    return "+".join(groups)


def _enumerate_unique_bipartitions(groups: Sequence[str]) -> list[tuple[np.ndarray, np.ndarray, int, int]]:
    """
    Enumerate unique non-empty bipartitions by fixing group-0 on the left side.
    This avoids duplicate mirror pairs (A|BCD vs BCD|A).
    """
    n_groups = len(groups)
    if n_groups < 2:
        return []

    full_mask = (1 << n_groups) - 1
    out: list[tuple[np.ndarray, np.ndarray, int, int]] = []
    for mask in range(1, full_mask):
        if mask == full_mask:
            continue
        if (mask & 1) == 0:
            continue
        left_idx = [i for i in range(n_groups) if (mask & (1 << i))]
        right_idx = [i for i in range(n_groups) if i not in left_idx]
        if not left_idx or not right_idx:
            continue
        left_mask = _groups_to_mask(left_idx)
        right_mask = full_mask ^ left_mask
        out.append(
            (
                np.asarray(left_idx, dtype=np.int64),
                np.asarray(right_idx, dtype=np.int64),
                left_mask,
                right_mask,
            )
        )
    return out


def _target_side_from_partition(
    *,
    left_idx: np.ndarray,
    right_idx: np.ndarray,
    left_mask: int,
    right_mask: int,
    delta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Resolve target side per row:
    - Prefer the smaller side (more specific cell-type set).
    - Tie-break by taking the hyper side.
    Returns:
      target_mask, hyper_mask, target_is_hyper (all length n_rows)
    """
    hyper_is_left = delta >= 0
    hyper_mask = np.where(hyper_is_left, left_mask, right_mask)

    if len(left_idx) < len(right_idx):
        target_mask = np.full(delta.shape[0], left_mask, dtype=np.int64)
        target_is_hyper = hyper_is_left
        return target_mask, hyper_mask.astype(np.int64), target_is_hyper

    if len(right_idx) < len(left_idx):
        target_mask = np.full(delta.shape[0], right_mask, dtype=np.int64)
        target_is_hyper = ~hyper_is_left
        return target_mask, hyper_mask.astype(np.int64), target_is_hyper

    target_mask = hyper_mask.astype(np.int64)
    target_is_hyper = np.ones(delta.shape[0], dtype=bool)
    return target_mask, hyper_mask.astype(np.int64), target_is_hyper


def _write_igv_bed9(dmr_df: pd.DataFrame, bed_out: str | Path) -> None:
    out_path = Path(bed_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    igv_cols = [
        "chr",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
    ]
    with out_path.open("w", encoding="utf-8") as f:
        f.write("#" + "\t".join(igv_cols) + "\n")

    if dmr_df.empty:
        return

    bed9 = dmr_df.copy()
    bed9["thickStart"] = bed9["start"]
    bed9["thickEnd"] = bed9["end"]
    bed9["itemRgb"] = 0
    bed9[igv_cols].to_csv(out_path, sep="\t", index=False, header=False, mode="a")


def call_ct_combination_dmrs(
    *,
    idx_df: pd.DataFrame,
    mean_by_group: Dict[str, pd.Series],
    diff_threshold: float = 0.40,
    min_rows: int = 2,
    min_cpgs: int = 3,
    min_bp: int = 0,
    direction: str = "both",
    max_gap_bp: int = 500,
    bed_out: Optional[str] = None,
    max_partitions: int = 4096,
) -> pd.DataFrame:
    """
    Find ctDMRs by evaluating bipartitions of cell types:
      e.g. A vs B|C|D, or A|B vs C|D.

    For each locus, we choose the partition with the largest absolute methylation
    margin and keep rows that pass `diff_threshold`.
    """
    required = {"chr", "start", "end", "startCpG", "endCpG"}
    if not required.issubset(idx_df.columns):
        missing = sorted(required - set(idx_df.columns))
        raise ValueError(f"idx_df missing required columns: {missing}")

    groups = list(mean_by_group.keys())
    if len(groups) < 2:
        raise ValueError("Need >= 2 groups in mean_by_group")

    n_rows = len(idx_df)
    for group, series in mean_by_group.items():
        if len(series) != n_rows:
            raise ValueError(f"Length mismatch for group '{group}': {len(series)} != {n_rows}")

    M = np.column_stack([mean_by_group[g].to_numpy() for g in groups])  # shape (N, G)
    n_groups = len(groups)
    full_mask = (1 << n_groups) - 1

    partitions = _enumerate_unique_bipartitions(groups)
    if not partitions:
        raise ValueError("No valid group bipartitions found (need >= 2 groups).")
    if len(partitions) > max_partitions:
        raise ValueError(
            f"Too many combinations ({len(partitions)}). "
            f"Reduce group count or raise max_partitions (current={max_partitions})."
        )

    best_score = np.full(n_rows, -np.inf, dtype=float)
    second_score = np.full(n_rows, -np.inf, dtype=float)
    winner_partition = np.full(n_rows, -1, dtype=np.int64)
    winner_target_mask = np.zeros(n_rows, dtype=np.int64)
    winner_hyper_mask = np.zeros(n_rows, dtype=np.int64)
    winner_target_is_hyper = np.zeros(n_rows, dtype=bool)

    for part_i, (left_idx, right_idx, left_mask, right_mask) in enumerate(partitions):
        left_mean = M[:, left_idx].mean(axis=1)
        right_mean = M[:, right_idx].mean(axis=1)
        delta = left_mean - right_mean
        margin_mag = np.abs(delta)

        target_mask, hyper_mask, target_is_hyper = _target_side_from_partition(
            left_idx=left_idx,
            right_idx=right_idx,
            left_mask=left_mask,
            right_mask=right_mask,
            delta=delta,
        )

        if direction == "both":
            pass_mask = margin_mag >= diff_threshold
        elif direction == "hyper":
            pass_mask = (margin_mag >= diff_threshold) & target_is_hyper
        elif direction == "hypo":
            pass_mask = (margin_mag >= diff_threshold) & (~target_is_hyper)
        else:
            raise ValueError("direction must be one of {'both','hyper','hypo'}")

        candidate_score = np.where(pass_mask, margin_mag, -np.inf)
        better = candidate_score > best_score
        if not np.any(better):
            continue

        second_score = np.where(
            better,
            best_score,
            np.maximum(second_score, candidate_score),
        )
        best_score = np.where(better, candidate_score, best_score)
        winner_partition[better] = part_i
        winner_target_mask[better] = target_mask[better]
        winner_hyper_mask[better] = hyper_mask[better]
        winner_target_is_hyper[better] = target_is_hyper[better]

    chr_arr = idx_df["chr"].to_numpy()
    start_arr = idx_df["start"].to_numpy()
    end_arr = idx_df["end"].to_numpy()
    scpg_arr = idx_df["startCpG"].to_numpy()
    ecpg_arr = idx_df["endCpG"].to_numpy()

    regions: list[tuple[int, int]] = []
    s_idx: Optional[int] = None
    last_i: Optional[int] = None

    def compatible(i: int, j: int) -> bool:
        if chr_arr[i] != chr_arr[j]:
            return False
        if winner_partition[i] != winner_partition[j]:
            return False
        if winner_target_mask[i] != winner_target_mask[j]:
            return False
        if winner_target_is_hyper[i] != winner_target_is_hyper[j]:
            return False
        if max_gap_bp < 0:
            return True
        return (start_arr[i] - end_arr[j]) <= max_gap_bp

    for i in range(n_rows):
        if winner_partition[i] >= 0:
            if s_idx is None:
                s_idx = i
            elif last_i is not None and not compatible(i, last_i):
                if (last_i - s_idx + 1) >= min_rows:
                    regions.append((s_idx, last_i))
                s_idx = i
        else:
            if s_idx is not None and last_i is not None and (last_i - s_idx + 1) >= min_rows:
                regions.append((s_idx, last_i))
            s_idx = None
        last_i = i

    if s_idx is not None and last_i is not None and (last_i - s_idx + 1) >= min_rows:
        regions.append((s_idx, last_i))

    cols = [
        "chr",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "n_rows",
        "n_cpgs",
        "bp_len",
        "best_group",
        "other_group",
        "best_dir",
        "mean_margin",
        "second_best_margin",
        "rest_std_mean",
        "mean_best_value",
        "mean_rest_value",
        "best_group_leaves",
        "other_group_leaves",
        "hyper_group_leaves",
        "hypo_group_leaves",
        "code_order",
    ] + [f"mean_{g}" for g in groups]

    if not regions:
        empty_df = pd.DataFrame(columns=cols)
        if bed_out:
            _write_igv_bed9(empty_df, bed_out)
            empty_df.to_csv(str(bed_out) + ".full.tsv", sep="\t", index=False)
        return empty_df

    out_rows = []
    for s, e in tqdm(regions, desc="Processing regions"):
        n_rows_region = int(e - s + 1)
        n_cpgs_region = int(ecpg_arr[e] - scpg_arr[s])
        bp_len_region = int(end_arr[e] - start_arr[s])
        if n_rows_region < min_rows or n_cpgs_region < min_cpgs or bp_len_region < min_bp:
            continue

        target_mask = int(winner_target_mask[s])
        hyper_mask = int(winner_hyper_mask[s])
        target_groups = _mask_to_groups(target_mask, groups)
        other_groups = _mask_to_groups(full_mask ^ target_mask, groups)
        hyper_groups = _mask_to_groups(hyper_mask, groups)
        hypo_groups = _mask_to_groups(full_mask ^ hyper_mask, groups)

        if not target_groups or not other_groups:
            continue

        target_idx = [i for i in range(n_groups) if (target_mask & (1 << i))]
        other_idx = [i for i in range(n_groups) if (full_mask ^ target_mask) & (1 << i)]
        region_M = M[s : e + 1, :]

        target_region_mean = region_M[:, target_idx].mean(axis=1)
        other_region_mean = region_M[:, other_idx].mean(axis=1)
        margin_mag = np.abs(target_region_mean - other_region_mean)
        mean_margin = float(np.mean(margin_mag))

        second_region = second_score[s : e + 1]
        finite_second = second_region[np.isfinite(second_region)]
        second_best_margin = float(np.mean(finite_second)) if finite_second.size else 0.0

        rest_std = (
            np.std(region_M[:, other_idx], axis=1, ddof=0)
            if len(other_idx) >= 2
            else np.zeros(region_M.shape[0], dtype=float)
        )
        rest_std_mean = float(np.mean(rest_std))

        per_group_means = [float(np.mean(region_M[:, i])) for i in range(n_groups)]
        mean_best_value = float(np.mean(region_M[:, target_idx]))
        mean_rest_value = float(np.mean(region_M[:, other_idx]))

        best_dir = "hyper" if bool(winner_target_is_hyper[s]) else "hypo"
        best_group = _group_label(target_groups)
        other_group = _group_label(other_groups)
        name = f"DMR_{best_group}_{best_dir}_vs_{other_group}"
        score = int(np.clip(round(mean_margin * 1000), 0, 1000))

        out_rows.append(
            [
                str(chr_arr[s]),
                int(start_arr[s]),
                int(end_arr[e]),
                name,
                score,
                ".",
                n_rows_region,
                n_cpgs_region,
                bp_len_region,
                best_group,
                other_group,
                best_dir,
                mean_margin,
                second_best_margin,
                rest_std_mean,
                mean_best_value,
                mean_rest_value,
                "|".join(target_groups),
                "|".join(other_groups),
                "|".join(hyper_groups),
                "|".join(hypo_groups),
                "|".join(groups),
                *per_group_means,
            ]
        )

    dmr_df = pd.DataFrame(out_rows, columns=cols)
    if not dmr_df.empty:
        dmr_df.sort_values(["chr", "start", "end", "best_group"], inplace=True, ignore_index=True)

    if bed_out:
        _write_igv_bed9(dmr_df, bed_out)
        dmr_df.to_csv(str(bed_out) + ".full.tsv", sep="\t", index=False)

    return dmr_df


def call_ct_specific_dmrs(
    *,
    idx_df: pd.DataFrame,
    mean_by_group: Dict[str, pd.Series],
    diff_threshold: float = 0.40,
    rest_std_threshold: float = 0.10,  # kept for backward compatibility
    min_rows: int = 2,
    min_cpgs: int = 3,
    min_bp: int = 0,
    direction: str = "both",
    max_gap_bp: int = 500,
    bed_out: Optional[str] = "DMR_celltype_specific.bed",
    per_group_bed_prefix: Optional[str] = None,  # kept for backward compatibility
) -> pd.DataFrame:
    """
    Backward-compatible wrapper.
    ctDMRs are now called with combinatorial bipartitions of cell types.
    """
    _ = rest_std_threshold
    _ = per_group_bed_prefix
    return call_ct_combination_dmrs(
        idx_df=idx_df,
        mean_by_group=mean_by_group,
        diff_threshold=diff_threshold,
        min_rows=min_rows,
        min_cpgs=min_cpgs,
        min_bp=min_bp,
        direction=direction,
        max_gap_bp=max_gap_bp,
        bed_out=bed_out,
    )
