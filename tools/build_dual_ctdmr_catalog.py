#!/usr/bin/env python3
"""Build a channel-aware brain ctDMR catalog around an immutable legacy backbone."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


COORD_COLS = ["chr", "start", "end"]


def read_ctdmr(path: str, modification: str, prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing = set(COORD_COLS + ["mean_Neuron", "mean_Oligodendrocyte"]) - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    frame = frame.copy()
    frame["chr"] = frame["chr"].astype(str)
    frame["start"] = pd.to_numeric(frame["start"], errors="raise").astype(np.int64)
    frame["end"] = pd.to_numeric(frame["end"], errors="raise").astype(np.int64)
    if np.any(frame["end"].to_numpy() <= frame["start"].to_numpy()):
        raise ValueError(f"{path} contains an interval with end <= start")
    frame["modification"] = modification
    frame["evidence_id"] = [f"{prefix}_{i + 1:07d}" for i in range(len(frame))]
    frame["neuron_minus_oligo"] = (
        pd.to_numeric(frame["mean_Neuron"], errors="coerce")
        - pd.to_numeric(frame["mean_Oligodendrocyte"], errors="coerce")
    )
    frame["effect_direction"] = np.where(
        frame["neuron_minus_oligo"] >= 0,
        "Neuron_hyper",
        "Oligodendrocyte_hyper",
    )
    frame["source_ctdmr"] = str(Path(path).resolve())
    return frame


def overlap_pairs(left: pd.DataFrame, right: pd.DataFrame) -> list[tuple[int, int, int]]:
    """Return positional row pairs and overlap bp for half-open BED intervals."""
    pairs: list[tuple[int, int, int]] = []
    right_by_chrom = {
        chrom: group.sort_values(["start", "end"], kind="stable")
        for chrom, group in right.groupby("chr", sort=False)
    }
    for chrom, left_chrom in left.groupby("chr", sort=False):
        right_chrom = right_by_chrom.get(chrom)
        if right_chrom is None or right_chrom.empty:
            continue
        left_sorted = left_chrom.sort_values(["start", "end"], kind="stable")
        right_rows = list(
            zip(
                right_chrom.index.to_numpy(dtype=np.int64),
                right_chrom["start"].to_numpy(dtype=np.int64),
                right_chrom["end"].to_numpy(dtype=np.int64),
                strict=False,
            )
        )
        add_at = 0
        active: list[tuple[int, int, int]] = []
        for left_idx, left_start, left_end in zip(
            left_sorted.index.to_numpy(dtype=np.int64),
            left_sorted["start"].to_numpy(dtype=np.int64),
            left_sorted["end"].to_numpy(dtype=np.int64),
            strict=False,
        ):
            while add_at < len(right_rows) and right_rows[add_at][1] < left_end:
                active.append(right_rows[add_at])
                add_at += 1
            active = [row for row in active if row[2] > left_start]
            for right_idx, right_start, right_end in active:
                overlap_bp = min(left_end, right_end) - max(left_start, right_start)
                if overlap_bp > 0:
                    pairs.append((int(left_idx), int(right_idx), int(overlap_bp)))
    return pairs


def add_overlap_annotations(
    target: pd.DataFrame,
    other: pd.DataFrame,
    pairs: list[tuple[int, int, int]],
    label: str,
) -> None:
    counts = np.zeros(len(target), dtype=np.int32)
    same = np.zeros(len(target), dtype=np.int32)
    opposite = np.zeros(len(target), dtype=np.int32)
    strongest = np.full(len(target), np.nan, dtype=np.float32)
    target_direction = target["effect_direction"].to_numpy()
    other_direction = other["effect_direction"].to_numpy()
    other_effect = np.abs(other["neuron_minus_oligo"].to_numpy(dtype=float))
    for target_idx, other_idx, _ in pairs:
        counts[target_idx] += 1
        if target_direction[target_idx] == other_direction[other_idx]:
            same[target_idx] += 1
        else:
            opposite[target_idx] += 1
        value = other_effect[other_idx]
        if np.isfinite(value) and (not np.isfinite(strongest[target_idx]) or value > strongest[target_idx]):
            strongest[target_idx] = value
    target[f"overlap_{label}_count"] = counts
    target[f"overlap_{label}_same_direction"] = same
    target[f"overlap_{label}_opposite_direction"] = opposite
    target[f"overlap_{label}_strongest_abs_effect"] = strongest


def links_frame(
    left: pd.DataFrame,
    right: pd.DataFrame,
    pairs: list[tuple[int, int, int]],
    relation: str,
) -> pd.DataFrame:
    rows = []
    for left_idx, right_idx, overlap_bp in pairs:
        left_row = left.iloc[left_idx]
        right_row = right.iloc[right_idx]
        rows.append(
            {
                "relation": relation,
                "left_evidence_id": left_row["evidence_id"],
                "right_evidence_id": right_row["evidence_id"],
                "chr": left_row["chr"],
                "left_start": int(left_row["start"]),
                "left_end": int(left_row["end"]),
                "right_start": int(right_row["start"]),
                "right_end": int(right_row["end"]),
                "overlap_bp": overlap_bp,
                "left_modification": left_row["modification"],
                "right_modification": right_row["modification"],
                "left_direction": left_row["effect_direction"],
                "right_direction": right_row["effect_direction"],
                "direction_relation": (
                    "same" if left_row["effect_direction"] == right_row["effect_direction"] else "opposite"
                ),
            }
        )
    return pd.DataFrame(rows)


def interval_count(frame: pd.DataFrame, column: str) -> int:
    return int((frame[column].to_numpy() > 0).sum())


def merged_component_count(frame: pd.DataFrame) -> int:
    """Count genomic components after BED-style merging of overlapping/book-ended rows."""
    count = 0
    for _, chrom in frame.groupby("chr", sort=False):
        ordered = chrom.sort_values(["start", "end"], kind="stable")
        current_end: int | None = None
        for start, end in zip(ordered["start"], ordered["end"], strict=False):
            if current_end is None or int(start) > current_end:
                count += 1
                current_end = int(end)
            else:
                current_end = max(current_end, int(end))
    return count


def summarize_effects(frame: pd.DataFrame, label: str) -> dict[str, object]:
    effect = np.abs(frame["neuron_minus_oligo"].to_numpy(dtype=float))
    return {
        "set": label,
        "ctdmrs": len(frame),
        "median_abs_neuron_minus_oligo": float(np.nanmedian(effect)),
        "q10_abs_neuron_minus_oligo": float(np.nanquantile(effect, 0.10)),
        "q90_abs_neuron_minus_oligo": float(np.nanquantile(effect, 0.90)),
        "median_cpgs": float(np.nanmedian(pd.to_numeric(frame.get("n_cpgs"), errors="coerce"))),
        "median_bp": float(np.nanmedian(frame["end"].to_numpy() - frame["start"].to_numpy())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", required=True)
    parser.add_argument("--five-mc", required=True)
    parser.add_argument("--five-hmc", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    legacy = read_ctdmr(args.legacy, "modifiedC", "legacy")
    five_mc = read_ctdmr(args.five_mc, "5mC", "5mc")
    five_hmc = read_ctdmr(args.five_hmc, "5hmC", "5hmc")

    legacy_mc_pairs = overlap_pairs(legacy, five_mc)
    legacy_hmc_pairs = overlap_pairs(legacy, five_hmc)
    add_overlap_annotations(legacy, five_mc, legacy_mc_pairs, "5mC")
    add_overlap_annotations(legacy, five_hmc, legacy_hmc_pairs, "5hmC")

    five_mc["overlaps_legacy"] = False
    five_hmc["overlaps_legacy"] = False
    for legacy_idx, ont_idx, _ in legacy_mc_pairs:
        five_mc.at[ont_idx, "overlaps_legacy"] = True
    for legacy_idx, ont_idx, _ in legacy_hmc_pairs:
        five_hmc.at[ont_idx, "overlaps_legacy"] = True

    novel_mc = five_mc.loc[~five_mc["overlaps_legacy"]].copy().reset_index(drop=True)
    novel_hmc = five_hmc.loc[~five_hmc["overlaps_legacy"]].copy().reset_index(drop=True)
    novel_cross_pairs = overlap_pairs(novel_mc, novel_hmc)
    add_overlap_annotations(novel_mc, novel_hmc, novel_cross_pairs, "novel_5hmC")
    reverse_cross_pairs = [(right_idx, left_idx, bp) for left_idx, right_idx, bp in novel_cross_pairs]
    add_overlap_annotations(novel_hmc, novel_mc, reverse_cross_pairs, "novel_5mC")

    legacy["catalog_role"] = "legacy_backbone"
    legacy["overlaps_legacy"] = True
    five_mc["catalog_role"] = np.where(five_mc["overlaps_legacy"], "legacy_annotation", "novel_addition")
    five_hmc["catalog_role"] = np.where(five_hmc["overlaps_legacy"], "legacy_annotation", "novel_addition")
    novel_mc["catalog_role"] = "novel_addition"
    novel_hmc["catalog_role"] = "novel_addition"

    all_evidence = pd.concat([legacy, five_mc, five_hmc], ignore_index=True, sort=False)
    dual_catalog = pd.concat([legacy, novel_mc, novel_hmc], ignore_index=True, sort=False)
    additional = pd.concat([novel_mc, novel_hmc], ignore_index=True, sort=False)

    link_frames = [
        links_frame(legacy, five_mc, legacy_mc_pairs, "legacy_to_5mC"),
        links_frame(legacy, five_hmc, legacy_hmc_pairs, "legacy_to_5hmC"),
        links_frame(novel_mc, novel_hmc, novel_cross_pairs, "novel_5mC_to_5hmC"),
    ]
    overlap_links = pd.concat(link_frames, ignore_index=True, sort=False)

    legacy.to_csv(output_dir / "brain_cereb.legacy.annotated.tsv.gz", sep="\t", index=False)
    additional.to_csv(output_dir / "brain_cereb.additional_ont_ctdmr.tsv.gz", sep="\t", index=False)
    dual_catalog.to_csv(output_dir / "brain_cereb.dual_catalog.tsv.gz", sep="\t", index=False)
    all_evidence.to_csv(output_dir / "brain_cereb.all_channel_evidence.tsv.gz", sep="\t", index=False)
    overlap_links.to_csv(output_dir / "brain_cereb.overlap_links.tsv.gz", sep="\t", index=False)
    dual_catalog[["chr", "start", "end", "evidence_id", "modification", "effect_direction"]].to_csv(
        output_dir / "brain_cereb.dual_catalog.bed",
        sep="\t",
        index=False,
        header=False,
    )

    summary_rows = [
        summarize_effects(legacy, "legacy_modifiedC_backbone"),
        summarize_effects(five_mc, "ONT_5mC_all"),
        summarize_effects(five_hmc, "ONT_5hmC_all"),
        summarize_effects(novel_mc, "ONT_5mC_novel_to_legacy"),
        summarize_effects(novel_hmc, "ONT_5hmC_novel_to_legacy"),
    ]
    pd.DataFrame(summary_rows).to_csv(output_dir / "brain_cereb.effect_quality_summary.tsv", sep="\t", index=False)

    cross_same = sum(
        novel_mc.iloc[left_idx]["effect_direction"] == novel_hmc.iloc[right_idx]["effect_direction"]
        for left_idx, right_idx, _ in novel_cross_pairs
    )
    summary = {
        "policy": {
            "legacy_backbone": "preserved exactly; never replaced or resegmented",
            "ONT_overlap_with_legacy": "retained as channel-specific annotation, not added as duplicate catalog rows",
            "ONT_novel": "added as channel-specific evidence; 5mC and 5hmC are never summed or averaged",
            "interval_overlap": "at least one base under BED half-open coordinates",
        },
        "counts": {
            "legacy_backbone": len(legacy),
            "5mC_all": len(five_mc),
            "5hmC_all": len(five_hmc),
            "5mC_overlapping_legacy": int(five_mc["overlaps_legacy"].sum()),
            "5hmC_overlapping_legacy": int(five_hmc["overlaps_legacy"].sum()),
            "5mC_novel_to_legacy": len(novel_mc),
            "5hmC_novel_to_legacy": len(novel_hmc),
            "legacy_supported_by_5mC": interval_count(legacy, "overlap_5mC_count"),
            "legacy_supported_by_5hmC": interval_count(legacy, "overlap_5hmC_count"),
            "legacy_supported_by_both": int(
                ((legacy["overlap_5mC_count"] > 0) & (legacy["overlap_5hmC_count"] > 0)).sum()
            ),
            "novel_5mC_overlapping_novel_5hmC": interval_count(novel_mc, "overlap_novel_5hmC_count"),
            "novel_5hmC_overlapping_novel_5mC": interval_count(novel_hmc, "overlap_novel_5mC_count"),
            "novel_cross_channel_overlap_links": len(novel_cross_pairs),
            "novel_cross_channel_same_direction_links": int(cross_same),
            "novel_cross_channel_opposite_direction_links": int(len(novel_cross_pairs) - cross_same),
            "dual_catalog_rows": len(dual_catalog),
            "dual_catalog_genomic_components": merged_component_count(dual_catalog),
            "additional_ont_genomic_components": merged_component_count(additional),
            "all_evidence_rows": len(all_evidence),
        },
        "inputs": {
            "legacy": str(Path(args.legacy).resolve()),
            "5mC": str(Path(args.five_mc).resolve()),
            "5hmC": str(Path(args.five_hmc).resolve()),
        },
    }
    with (output_dir / "brain_cereb.dual_catalog_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    print(json.dumps(summary["counts"], indent=2))


if __name__ == "__main__":
    main()
