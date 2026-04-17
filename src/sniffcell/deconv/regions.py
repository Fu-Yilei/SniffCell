from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pysam


@dataclass(frozen=True)
class TargetRegion:
    chrom: str
    start: int
    end: int
    name: str | None = None


@dataclass(frozen=True)
class ResolvedRegionalInputs:
    regions_arg: str
    selected_ctdmr_count: int
    subset_bed_path: str
    subset_bam_path: str
    targets_bed_path: str
    expanded_bed_path: str
    manifest_path: str
    target_regions: tuple[TargetRegion, ...]
    expanded_regions: tuple[TargetRegion, ...]
    subset_bam_read_count: int


def _normalize_chrom(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    lowered = text.lower()
    return lowered[3:] if lowered.startswith("chr") else lowered


def parse_region_spec(value: str) -> TargetRegion:
    text = str(value).strip().replace(",", "")
    match = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", text)
    if match is None:
        raise ValueError(f"Invalid region string: {value!r}")
    chrom, start_text, end_text = match.groups()
    start = int(start_text)
    end = int(end_text)
    if start < 0 or end < 0:
        raise ValueError("Region coordinates must be non-negative")
    if end <= start:
        raise ValueError("Region end must be greater than start")
    return TargetRegion(chrom=chrom, start=start, end=end, name=None)


def load_target_regions(regions_arg: str) -> list[TargetRegion]:
    text = str(regions_arg).strip()
    if not text:
        raise ValueError("--regions must not be empty")

    candidate = Path(text).expanduser()
    if candidate.exists():
        targets: list[TargetRegion] = []
        with candidate.open(encoding="utf-8") as handle:
            for line_no, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    raise ValueError(
                        f"BED row {line_no} in {candidate} must contain at least 3 tab-delimited columns"
                    )
                chrom = fields[0].strip()
                start = int(fields[1])
                end = int(fields[2])
                if end <= start:
                    raise ValueError(f"BED row {line_no} in {candidate} has end <= start")
                name = fields[3].strip() if len(fields) >= 4 and fields[3].strip() else None
                targets.append(TargetRegion(chrom=chrom, start=start, end=end, name=name))
        if not targets:
            raise ValueError(f"No usable target regions were found in BED: {candidate}")
        return targets

    return [parse_region_spec(text)]


def select_ctdmrs_for_target(
    ctdmr_df: pd.DataFrame,
    *,
    target: TargetRegion,
    left_ctdmrs: int,
    right_ctdmrs: int,
) -> pd.DataFrame:
    same_chrom = ctdmr_df.loc[
        ctdmr_df["chr"].map(_normalize_chrom).eq(_normalize_chrom(target.chrom))
    ].copy()
    if same_chrom.empty:
        raise ValueError(f"No ctDMRs found on chromosome {target.chrom}")

    same_chrom = same_chrom.sort_values(["start", "end"], kind="stable").reset_index(drop=True)
    overlaps = same_chrom.loc[
        same_chrom["start"].lt(target.end) & same_chrom["end"].gt(target.start)
    ].copy()
    left = same_chrom.loc[same_chrom["end"].le(target.start)].copy()
    right = same_chrom.loc[same_chrom["start"].ge(target.end)].copy()

    left_selected = left.iloc[0:0].copy()
    right_selected = right.iloc[0:0].copy()
    if left_ctdmrs > 0 and not left.empty:
        left = left.assign(_distance=target.start - left["end"])
        left_selected = left.nsmallest(left_ctdmrs, "_distance").drop(columns="_distance")
    if right_ctdmrs > 0 and not right.empty:
        right = right.assign(_distance=right["start"] - target.end)
        right_selected = right.nsmallest(right_ctdmrs, "_distance").drop(columns="_distance")

    selected = (
        pd.concat([overlaps, left_selected, right_selected], ignore_index=True)
        .drop_duplicates()
        .sort_values(["start", "end"], kind="stable")
        .reset_index(drop=True)
    )
    if selected.empty:
        raise ValueError(
            f"No ctDMRs were selected near {target.chrom}:{target.start}-{target.end}. "
            "Increase the region flank ctDMR counts or verify the ctDMR catalog."
        )
    return selected


def expand_target_from_ctdmrs(selected_ctdmrs: pd.DataFrame) -> TargetRegion:
    if selected_ctdmrs.empty:
        raise ValueError("selected_ctdmrs must not be empty")
    chrom = str(selected_ctdmrs.iloc[0]["chr"])
    return TargetRegion(
        chrom=chrom,
        start=int(selected_ctdmrs["start"].min()),
        end=int(selected_ctdmrs["end"].max()),
        name=None,
    )


def _merge_regions(regions: list[TargetRegion]) -> list[TargetRegion]:
    if not regions:
        return []

    merged: list[TargetRegion] = []
    for region in sorted(regions, key=lambda row: (_normalize_chrom(row.chrom), row.start, row.end)):
        if not merged:
            merged.append(region)
            continue
        prev = merged[-1]
        if _normalize_chrom(prev.chrom) == _normalize_chrom(region.chrom) and region.start <= prev.end:
            merged[-1] = TargetRegion(
                chrom=prev.chrom,
                start=prev.start,
                end=max(prev.end, region.end),
                name=prev.name,
            )
            continue
        merged.append(region)
    return merged


def _write_regions_bed(path: Path, regions: list[TargetRegion]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for idx, region in enumerate(regions, start=1):
            name = region.name or f"region_{idx}"
            handle.write(f"{region.chrom}\t{region.start}\t{region.end}\t{name}\n")


def _write_subset_bam(
    *,
    input_bam: str,
    output_bam: Path,
    regions: list[TargetRegion],
) -> int:
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    n_reads = 0
    with pysam.AlignmentFile(str(input_bam), "rb") as in_bam:
        with pysam.AlignmentFile(str(output_bam), "wb", template=in_bam) as out_bam:
            for region in regions:
                for read in in_bam.fetch(str(region.chrom), int(region.start), int(region.end)):
                    if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                        continue
                    if int(read.reference_end) <= int(region.start) or int(read.reference_start) >= int(region.end):
                        continue
                    out_bam.write(read)
                    n_reads += 1
    pysam.index(str(output_bam))
    return int(n_reads)


def resolve_regional_inputs(
    *,
    ctdmr_df: pd.DataFrame,
    input_bam: str,
    output_dir: str,
    regions_arg: str,
    left_ctdmrs: int,
    right_ctdmrs: int,
) -> ResolvedRegionalInputs:
    targets = load_target_regions(regions_arg)
    selected_frames: list[pd.DataFrame] = []
    expanded_regions: list[TargetRegion] = []
    for target in targets:
        selected = select_ctdmrs_for_target(
            ctdmr_df,
            target=target,
            left_ctdmrs=int(left_ctdmrs),
            right_ctdmrs=int(right_ctdmrs),
        )
        selected_frames.append(selected)
        expanded_regions.append(expand_target_from_ctdmrs(selected))

    selected_ctdmrs = (
        pd.concat(selected_frames, ignore_index=True)
        .drop_duplicates()
        .sort_values(["chr", "start", "end"], kind="stable")
        .reset_index(drop=True)
    )
    merged_expanded = _merge_regions(expanded_regions)

    region_dir = Path(output_dir) / "deconv_regions"
    region_dir.mkdir(parents=True, exist_ok=True)
    targets_bed = region_dir / "targets.bed"
    expanded_bed = region_dir / "expanded_regions.bed"
    subset_bed = region_dir / "ctdmr_subset.tsv"
    subset_bam = region_dir / "region_subset.bam"
    manifest = region_dir / "region_manifest.json"

    _write_regions_bed(targets_bed, targets)
    _write_regions_bed(expanded_bed, merged_expanded)
    selected_ctdmrs.to_csv(subset_bed, sep="\t", index=False)
    subset_bam_read_count = _write_subset_bam(
        input_bam=input_bam,
        output_bam=subset_bam,
        regions=merged_expanded,
    )

    payload = {
        "regions_arg": str(regions_arg),
        "left_ctdmrs": int(left_ctdmrs),
        "right_ctdmrs": int(right_ctdmrs),
        "target_regions": [
            {"chrom": row.chrom, "start": int(row.start), "end": int(row.end), "name": row.name}
            for row in targets
        ],
        "expanded_regions": [
            {"chrom": row.chrom, "start": int(row.start), "end": int(row.end), "name": row.name}
            for row in merged_expanded
        ],
        "selected_ctdmr_count": int(len(selected_ctdmrs)),
        "subset_bam_read_count": int(subset_bam_read_count),
        "subset_bed": str(subset_bed.resolve()),
        "subset_bam": str(subset_bam.resolve()),
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return ResolvedRegionalInputs(
        regions_arg=str(regions_arg),
        selected_ctdmr_count=int(len(selected_ctdmrs)),
        subset_bed_path=str(subset_bed.resolve()),
        subset_bam_path=str(subset_bam.resolve()),
        targets_bed_path=str(targets_bed.resolve()),
        expanded_bed_path=str(expanded_bed.resolve()),
        manifest_path=str(manifest.resolve()),
        target_regions=tuple(targets),
        expanded_regions=tuple(merged_expanded),
        subset_bam_read_count=int(subset_bam_read_count),
    )
