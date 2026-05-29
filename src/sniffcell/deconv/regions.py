from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pysam

logger = logging.getLogger(__name__)


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


@dataclass(frozen=True)
class ResolvedRegionPlan:
    regions_arg: str
    selected_ctdmr_count: int
    subset_bed_path: str
    subset_regions_bed_path: str
    targets_bed_path: str
    region_summary_path: str
    selected_summary_path: str
    manifest_path: str
    target_regions: tuple[TargetRegion, ...]
    expanded_regions: tuple[TargetRegion, ...]


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


def _total_region_bp(regions: list[TargetRegion]) -> int:
    return int(sum(max(0, int(region.end) - int(region.start)) for region in regions))


def load_ctdmr_table(bed_path: str) -> pd.DataFrame:
    ctdmrs = pd.read_csv(bed_path, sep="\t")
    if not ctdmrs.empty and isinstance(ctdmrs.columns[0], str) and ctdmrs.columns[0].startswith("#"):
        ctdmrs.rename(columns={ctdmrs.columns[0]: ctdmrs.columns[0].lstrip("#")}, inplace=True)
    required = ["chr", "start", "end"]
    missing = [col for col in required if col not in ctdmrs.columns]
    if missing:
        raise ValueError(
            f"ctDMR table missing required columns: {missing}. "
            "Use the main sniffcell find TSV output, not the headerless .igv.bed companion."
        )
    ctdmrs = ctdmrs.drop_duplicates(ignore_index=True)
    ctdmrs = ctdmrs.sort_values(["chr", "start", "end"], kind="stable", ignore_index=True)
    return ctdmrs


def _relation_to_target(ctdmr_row: pd.Series, target: TargetRegion) -> tuple[str, int]:
    start = int(ctdmr_row["start"])
    end = int(ctdmr_row["end"])
    if start < int(target.end) and end > int(target.start):
        return "overlap", 0
    if end <= int(target.start):
        return "left_flank", int(target.start) - end
    return "right_flank", start - int(target.end)


def _join_unique(values: pd.Series) -> str:
    out: list[str] = []
    seen: set[str] = set()
    for value in values.dropna():
        text = str(value).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return "|".join(out)


def _selected_summary_extra_columns(ctdmr_df: pd.DataFrame) -> list[str]:
    preferred = [
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
    ]
    cols = [col for col in preferred if col in ctdmr_df.columns]
    cols.extend(
        col for col in ctdmr_df.columns
        if col.startswith("mean_") and col not in cols
    )
    return cols


def _build_selected_summary(
    *,
    target_index: int,
    target: TargetRegion,
    selected: pd.DataFrame,
) -> pd.DataFrame:
    extra_cols = _selected_summary_extra_columns(selected)
    rows: list[dict[str, object]] = []
    for _, row in selected.iterrows():
        relation, distance_bp = _relation_to_target(row, target)
        out = {
            "target_index": int(target_index),
            "target_name": target.name or f"target_{target_index}",
            "target_chrom": target.chrom,
            "target_start": int(target.start),
            "target_end": int(target.end),
            "ctdmr_chr": row["chr"],
            "ctdmr_start": int(row["start"]),
            "ctdmr_end": int(row["end"]),
            "relation": relation,
            "distance_bp": int(distance_bp),
        }
        for col in extra_cols:
            out[col] = row.get(col, pd.NA)
        rows.append(out)
    return pd.DataFrame(rows)


def _build_region_summary_row(
    *,
    target_index: int,
    target: TargetRegion,
    expanded: TargetRegion,
    selected_summary: pd.DataFrame,
) -> dict[str, object]:
    relation_counts = selected_summary["relation"].value_counts().to_dict()
    distance = selected_summary["distance_bp"]
    nonoverlap_distance = distance[selected_summary["relation"] != "overlap"]
    row: dict[str, object] = {
        "target_index": int(target_index),
        "target_name": target.name or f"target_{target_index}",
        "target_chrom": target.chrom,
        "target_start": int(target.start),
        "target_end": int(target.end),
        "subset_chrom": expanded.chrom,
        "subset_start": int(expanded.start),
        "subset_end": int(expanded.end),
        "subset_bp": int(expanded.end) - int(expanded.start),
        "selected_ctdmr_count": int(len(selected_summary)),
        "overlap_ctdmr_count": int(relation_counts.get("overlap", 0)),
        "left_flank_ctdmr_count": int(relation_counts.get("left_flank", 0)),
        "right_flank_ctdmr_count": int(relation_counts.get("right_flank", 0)),
        "nearest_flank_distance_bp": (
            int(nonoverlap_distance.min()) if not nonoverlap_distance.empty else pd.NA
        ),
    }
    for col in ("best_group", "other_group", "best_dir", "code_order"):
        if col in selected_summary.columns:
            row[col + "s"] = _join_unique(selected_summary[col])
    return row


def _regions_to_manifest_rows(regions: list[TargetRegion]) -> list[dict[str, object]]:
    return [
        {"chrom": row.chrom, "start": int(row.start), "end": int(row.end), "name": row.name}
        for row in regions
    ]


def resolve_region_plan(
    *,
    ctdmr_df: pd.DataFrame,
    output_dir: str,
    regions_arg: str,
    left_ctdmrs: int,
    right_ctdmrs: int,
    ctdmr_path: str | None = None,
    subset_regions_filename: str = "subset_regions.bed",
) -> ResolvedRegionPlan:
    targets = load_target_regions(regions_arg)
    selected_frames: list[pd.DataFrame] = []
    selected_summary_frames: list[pd.DataFrame] = []
    region_summary_rows: list[dict[str, object]] = []
    expanded_regions: list[TargetRegion] = []

    for target_index, target in enumerate(targets, start=1):
        selected = select_ctdmrs_for_target(
            ctdmr_df,
            target=target,
            left_ctdmrs=int(left_ctdmrs),
            right_ctdmrs=int(right_ctdmrs),
        )
        expanded = expand_target_from_ctdmrs(selected)
        selected_summary = _build_selected_summary(
            target_index=target_index,
            target=target,
            selected=selected,
        )

        selected_frames.append(selected)
        selected_summary_frames.append(selected_summary)
        expanded_regions.append(expanded)
        region_summary_rows.append(
            _build_region_summary_row(
                target_index=target_index,
                target=target,
                expanded=expanded,
                selected_summary=selected_summary,
            )
        )

    selected_ctdmrs = (
        pd.concat(selected_frames, ignore_index=True)
        .drop_duplicates()
        .sort_values(["chr", "start", "end"], kind="stable")
        .reset_index(drop=True)
    )
    selected_summary = pd.concat(selected_summary_frames, ignore_index=True)
    region_summary = pd.DataFrame(region_summary_rows)
    merged_expanded = _merge_regions(expanded_regions)

    region_dir = Path(output_dir)
    region_dir.mkdir(parents=True, exist_ok=True)
    targets_bed = region_dir / "targets.bed"
    subset_regions_bed = region_dir / subset_regions_filename
    subset_bed = region_dir / "ctdmr_subset.tsv"
    region_summary_path = region_dir / "ctdmr_region_summary.tsv"
    selected_summary_path = region_dir / "ctdmr_selected_summary.tsv"
    manifest = region_dir / "region_manifest.json"

    _write_regions_bed(targets_bed, targets)
    _write_regions_bed(subset_regions_bed, merged_expanded)
    selected_ctdmrs.to_csv(subset_bed, sep="\t", index=False)
    region_summary.to_csv(region_summary_path, sep="\t", index=False)
    selected_summary.to_csv(selected_summary_path, sep="\t", index=False)

    payload: dict[str, object] = {
        "command": "regions",
        "version": "v1",
        "regions_arg": str(regions_arg),
        "left_ctdmrs": int(left_ctdmrs),
        "right_ctdmrs": int(right_ctdmrs),
        "target_regions": _regions_to_manifest_rows(targets),
        "expanded_regions": _regions_to_manifest_rows(merged_expanded),
        "expanded_region_bp": _total_region_bp(merged_expanded),
        "selected_ctdmr_count": int(len(selected_ctdmrs)),
        "inputs": {
            "ctdmr_bed": str(Path(ctdmr_path).expanduser().resolve()) if ctdmr_path else None,
        },
        "outputs": {
            "subset_bed": str(subset_bed.resolve()),
            "subset_regions_bed": str(subset_regions_bed.resolve()),
            "expanded_bed": str(subset_regions_bed.resolve()),
            "targets_bed": str(targets_bed.resolve()),
            "region_summary": str(region_summary_path.resolve()),
            "selected_summary": str(selected_summary_path.resolve()),
        },
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return ResolvedRegionPlan(
        regions_arg=str(regions_arg),
        selected_ctdmr_count=int(len(selected_ctdmrs)),
        subset_bed_path=str(subset_bed.resolve()),
        subset_regions_bed_path=str(subset_regions_bed.resolve()),
        targets_bed_path=str(targets_bed.resolve()),
        region_summary_path=str(region_summary_path.resolve()),
        selected_summary_path=str(selected_summary_path.resolve()),
        manifest_path=str(manifest.resolve()),
        target_regions=tuple(targets),
        expanded_regions=tuple(merged_expanded),
    )


def _write_subset_bam(
    *,
    input_bam: str,
    output_bam: Path,
    regions: list[TargetRegion],
    regions_bed: Path | None = None,
    threads: int = 1,
) -> int:
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    unsorted_bam = output_bam.with_name(f"{output_bam.stem}.unsorted{output_bam.suffix}")
    for bam_path in (output_bam, unsorted_bam):
        for path in (bam_path, Path(str(bam_path) + ".bai"), bam_path.with_suffix(".bai")):
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    extra_threads = max(0, int(threads) - 1)
    if regions_bed is not None:
        logger.info(
            "Extracting regional BAM with samtools view: regions=%d bp=%d bed=%s input=%s output=%s threads=%d",
            len(regions),
            _total_region_bp(regions),
            regions_bed,
            input_bam,
            unsorted_bam,
            max(1, int(threads)),
        )
        view_args = ["-M", "-b", "-L", str(regions_bed), "-o", str(unsorted_bam)]
        if extra_threads > 0:
            view_args = ["-@", str(extra_threads)] + view_args
        view_args.append(str(input_bam))
        pysam.view(*view_args, catch_stdout=False)
    else:
        with pysam.AlignmentFile(str(input_bam), "rb") as in_bam:
            with pysam.AlignmentFile(str(unsorted_bam), "wb", template=in_bam) as out_bam:
                for region in regions:
                    for read in in_bam.fetch(str(region.chrom), int(region.start), int(region.end)):
                        if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                            continue
                        if int(read.reference_end) <= int(region.start) or int(read.reference_start) >= int(region.end):
                            continue
                        out_bam.write(read)

    sort_args = []
    if extra_threads > 0:
        sort_args.extend(["-@", str(extra_threads)])
    sort_args.extend(["-o", str(output_bam), str(unsorted_bam)])
    logger.info("Sorting regional BAM: input=%s output=%s", unsorted_bam, output_bam)
    pysam.sort(*sort_args)

    index_args = [str(output_bam)]
    if extra_threads > 0:
        index_args = ["-@", str(extra_threads), str(output_bam)]
    logger.info("Indexing regional BAM: %s", output_bam)
    pysam.index(*index_args)

    idxstats = pysam.idxstats(str(output_bam)).splitlines()
    n_reads = 0
    for row in idxstats:
        fields = row.split("\t")
        if len(fields) >= 4:
            n_reads += int(fields[2]) + int(fields[3])

    for path in (unsorted_bam, Path(str(unsorted_bam) + ".bai"), unsorted_bam.with_suffix(".bai")):
        try:
            path.unlink()
        except FileNotFoundError:
            pass
    return int(n_reads)


def resolve_regional_inputs(
    *,
    ctdmr_df: pd.DataFrame,
    input_bam: str,
    output_dir: str,
    regions_arg: str,
    left_ctdmrs: int,
    right_ctdmrs: int,
    threads: int = 1,
) -> ResolvedRegionalInputs:
    region_dir = Path(output_dir) / "deconv_regions"
    plan = resolve_region_plan(
        ctdmr_df=ctdmr_df,
        output_dir=str(region_dir),
        regions_arg=regions_arg,
        left_ctdmrs=left_ctdmrs,
        right_ctdmrs=right_ctdmrs,
        subset_regions_filename="expanded_regions.bed",
    )

    subset_bam = region_dir / "region_subset.bam"
    subset_bam_read_count = _write_subset_bam(
        input_bam=input_bam,
        output_bam=subset_bam,
        regions=list(plan.expanded_regions),
        regions_bed=Path(plan.subset_regions_bed_path),
        threads=threads,
    )

    manifest = Path(plan.manifest_path)
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    outputs = payload.setdefault("outputs", {})
    if isinstance(outputs, dict):
        outputs["subset_bam"] = str(subset_bam.resolve())
    payload["subset_bam_read_count"] = int(subset_bam_read_count)
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return ResolvedRegionalInputs(
        regions_arg=str(regions_arg),
        selected_ctdmr_count=int(plan.selected_ctdmr_count),
        subset_bed_path=plan.subset_bed_path,
        subset_bam_path=str(subset_bam.resolve()),
        targets_bed_path=plan.targets_bed_path,
        expanded_bed_path=plan.subset_regions_bed_path,
        manifest_path=plan.manifest_path,
        target_regions=plan.target_regions,
        expanded_regions=plan.expanded_regions,
        subset_bam_read_count=int(subset_bam_read_count),
    )


def _resolve_flank_counts(args) -> tuple[int, int, int]:
    regions_ctdmrs = int(getattr(args, "regions_ctdmrs", 10))
    left = getattr(args, "regions_left_ctdmrs", None)
    right = getattr(args, "regions_right_ctdmrs", None)
    left_ctdmrs = regions_ctdmrs if left is None else int(left)
    right_ctdmrs = regions_ctdmrs if right is None else int(right)
    if regions_ctdmrs < 0:
        raise ValueError("regions_ctdmrs must be >= 0")
    if left_ctdmrs < 0:
        raise ValueError("regions_left_ctdmrs must be >= 0")
    if right_ctdmrs < 0:
        raise ValueError("regions_right_ctdmrs must be >= 0")
    return regions_ctdmrs, left_ctdmrs, right_ctdmrs


def regions_main(args) -> ResolvedRegionPlan:
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    _regions_ctdmrs, left_ctdmrs, right_ctdmrs = _resolve_flank_counts(args)
    ctdmr_df = load_ctdmr_table(str(args.bed))
    plan = resolve_region_plan(
        ctdmr_df=ctdmr_df,
        output_dir=str(args.output),
        regions_arg=str(args.regions),
        left_ctdmrs=left_ctdmrs,
        right_ctdmrs=right_ctdmrs,
        ctdmr_path=str(args.bed),
        subset_regions_filename="subset_regions.bed",
    )
    logger.info("Wrote target BED: %s", plan.targets_bed_path)
    logger.info("Wrote subset regions BED: %s", plan.subset_regions_bed_path)
    logger.info("Wrote ctDMR subset: %s", plan.subset_bed_path)
    logger.info("Wrote ctDMR summaries: %s, %s", plan.region_summary_path, plan.selected_summary_path)
    logger.info("Selected ctDMRs: %d", plan.selected_ctdmr_count)
    return plan
