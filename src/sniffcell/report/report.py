from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd

from . import full_report


LOGGER = logging.getLogger("sniffcell_lite.report")


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "code": "string",
            "majority_code": "string",
            "assigned_code": "string",
            "intersection_code": "string",
        },
    )


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _load_lite_manifest(anno_output: Path) -> tuple[dict[str, Any], Path]:
    for name in ("anno_batch_manifest.json", "anno_compact_manifest.json"):
        path = anno_output / name
        if path.exists():
            return _load_json(path), path
    raise FileNotFoundError(
        f"Could not find anno_batch_manifest.json or anno_compact_manifest.json under: {anno_output}"
    )


def _load_excluded_variant_ids(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"variant exclusion file does not exist: {path}")

    try:
        table = pd.read_csv(path, sep=None, engine="python", dtype=str)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError):
        table = pd.DataFrame()

    if "id" in table.columns:
        values = table["id"].dropna().astype(str).tolist()
    else:
        values = []
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                value = re.split(r"[\t,]", line.strip(), maxsplit=1)[0].strip()
                if value and value.lower() != "id":
                    values.append(value)

    excluded = list(dict.fromkeys(value.strip() for value in values if value.strip()))
    if not excluded:
        raise ValueError(f"variant exclusion file contains no variant IDs: {path}")
    return excluded


def _parse_variant_location(value: Any) -> tuple[str, int, int]:
    text = "" if value is None or pd.isna(value) else str(value).strip()
    match = re.fullmatch(r"([^:\s]+):([0-9,]+)(?:-([0-9,]+))?", text)
    if not match:
        return "", 0, 1
    chrom = match.group(1)
    start_1 = int(match.group(2).replace(",", ""))
    end_1 = int(match.group(3).replace(",", "")) if match.group(3) else start_1
    return chrom, max(0, start_1 - 1), max(start_1, end_1)


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
    except TypeError:
        pass
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _parse_supporting_reads(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        values = value
    else:
        text = "" if value is None or pd.isna(value) else str(value).strip()
        if not text or text in {".", "NA", "None", "null", "[]"}:
            return []
        try:
            parsed = json.loads(text)
        except Exception:
            parsed = None
        values = parsed if isinstance(parsed, list) else re.split(r"[,|;\s]+", text)
    out = []
    for item in values:
        text = str(item).strip()
        if "::" in text:
            text = text.split("::", 1)[1]
        if text:
            out.append(text)
    return list(dict.fromkeys(out))


def _raw_read_assignment_table(path: Path) -> pd.DataFrame:
    reads = _read_table(path)
    if reads.empty:
        return reads
    reads = reads.copy()
    first_col = reads.columns[0]
    if first_col != "readname":
        reads = reads.rename(columns={first_col: "readname"})
    reads["readname"] = reads["readname"].astype(str).map(lambda x: x.split("::", 1)[1] if "::" in x else x)
    return reads


def _load_batch_table(manifest: dict[str, Any]) -> pd.DataFrame:
    inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}
    batch_text = str(inputs.get("batch", "")).strip()
    if not batch_text:
        return pd.DataFrame()
    batch_path = Path(batch_text)
    if not batch_path.exists():
        return pd.DataFrame()
    sep = "," if batch_path.suffix.lower() == ".csv" else "\t"
    batch = pd.read_csv(batch_path, sep=sep)
    if "variant_name" not in batch.columns:
        return pd.DataFrame()
    batch = batch.copy()
    if "variant_location" in batch.columns:
        parsed = batch["variant_location"].map(_parse_variant_location)
        batch["batch_chrom"] = [x[0] for x in parsed]
        batch["batch_start"] = [x[1] for x in parsed]
        batch["batch_end"] = [x[2] for x in parsed]
    return batch


def _with_batch_metadata(assignments: pd.DataFrame, batch: pd.DataFrame) -> pd.DataFrame:
    if batch.empty or "variant_name" not in batch.columns:
        return assignments.copy()
    keep = [
        col
        for col in [
            "variant_name",
            "variant_location",
            "catalog",
            "bam",
            "reference",
            "donor",
            "tissue_code",
            "tissue_name",
            "gcc",
            "platform",
            "batch_chrom",
            "batch_start",
            "batch_end",
        ]
        if col in batch.columns
    ]
    return assignments.merge(batch[keep].drop_duplicates("variant_name"), left_on="id", right_on="variant_name", how="left")


def _support_reads_by_variant(mappings: pd.DataFrame, assignments: pd.DataFrame) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    if not mappings.empty and {"variant_id", "readname"}.issubset(mappings.columns):
        for variant_id, sub in mappings.groupby(mappings["variant_id"].astype(str), sort=False):
            reads = [str(v).strip() for v in sub["readname"].dropna().tolist() if str(v).strip()]
            out[str(variant_id)] = list(dict.fromkeys(reads))
    for row in assignments.to_dict(orient="records"):
        variant_id = str(row.get("id", "")).strip()
        if variant_id and variant_id not in out:
            out[variant_id] = _parse_supporting_reads(row.get("group_a_read_names", "[]"))
    return out


def _normalize_assignment_for_full_report(assignments: pd.DataFrame) -> pd.DataFrame:
    out = assignments.copy()
    if "variant_subtype" not in out.columns:
        out["variant_subtype"] = "VAR"
    out["variant_subtype"] = out["variant_subtype"].fillna("").astype(str)
    out.loc[out["variant_subtype"].str.strip().eq(""), "variant_subtype"] = "VAR"
    if "sv_type" not in out.columns:
        out["sv_type"] = out["variant_subtype"]
    out["sv_type"] = out["sv_type"].fillna("").astype(str)
    out.loc[out["sv_type"].str.strip().eq(""), "sv_type"] = out.loc[out["sv_type"].str.strip().eq(""), "variant_subtype"]
    if "variant_class" not in out.columns:
        out["variant_class"] = "VAR"
    out["variant_class"] = out["variant_class"].fillna("").astype(str)
    out.loc[out["variant_class"].str.strip().eq(""), "variant_class"] = "VAR"
    if "sv_len" not in out.columns:
        out["sv_len"] = pd.NA
    out["sv_len"] = pd.to_numeric(out["sv_len"], errors="coerce")
    if "change_size_bp" in out.columns:
        change_size = pd.to_numeric(out["change_size_bp"], errors="coerce").abs()
        out["sv_len"] = out["sv_len"].where(out["sv_len"].notna(), change_size)
    out["sv_len"] = out["sv_len"].fillna(1)
    if "vaf" not in out.columns:
        out["vaf"] = pd.NA
    out["vaf"] = pd.to_numeric(out["vaf"], errors="coerce")
    return out


def _build_variant_table(assignments: pd.DataFrame, support_reads: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for row in assignments.to_dict(orient="records"):
        variant_id = str(row.get("id", "")).strip()
        chrom = str(row.get("sv_chr", "")).strip()
        start = _safe_int(row.get("sv_pos"), 0) - 1
        if start < 0:
            start = _safe_int(row.get("batch_start"), 0)
        if not chrom:
            chrom = str(row.get("batch_chrom", "")).strip()
        if not chrom and "variant_location" in row:
            chrom, start, parsed_end = _parse_variant_location(row.get("variant_location"))
        else:
            parsed_end = _safe_int(row.get("batch_end"), start + 1)
        size = abs(_safe_int(row.get("sv_len", row.get("change_size_bp", 1)), 1))
        if size <= 0:
            size = max(1, _safe_int(row.get("change_size_bp"), 1))
        end = max(start + 1, parsed_end, start + size)
        reads = support_reads.get(variant_id, [])
        rows.append(
            {
                "variant_id": variant_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "variant_class": str(row.get("variant_class", "SV") or "SV"),
                "variant_subtype": str(row.get("variant_subtype", row.get("sv_type", "VAR")) or "VAR"),
                "change_size_bp": max(1, end - start),
                "group_a_read_names": json.dumps(reads),
                "group_b_read_names": json.dumps([]),
                "supporting_reads": json.dumps(reads),
            }
        )
    return pd.DataFrame(rows)


def _build_runtime_table(assignments: pd.DataFrame, manifest: dict[str, Any]) -> pd.DataFrame:
    inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}
    default_reference = str(inputs.get("reference", "")).strip()
    rows = []
    for row in assignments.to_dict(orient="records"):
        rows.append(
            {
                "id": str(row.get("id", "")).strip(),
                "bam": str(row.get("bam", "")).strip(),
                "bed": str(row.get("catalog", "")).strip(),
                "reference": str(row.get("reference", "")).strip() or default_reference,
            }
        )
    return pd.DataFrame(rows)


def _prepare_full_report_anno_output(
    lite_anno_output: Path,
    excluded_variant_ids: list[str] | None = None,
    exclusion_source: Path | None = None,
    included_variant_ids: list[str] | None = None,
    inclusion_source: Path | None = None,
) -> Path:
    manifest, manifest_path = _load_lite_manifest(lite_anno_output)
    excluded_variant_ids = excluded_variant_ids or []
    included_variant_ids = included_variant_ids or []
    excluded_set = set(excluded_variant_ids)
    included_set = set(included_variant_ids)
    if excluded_set and included_set:
        raise ValueError("include and exclude variant lists cannot be used together")
    if included_set:
        inclusion_digest = hashlib.sha256(
            "\n".join(sorted(included_set)).encode("utf-8")
        ).hexdigest()[:12]
        compat_dir = lite_anno_output / f".sniffcell_full_report_compat_include_{inclusion_digest}"
    elif excluded_set:
        exclusion_digest = hashlib.sha256(
            "\n".join(sorted(excluded_set)).encode("utf-8")
        ).hexdigest()[:12]
        compat_dir = lite_anno_output / f".sniffcell_full_report_compat_exclude_{exclusion_digest}"
    else:
        compat_dir = lite_anno_output / ".sniffcell_full_report_compat"
    compat_dir.mkdir(parents=True, exist_ok=True)

    assignments = _read_table(lite_anno_output / "variant_assignment.tsv")
    if assignments.empty:
        raise ValueError(f"No assignments found in {lite_anno_output / 'variant_assignment.tsv'}")
    if "id" not in assignments.columns:
        raise ValueError("variant_assignment.tsv must contain an id column")
    assignment_ids = assignments["id"].fillna("").astype(str).str.strip()
    matched_excluded_ids = sorted(excluded_set.intersection(assignment_ids))
    missing_excluded_ids = sorted(excluded_set.difference(assignment_ids))
    matched_included_ids = sorted(included_set.intersection(assignment_ids))
    missing_included_ids = sorted(included_set.difference(assignment_ids))
    assignments_before_exclusion = len(assignments)
    if included_set:
        assignments = assignments.loc[assignment_ids.isin(included_set)].copy()
        if assignments.empty:
            raise ValueError("none of the requested variant IDs matched the assignments")
    elif excluded_set:
        assignments = assignments.loc[~assignment_ids.isin(excluded_set)].copy()
        if assignments.empty:
            raise ValueError("variant exclusion removed every assignment; no report can be generated")
    assignments = _with_batch_metadata(assignments, _load_batch_table(manifest))
    assignments = _normalize_assignment_for_full_report(assignments)
    mappings = _read_table(lite_anno_output / "support_read_mappings.tsv")
    support_reads = _support_reads_by_variant(mappings, assignments)

    variant_table = _build_variant_table(assignments, support_reads)
    runtime_table = _build_runtime_table(assignments, manifest)
    raw_reads = _raw_read_assignment_table(lite_anno_output / "reads_classification.tsv")

    variant_assignment_path = compat_dir / "variant_assignment.tsv"
    sv_assignment_path = compat_dir / "sv_assignment.tsv"
    variant_table_path = compat_dir / "lite_variants.tsv"
    runtime_table_path = compat_dir / "lite_variant_runtime.tsv"
    reads_path = compat_dir / "reads_classification.tsv"
    excluded_variants_path = compat_dir / "excluded_variants.tsv"
    included_variants_path = compat_dir / "included_variants.tsv"

    assignments.to_csv(variant_assignment_path, sep="\t", index=False)
    assignments.to_csv(sv_assignment_path, sep="\t", index=False)
    variant_table.to_csv(variant_table_path, sep="\t", index=False)
    runtime_table.to_csv(runtime_table_path, sep="\t", index=False)
    raw_reads.to_csv(reads_path, sep="\t", index=False)
    if excluded_set:
        pd.DataFrame({"id": sorted(excluded_set)}).to_csv(excluded_variants_path, sep="\t", index=False)
    elif excluded_variants_path.exists():
        excluded_variants_path.unlink()
    if included_set:
        pd.DataFrame({"id": sorted(included_set)}).to_csv(included_variants_path, sep="\t", index=False)
    elif included_variants_path.exists():
        included_variants_path.unlink()

    first_runtime = runtime_table.iloc[0].to_dict() if not runtime_table.empty else {}
    full_manifest = {
        "command": "sniffcell-lite report adapter",
        "source": "sniffcell-lite",
        "inputs": {
            "lite_anno_output": str(lite_anno_output.resolve()),
            "lite_manifest": str(manifest_path.resolve()),
            "vcf": str(variant_table_path.resolve()),
            "variants": str(variant_table_path.resolve()),
            "bam": str(first_runtime.get("bam", "")),
            "bed": str(first_runtime.get("bed", "")),
            "reference": str(first_runtime.get("reference", "")),
            "exclude_variants": (
                str(exclusion_source.resolve()) if exclusion_source is not None else ""
            ),
            "include_variants": (
                str(inclusion_source.resolve()) if inclusion_source is not None else ""
            ),
        },
        "runtime": {
            "window": int(manifest.get("runtime", {}).get("window", 10000)) if isinstance(manifest.get("runtime", {}), dict) else 10000,
            "exclusions": {
                "requested_unique": len(excluded_set),
                "matched": len(matched_excluded_ids),
                "missing": len(missing_excluded_ids),
                "assignments_before": assignments_before_exclusion,
                "assignments_after": len(assignments),
                "missing_ids": missing_excluded_ids,
            },
            "inclusions": {
                "requested_unique": len(included_set),
                "matched": len(matched_included_ids),
                "missing": len(missing_included_ids),
                "assignments_before": assignments_before_exclusion,
                "assignments_after": len(assignments),
                "missing_ids": missing_included_ids,
            },
        },
        "outputs": {
            "variant_assignment": str(variant_assignment_path.resolve()),
            "sv_assignment": str(sv_assignment_path.resolve()),
            "reads_classification": str(reads_path.resolve()),
            "lite_variant_runtime": str(runtime_table_path.resolve()),
            "excluded_variants": (
                str(excluded_variants_path.resolve()) if excluded_set else ""
            ),
            "included_variants": (
                str(included_variants_path.resolve()) if included_set else ""
            ),
        },
    }
    (compat_dir / "anno_run_manifest.json").write_text(json.dumps(full_manifest, indent=2, sort_keys=True), encoding="utf-8")
    return compat_dir


def _add_exclusion_provenance_to_report(
    args,
    compat_anno_output: Path,
    exclusion_source: Path,
) -> None:
    report_dir, _, _, archive_path = full_report._resolve_report_paths(
        compat_anno_output, getattr(args, "output", None)
    )
    source_manifest = _load_json(compat_anno_output / "anno_run_manifest.json")
    exclusions = source_manifest.get("runtime", {}).get("exclusions", {})
    excluded_variants = compat_anno_output / "excluded_variants.tsv"
    report_excluded_variants = report_dir / "excluded_variants.tsv"
    report_excluded_variants.write_bytes(excluded_variants.read_bytes())

    report_manifest_path = report_dir / "report_manifest.json"
    report_manifest = _load_json(report_manifest_path)
    report_manifest["exclusions"] = {
        **exclusions,
        "source": str(exclusion_source.resolve()),
        "normalized_ids": str(report_excluded_variants.resolve()),
    }
    report_manifest.setdefault("outputs", {})["excluded_variants_tsv"] = str(
        report_excluded_variants.resolve()
    )
    report_manifest_path.write_text(
        json.dumps(report_manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    if archive_path is not None:
        full_report._write_report_archive(report_dir, archive_path)


def _add_inclusion_provenance_to_report(
    args,
    compat_anno_output: Path,
    inclusion_source: Path,
) -> None:
    report_dir, _, _, archive_path = full_report._resolve_report_paths(
        compat_anno_output, getattr(args, "output", None)
    )
    source_manifest = _load_json(compat_anno_output / "anno_run_manifest.json")
    inclusions = source_manifest.get("runtime", {}).get("inclusions", {})
    included_variants = compat_anno_output / "included_variants.tsv"
    report_included_variants = report_dir / "included_variants.tsv"
    report_included_variants.write_bytes(included_variants.read_bytes())

    report_manifest_path = report_dir / "report_manifest.json"
    report_manifest = _load_json(report_manifest_path)
    report_manifest["inclusions"] = {
        **inclusions,
        "source": str(inclusion_source.resolve()),
        "normalized_ids": str(report_included_variants.resolve()),
    }
    report_manifest.setdefault("outputs", {})["included_variants_tsv"] = str(
        report_included_variants.resolve()
    )
    report_manifest_path.write_text(
        json.dumps(report_manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    if archive_path is not None:
        full_report._write_report_archive(report_dir, archive_path)


def _full_report_args(args, compat_anno_output: Path) -> SimpleNamespace:
    max_sv = getattr(args, "max_sv", None)
    if max_sv is None:
        max_sv = getattr(args, "max_variants", 0)
    fmt = getattr(args, "format", None) or getattr(args, "figure_format", "png")
    return SimpleNamespace(
        anno_output=str(compat_anno_output),
        output=getattr(args, "output", None),
        min_overlap_pct=float(getattr(args, "min_overlap_pct", 0.0)),
        overlap_filter_mode=str(getattr(args, "overlap_filter_mode", "hard_clip")),
        overlap_gradient_exponent=float(getattr(args, "overlap_gradient_exponent", 0.5)),
        min_majority_pct=float(getattr(args, "min_majority_pct", 0.0)),
        include_unassigned=bool(getattr(args, "include_unassigned", False)),
        allow_hard_conflict=bool(getattr(args, "allow_hard_conflict", False)),
        max_sv=int(max_sv),
        with_figures=bool(getattr(args, "with_figures", False)),
        window=int(getattr(args, "window", 5000)),
        max_reads=int(getattr(args, "max_reads", 250)),
        format=str(fmt),
        figure_profile=str(getattr(args, "figure_profile", "full")),
        figure_dpi=int(getattr(args, "figure_dpi", 160)),
        reuse_existing_viz=bool(getattr(args, "reuse_existing_viz", False)),
        figure_threads=int(getattr(args, "figure_threads", 1)),
        with_igvviz=bool(getattr(args, "with_igvviz", False)),
        igv_bams=getattr(args, "igv_bams", None),
        igv_cmd=str(getattr(args, "igv_cmd", "igv.sh")),
        igv_snapshot_format=str(getattr(args, "igv_snapshot_format", "png")),
        igv_snapshot_width=int(getattr(args, "igv_snapshot_width", 3600)),
        igv_snapshot_height=int(getattr(args, "igv_snapshot_height", 1600)),
        reuse_existing_igvviz=bool(getattr(args, "reuse_existing_igvviz", False)),
        with_igvreport=bool(getattr(args, "with_igvreport", False)),
    )


def report_main(args) -> None:
    if not LOGGER.handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    lite_anno_output = Path(args.anno_output)
    if not lite_anno_output.exists():
        raise FileNotFoundError(f"anno output does not exist: {lite_anno_output}")

    exclusion_source_text = str(getattr(args, "exclude_variants", "") or "").strip()
    exclusion_source = Path(exclusion_source_text) if exclusion_source_text else None
    inclusion_source_text = str(getattr(args, "include_variants", "") or "").strip()
    inclusion_source = Path(inclusion_source_text) if inclusion_source_text else None
    if exclusion_source is not None and inclusion_source is not None:
        raise ValueError("--exclude_variants and --include_variants cannot be used together")
    excluded_variant_ids = (
        _load_excluded_variant_ids(exclusion_source) if exclusion_source is not None else []
    )
    included_variant_ids = (
        _load_excluded_variant_ids(inclusion_source) if inclusion_source is not None else []
    )
    compat_anno_output = _prepare_full_report_anno_output(
        lite_anno_output,
        excluded_variant_ids=excluded_variant_ids,
        exclusion_source=exclusion_source,
        included_variant_ids=included_variant_ids,
        inclusion_source=inclusion_source,
    )
    LOGGER.info("Prepared SniffCell report-compatible lite anno folder: %s", compat_anno_output)
    full_report.report_main(_full_report_args(args, compat_anno_output))
    if exclusion_source is not None:
        _add_exclusion_provenance_to_report(args, compat_anno_output, exclusion_source)
    if inclusion_source is not None:
        _add_inclusion_provenance_to_report(args, compat_anno_output, inclusion_source)
