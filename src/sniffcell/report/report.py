from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import html
import json
import logging
import re
import shlex
import tarfile
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from sniffcell.anno.vcf_to_df import read_vcf_to_df
from sniffcell.viz import viz as viz_module
from sniffcell.viz import igvviz as igvviz_module
from . import igvreport as igvreport_module

try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover - optional dependency
    tqdm = None


def _has_text(value: object) -> bool:
    if pd.isna(value):
        return False
    return bool(str(value).strip())


def _parse_bool(value: object) -> object:
    if pd.isna(value):
        return pd.NA
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    return pd.NA


def _fmt_float(value: object, *, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "NA"


def _fmt_text(value: object) -> str:
    try:
        if pd.isna(value):
            return "NA"
    except TypeError:
        pass
    text = str(value).strip()
    return text if text else "NA"


def _to_json_scalar(value: object) -> object:
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    return value


def _load_json_dict(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _build_dashboard_records(selected_report_df: pd.DataFrame) -> list[dict[str, object]]:
    if selected_report_df.empty:
        return []

    cols = [
        "id",
        "sv_chr",
        "sv_pos",
        "sv_type",
        "sv_len",
        "vaf",
        "overlap_pct",
        "majority_pct",
        "n_supporting",
        "n_overlapped",
        "primary_celltype",
        "linked_celltypes",
        "has_hard_conflict",
        "viz_status",
    ]
    available_cols = [col for col in cols if col in selected_report_df.columns]
    base_records = selected_report_df[available_cols].to_dict(orient="records")
    records: list[dict[str, object]] = []
    for row in base_records:
        payload = {col: None for col in cols}
        for col, value in row.items():
            payload[col] = _to_json_scalar(value)
        records.append(payload)
    return records


def _safe_slug(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    slug = slug.strip("._")
    return slug or "sv"


def _normalize_review_status(value: object) -> str:
    try:
        if pd.isna(value):
            return "undecided"
    except TypeError:
        pass
    text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if text == "real":
        return "real"
    if text in {"not_real", "notreal"}:
        return "not_real"
    return "undecided"


def _review_state_class(status: str) -> str:
    normalized = _normalize_review_status(status)
    if normalized == "real":
        return "review-state-real"
    if normalized == "not_real":
        return "review-state-not-real"
    return "review-state-undecided"


def _review_state_label(status: str) -> str:
    normalized = _normalize_review_status(status)
    if normalized == "real":
        return "real"
    if normalized == "not_real":
        return "not real"
    return "undecided"


def _load_sv_assignment(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "id" not in df.columns:
        raise ValueError(f"sv_assignment file is missing required column 'id': {path}")

    out = df.copy()
    for col in (
        "assigned_code",
        "linked_celltypes",
        "primary_celltype",
        "majority_code",
    ):
        if col not in out.columns:
            out[col] = ""
        out[col] = out[col].astype("string")

    for col in ("overlap_pct", "majority_pct"):
        if col not in out.columns:
            out[col] = pd.NA
        out[col] = pd.to_numeric(out[col], errors="coerce")

    for col in ("n_supporting", "n_overlapped"):
        if col not in out.columns:
            out[col] = pd.NA
        out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")

    if "sv_len" not in out.columns:
        out["sv_len"] = pd.NA
    out["sv_len"] = pd.to_numeric(out["sv_len"], errors="coerce").astype("Int64")

    if "vaf" not in out.columns:
        out["vaf"] = pd.NA
    out["vaf"] = pd.to_numeric(out["vaf"], errors="coerce")

    if "sv_type" not in out.columns:
        out["sv_type"] = pd.Series("", index=out.index, dtype="string")
    else:
        out["sv_type"] = out["sv_type"].astype("string")

    if "has_hard_conflict" not in out.columns:
        out["has_hard_conflict"] = pd.Series(pd.array([pd.NA] * len(out), dtype="boolean"))
    else:
        out["has_hard_conflict"] = out["has_hard_conflict"].map(_parse_bool).astype("boolean")

    out["id"] = out["id"].astype("string")
    return out


def _backfill_sv_fields_from_manifest_vcf(
    sv_df: pd.DataFrame,
    manifest_payload: dict[str, object],
    logger: logging.Logger,
) -> pd.DataFrame:
    if sv_df.empty or not isinstance(manifest_payload, dict):
        return sv_df

    inputs = manifest_payload.get("inputs", {})
    if not isinstance(inputs, dict):
        return sv_df

    vcf_text = str(inputs.get("vcf", "")).strip()
    if not vcf_text:
        return sv_df

    needs_sv_type = sv_df["sv_type"].fillna("").astype(str).str.strip().eq("").any()
    needs_vaf = sv_df["vaf"].isna().any()
    needs_sv_len = sv_df["sv_len"].isna().any()
    if not (needs_sv_type or needs_vaf or needs_sv_len):
        return sv_df

    vcf_path = Path(vcf_text)
    if not vcf_path.exists():
        logger.warning("Could not backfill report SV fields because manifest VCF does not exist: %s", vcf_path)
        return sv_df

    try:
        vcf_df = read_vcf_to_df(str(vcf_path))
    except Exception as exc:
        logger.warning("Could not backfill report SV fields from manifest VCF %s: %s", vcf_path, exc)
        return sv_df

    if vcf_df.empty or "id" not in vcf_df.columns:
        return sv_df

    vcf_cols = [col for col in ("id", "sv_type", "vaf", "sv_len") if col in vcf_df.columns]
    if len(vcf_cols) <= 1:
        return sv_df

    merged = sv_df.merge(
        vcf_df[vcf_cols].drop_duplicates(subset=["id"]),
        on="id",
        how="left",
        suffixes=("", "__vcf"),
    )

    if "sv_type__vcf" in merged.columns:
        missing_sv_type = merged["sv_type"].fillna("").astype(str).str.strip().eq("")
        merged.loc[missing_sv_type, "sv_type"] = merged.loc[missing_sv_type, "sv_type__vcf"]
        merged["sv_type"] = merged["sv_type"].astype("string")

    for col in ("vaf", "sv_len"):
        extra_col = f"{col}__vcf"
        if extra_col not in merged.columns:
            continue
        current = pd.to_numeric(merged[col], errors="coerce")
        recovered = pd.to_numeric(merged[extra_col], errors="coerce")
        merged[col] = current.where(current.notna(), recovered)

    drop_cols = [col for col in merged.columns if col.endswith("__vcf")]
    if drop_cols:
        merged = merged.drop(columns=drop_cols)

    merged["sv_len"] = pd.to_numeric(merged["sv_len"], errors="coerce").astype("Int64")
    merged["vaf"] = pd.to_numeric(merged["vaf"], errors="coerce")
    merged["sv_type"] = merged["sv_type"].fillna("").astype("string")
    return merged


def _select_high_confidence_svs(
    sv_df: pd.DataFrame,
    *,
    min_overlap_pct: float,
    min_majority_pct: float,
    include_unassigned: bool,
    allow_hard_conflict: bool,
    max_sv: int,
) -> pd.DataFrame:
    if not (0.0 <= float(min_overlap_pct) <= 1.0):
        raise ValueError("min_overlap_pct must be in [0, 1]")
    if not (0.0 <= float(min_majority_pct) <= 1.0):
        raise ValueError("min_majority_pct must be in [0, 1]")
    if int(max_sv) < 0:
        raise ValueError("max_sv must be >= 0")

    selected = sv_df.copy()
    if "assigned_code" not in selected.columns:
        selected["assigned_code"] = pd.Series("", index=selected.index, dtype="string")
    if "linked_celltypes" not in selected.columns:
        selected["linked_celltypes"] = pd.Series("", index=selected.index, dtype="string")
    if "has_hard_conflict" not in selected.columns:
        selected["has_hard_conflict"] = pd.Series(pd.array([pd.NA] * len(selected), dtype="boolean"))
    if "overlap_pct" not in selected.columns:
        selected["overlap_pct"] = pd.Series(pd.NA, index=selected.index, dtype="Float64")
    if "majority_pct" not in selected.columns:
        selected["majority_pct"] = pd.Series(pd.NA, index=selected.index, dtype="Float64")
    if "n_overlapped" not in selected.columns:
        selected["n_overlapped"] = pd.Series(pd.NA, index=selected.index, dtype="Int64")
    if "id" not in selected.columns:
        selected["id"] = pd.Series(pd.NA, index=selected.index, dtype="string")
    linked_mask = selected["linked_celltypes"].map(_has_text)
    linked_mask = pd.Series(linked_mask, index=selected.index).fillna(False).astype(bool)
    selected = selected.loc[linked_mask]

    if not include_unassigned:
        assigned_mask = selected["assigned_code"].map(_has_text)
        assigned_mask = pd.Series(assigned_mask, index=selected.index).fillna(False).astype(bool)
        selected = selected.loc[assigned_mask]
    if not allow_hard_conflict:
        keep_mask = ~selected["has_hard_conflict"].fillna(False).astype(bool)
        selected = selected.loc[keep_mask]

    selected = selected[selected["overlap_pct"].fillna(0.0) >= float(min_overlap_pct)]
    selected = selected[selected["majority_pct"].fillna(0.0) >= float(min_majority_pct)]
    selected = selected.sort_values(
        ["majority_pct", "overlap_pct", "n_overlapped", "id"],
        ascending=[False, False, False, True],
        kind="stable",
        ignore_index=True,
    )

    if int(max_sv) > 0:
        selected = selected.iloc[: int(max_sv)].copy()
    return selected


def _report_dir_from_archive_path(archive_path: Path) -> Path:
    name = archive_path.name
    lowered = name.lower()
    if lowered.endswith(".tar.gz"):
        stem = name[:-7]
    elif lowered.endswith(".tgz"):
        stem = name[:-4]
    elif lowered.endswith(".gz"):
        stem = name[:-3]
    else:
        stem = archive_path.stem
    stem = stem.strip()
    if not stem:
        stem = "report"
    return archive_path.with_name(stem)


def _resolve_report_paths(anno_output: Path, output: str | None) -> tuple[Path, Path, Path, Path | None]:
    archive_path: Path | None = None
    if output is None:
        report_dir = anno_output / "report"
        html_path = report_dir / "index.html"
    else:
        out = Path(output)
        lower_name = out.name.lower()
        if out.suffix.lower() == ".html":
            html_path = out
            report_dir = out.parent if str(out.parent) else Path(".")
        elif lower_name.endswith(".tar.gz") or lower_name.endswith(".tgz") or lower_name.endswith(".gz"):
            archive_path = out
            report_dir = _report_dir_from_archive_path(out)
            html_path = report_dir / "index.html"
        else:
            report_dir = out
            html_path = report_dir / "index.html"
    figure_dir = report_dir / "figures"
    return report_dir, figure_dir, html_path, archive_path


def _write_report_archive(report_dir: Path, archive_path: Path) -> Path:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode="w:gz") as tar:
        tar.add(str(report_dir), arcname=report_dir.name)
    return archive_path


def _build_viz_cli_command(
    *,
    anno_output: Path,
    sv_id: str,
    output_path: Path,
    window: int,
    max_reads: int,
    fmt: str,
    dpi: int,
    exact_window: bool,
    skip_methylation_overlay: bool,
) -> str:
    parts = [
        "sniffcell",
        "viz",
        "--anno_output",
        str(anno_output),
        "-s",
        str(sv_id),
        "-o",
        str(output_path),
        "-f",
        str(fmt),
        "--dpi",
        str(dpi),
        "-w",
        str(window),
        "-m",
        str(max_reads),
    ]
    if exact_window:
        parts.append("--exact_window")
    if skip_methylation_overlay:
        parts.append("--skip_methylation_overlay")
    parts.append("--export_tables")
    return " ".join(shlex.quote(p) for p in parts)


def _render_one_viz_panel(
    *,
    anno_output: Path,
    sv_id: str,
    figure_path: Path,
    window: int,
    max_reads: int,
    fmt: str,
    dpi: int,
    exact_window: bool,
    skip_methylation_overlay: bool,
    reuse_existing_viz: bool,
) -> tuple[str, str]:
    if reuse_existing_viz and figure_path.exists():
        return "reused", ""

    try:
        viz_args = SimpleNamespace(
            anno_output=str(anno_output),
            sv_id=sv_id,
            input=None,
            vcf=None,
            reference=None,
            bed=None,
            read_assignment=None,
            kanpig_read_names=None,
            window=int(window),
            max_reads=int(max_reads),
            format=fmt,
            dpi=int(dpi),
            exact_window=bool(exact_window),
            skip_methylation_overlay=bool(skip_methylation_overlay),
            export_tables=True,
            output=str(figure_path),
        )
        viz_module.viz_main(viz_args)
    except Exception as exc:
        return "failed", str(exc)

    if figure_path.exists():
        return "rendered", ""
    return "failed_no_output", "viz completed but output figure was not found."


def _build_igvviz_cli_command(
    *,
    anno_output: Path,
    sv_id: str,
    output_dir: Path,
    window: int,
    igv_bams: list[str] | None,
    igv_cmd: str,
    snapshot_format: str,
    snapshot_width: int,
    snapshot_height: int,
) -> str:
    parts = [
        "sniffcell",
        "igvviz",
        "--anno_output",
        str(anno_output),
        "-s",
        str(sv_id),
        "-o",
        str(output_dir),
        "-w",
        str(window),
        "--igv_cmd",
        str(igv_cmd),
        "--keep_intermediates",
    ]
    parts += [
        "--snapshot_format",
        str(snapshot_format),
        "--snapshot_width",
        str(snapshot_width),
        "--snapshot_height",
        str(snapshot_height),
    ]
    if igv_bams:
        parts.extend(["-i", *[str(x) for x in igv_bams]])
    return " ".join(shlex.quote(p) for p in parts)


def _read_igvviz_snapshot_stats(manifest_path: Path) -> tuple[int, int]:
    if not manifest_path.exists():
        return 0, 0
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception:
        return 0, 0
    jobs = payload.get("jobs", [])
    if not isinstance(jobs, list):
        return 0, 0
    total = 0
    existing = 0
    for job in jobs:
        if not isinstance(job, dict):
            continue
        snap = str(job.get("snapshot", "")).strip()
        if not snap:
            continue
        total += 1
        if Path(snap).exists():
            existing += 1
    return total, existing


def _load_igvviz_snapshot_rows(manifest_path: Path, html_parent: Path) -> list[dict[str, object]]:
    if not manifest_path.exists():
        return []
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception:
        return []
    jobs = payload.get("jobs", [])
    if not isinstance(jobs, list):
        return []

    out: list[dict[str, object]] = []
    for idx, job in enumerate(jobs, start=1):
        if not isinstance(job, dict):
            continue
        snap_text = str(job.get("snapshot", "")).strip()
        if not snap_text:
            continue
        bam_text = str(job.get("bam", "")).strip()
        bam_label = Path(bam_text).name if bam_text else f"bam_{idx:02d}"
        snap_path = Path(snap_text)
        exists = snap_path.exists()
        snap_rel = ""
        if exists:
            try:
                snap_rel = snap_path.relative_to(html_parent).as_posix()
            except ValueError:
                snap_rel = str(snap_path.resolve())
        out.append(
            {
                "bam": bam_text,
                "bam_label": bam_label,
                "snapshot": str(snap_path),
                "snapshot_rel": snap_rel,
                "exists": bool(exists),
            }
        )
    return out


def _render_one_igvviz_bundle(
    *,
    anno_output: Path,
    sv_id: str,
    output_dir: Path,
    expected_manifest_path: Path,
    window: int,
    igv_bams: list[str] | None,
    igv_cmd: str,
    snapshot_format: str,
    snapshot_width: int,
    snapshot_height: int,
    reuse_existing_igvviz: bool,
) -> tuple[str, str]:
    if reuse_existing_igvviz and expected_manifest_path.exists():
        total, existing = _read_igvviz_snapshot_stats(expected_manifest_path)
        if total > 0 and existing == total:
            return "reused", ""

    try:
        igv_args = SimpleNamespace(
            anno_output=str(anno_output),
            sv_id=sv_id,
            input=(list(igv_bams) if igv_bams else None),
            vcf=None,
            reference=None,
            bed=None,
            kanpig_read_names=None,
            window=int(window),
            visibility_window=int(window),
            phase_tag="HP",
            support_tag="SC",
            igv_cmd=str(igv_cmd),
            batch_only=False,
            keep_intermediates=True,
            snapshot_format=str(snapshot_format),
            snapshot_width=int(snapshot_width),
            snapshot_height=int(snapshot_height),
            output=str(output_dir),
        )
        igvviz_module.igvviz_main(igv_args)
    except Exception as exc:
        return "failed", str(exc)

    if not expected_manifest_path.exists():
        return "failed_no_output", "igvviz completed but output manifest was not found."
    total, existing = _read_igvviz_snapshot_stats(expected_manifest_path)
    if total == 0:
        return "failed_no_output", "igvviz manifest exists but no snapshot jobs were recorded."
    if existing == total:
        return "rendered", ""
    return "rendered_partial", f"Missing {total-existing}/{total} IGV snapshots."


def _build_report_html(
    *,
    generated_at: str,
    anno_output: Path,
    sv_assignment_path: Path,
    review_storage_key: str,
    filters: dict[str, object],
    viz: dict[str, object],
    igvreport: dict[str, object],
    total_sv: int,
    selected_count: int,
    rendered_count: int,
    failed_count: int,
    rows: list[dict[str, object]],
    dashboard_records: list[dict[str, object]],
) -> str:
    page: list[str] = []
    page.append("<!doctype html>")
    page.append("<html lang=\"en\">")
    page.append("<head>")
    page.append("<meta charset=\"utf-8\">")
    page.append("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    page.append("<title>SniffCell SV Report</title>")
    page.append("<style>")
    page.append(
        "body{font-family:Helvetica,Arial,sans-serif;background:#f3f5f7;color:#13212c;margin:0;padding:24px;}"
        ".wrap{max-width:1200px;margin:0 auto;}"
        "h1{margin:0 0 8px 0;font-size:28px;}"
        ".meta{color:#44525f;font-size:14px;margin-bottom:16px;}"
        ".stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:16px 0 20px;}"
        ".card{background:#ffffff;border-radius:10px;padding:12px 14px;box-shadow:0 1px 2px rgba(0,0,0,0.08);}"
        ".label{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:#5d6b77;}"
        ".value{font-size:22px;font-weight:700;color:#0f2233;}"
        ".filter{background:#ffffff;border-radius:10px;padding:10px 14px;margin-bottom:18px;font-size:13px;line-height:1.5;}"
        ".sv{background:#ffffff;border-radius:10px;padding:14px 14px 18px;margin:14px 0;"
        "box-shadow:0 1px 2px rgba(0,0,0,0.08);}"
        ".sv h2{margin:0 0 8px 0;font-size:20px;}"
        ".kv{font-size:14px;color:#1c2e3b;margin:2px 0;}"
        ".err{color:#9b1c1c;font-weight:600;}"
        ".cmd{margin-top:8px;padding:8px;background:#f8fafc;border:1px solid #d8dee4;border-radius:8px;overflow:auto;}"
        ".copy{margin-top:8px;border:0;border-radius:6px;padding:6px 10px;background:#1f6feb;color:#fff;cursor:pointer;font-size:13px;}"
        ".copy:hover{background:#1856b0;}"
        ".dash{background:#ffffff;border-radius:10px;padding:14px 14px 18px;margin:14px 0;"
        "box-shadow:0 1px 2px rgba(0,0,0,0.08);}"
        ".dash h2{margin:0 0 10px 0;font-size:20px;}"
        ".dash-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:12px;}"
        ".chart{border:1px solid #d8dee4;border-radius:8px;padding:4px;background:#fff;min-height:280px;}"
        ".chart-wide{grid-column:1/-1;min-height:320px;}"
        "img{max-width:100%;height:auto;border:1px solid #d8dee4;border-radius:8px;background:#fff;}"
        ".empty{background:#ffffff;border-radius:10px;padding:18px;font-size:15px;}"
        "code{background:#edf2f7;padding:1px 4px;border-radius:4px;}"
        ".review-controls{background:#ffffff;border-radius:10px;padding:14px;margin:14px 0;"
        "box-shadow:0 1px 2px rgba(0,0,0,0.08);}"
        ".alt-report{background:#ffffff;border-radius:10px;padding:14px;margin:14px 0;"
        "box-shadow:0 1px 2px rgba(0,0,0,0.08);}"
        ".alt-report h2{margin:0 0 10px 0;font-size:20px;}"
        ".review-controls h2{margin:0 0 10px 0;font-size:20px;}"
        ".review-row{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-bottom:8px;}"
        ".filter-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:12px;margin:12px 0;}"
        ".range-filter{border:1px solid #d8dee4;border-radius:8px;padding:10px;background:#f8fafc;}"
        ".range-filter.is-disabled{opacity:0.65;}"
        ".range-header{display:flex;justify-content:space-between;gap:8px;align-items:center;margin-bottom:6px;font-size:13px;}"
        ".range-label{font-weight:700;color:#14212c;}"
        ".range-readout{font-variant-numeric:tabular-nums;color:#405160;}"
        ".range-sliders{display:grid;gap:6px;margin:8px 0;}"
        ".range-sliders label{display:grid;gap:4px;font-size:12px;color:#405160;}"
        ".range-sliders input[type=range]{width:100%;}"
        ".range-values{display:grid;grid-template-columns:1fr 1fr auto;gap:8px;align-items:end;}"
        ".range-values label{display:grid;gap:4px;font-size:12px;color:#405160;}"
        ".range-values input{border:1px solid #c8d1da;border-radius:6px;padding:6px 8px;background:#fff;color:#14212c;font-size:13px;}"
        ".range-reset{border:1px solid #c2ccd6;border-radius:6px;padding:6px 10px;background:#fff;color:#243746;cursor:pointer;font-size:12px;}"
        ".range-reset:hover{border-color:#8fa4b6;}"
        ".review-select{border:1px solid #c8d1da;border-radius:6px;padding:6px 10px;background:#fff;color:#14212c;font-size:14px;}"
        ".review-summary{font-size:13px;color:#405160;margin:4px 0 10px 0;}"
        ".review-buttons{display:flex;flex-wrap:wrap;gap:8px;margin:8px 0 10px 0;}"
        ".review-btn{border:1px solid #c2ccd6;border-radius:999px;padding:4px 10px;background:#fff;color:#243746;cursor:pointer;font-size:12px;}"
        ".review-btn:hover{border-color:#8fa4b6;}"
        ".review-btn.is-active{color:#fff;border-color:transparent;}"
        ".review-btn.review-real.is-active{background:#1a7f37;}"
        ".review-btn.review-not-real.is-active{background:#b42318;}"
        ".review-btn.review-undecided.is-active{background:#475467;}"
        ".review-state-real{border-left:6px solid #1a7f37;}"
        ".review-state-not-real{border-left:6px solid #b42318;}"
        ".review-state-undecided{border-left:6px solid #98a2b3;}"
        ".review-badge{text-transform:lowercase;font-weight:700;}"
        ".igv-review-layout{display:flex;flex-direction:column;gap:10px;margin-top:10px;}"
        ".igv-side-actions{background:#f8fafc;border:1px solid #d8dee4;border-radius:8px;padding:8px;}"
        ".igv-side-actions .review-buttons{margin-top:6px;margin-bottom:0;}"
        ".igv-grid{display:flex;flex-direction:column;gap:10px;}"
        ".igv-card{background:#f8fafc;border:1px solid #d8dee4;border-radius:8px;padding:8px;}"
        ".igv-card img{display:block;width:100%;}"
    )
    page.append("</style>")
    page.append("</head>")
    page.append("<body>")
    page.append("<div class=\"wrap\">")
    page.append("<h1>SniffCell High-Confidence SV Report</h1>")
    page.append(
        "<div class=\"meta\">"
        f"Generated: {html.escape(generated_at)}"
        f" | anno_output: <code>{html.escape(str(anno_output))}</code>"
        f" | sv_assignment: <code>{html.escape(str(sv_assignment_path))}</code>"
        "</div>"
    )
    page.append("<div class=\"stats\">")
    page.append(f"<div class=\"card\"><div class=\"label\">Total SVs</div><div class=\"value\">{total_sv}</div></div>")
    page.append(f"<div class=\"card\"><div class=\"label\">Selected SVs</div><div class=\"value\">{selected_count}</div></div>")
    page.append(f"<div class=\"card\"><div class=\"label\">Rendered Viz</div><div class=\"value\">{rendered_count}</div></div>")
    page.append(f"<div class=\"card\"><div class=\"label\">Viz Failures</div><div class=\"value\">{failed_count}</div></div>")
    page.append("</div>")
    page.append(
        "<div class=\"filter\">"
        f"Filters: min_overlap_pct={filters['min_overlap_pct']}, "
        f"min_majority_pct={filters['min_majority_pct']}, "
        f"include_unassigned={filters['include_unassigned']}, "
        f"allow_hard_conflict={filters['allow_hard_conflict']}, "
        f"max_sv={filters['max_sv']}, "
        f"with_figures={viz['with_figures']}, "
        f"figure_threads={viz['figure_threads']}, "
        f"figure_profile={viz['figure_profile']}, "
        f"figure_dpi={viz['figure_dpi']}, "
        f"skip_methylation_overlay={viz['skip_methylation_overlay']}, "
        f"exact_window={viz['exact_window']}, "
        f"with_igvreport={igvreport['enabled']}"
        "</div>"
    )
    igvreport_status = str(igvreport.get("status", "")).strip()
    if bool(igvreport.get("enabled")) or igvreport_status in {"existing", "rendered", "reused", "failed", "failed_no_output"}:
        page.append("<section class=\"alt-report\">")
        page.append("<h2>Alternate IGV Report</h2>")
        page.append(f"<div class=\"kv\"><b>igvreport status:</b> {html.escape(igvreport_status or 'not_rendered')}</div>")
        igvreport_error = str(igvreport.get("error", "")).strip()
        if igvreport_error:
            page.append(f"<div class=\"kv err\">igvreport error: {html.escape(igvreport_error)}</div>")
        igvreport_html_rel = str(igvreport.get("html_rel", "")).strip()
        if igvreport_html_rel:
            page.append(
                f"<div class=\"kv\"><b>igvreport HTML:</b> "
                f"<a href=\"{html.escape(igvreport_html_rel)}\" target=\"_blank\" rel=\"noopener\">"
                f"{html.escape(igvreport_html_rel)}</a></div>"
            )
        igvreport_manifest_rel = str(igvreport.get("manifest_rel", "")).strip()
        if igvreport_manifest_rel:
            page.append(
                f"<div class=\"kv\"><b>igvreport manifest:</b> "
                f"<a href=\"{html.escape(igvreport_manifest_rel)}\" target=\"_blank\" rel=\"noopener\">"
                f"{html.escape(igvreport_manifest_rel)}</a></div>"
            )
        igvreport_command = str(igvreport.get("command", "")).strip()
        if igvreport_command:
            escaped_igvreport_command = html.escape(igvreport_command, quote=True)
            page.append(f"<div class=\"cmd\"><code>{html.escape(igvreport_command)}</code></div>")
            page.append(
                f"<button class=\"copy\" type=\"button\" data-cmd=\"{escaped_igvreport_command}\" "
                "onclick=\"copyVizCommand(this)\">Copy igvreport command</button>"
            )
        page.append("</section>")
    page.append("<section class=\"review-controls\">")
    page.append("<h2>SV Review Controls</h2>")
    page.append("<div class=\"review-row\">")
    page.append("<label for=\"review-filter\"><b>Show SVs:</b></label>")
    page.append(
        "<select id=\"review-filter\" class=\"review-select\" onchange=\"applyReviewFilter()\">"
        "<option value=\"all\">All</option>"
        "<option value=\"real\">Real</option>"
        "<option value=\"not_real\">Not real</option>"
        "<option value=\"undecided\">Undecided</option>"
        "</select>"
    )
    page.append("<label for=\"celltype-filter\"><b>Assigned cell type:</b></label>")
    page.append(
        "<select id=\"celltype-filter\" class=\"review-select\" onchange=\"applyReviewFilter()\">"
        "<option value=\"all\">All assigned cell types</option>"
        "</select>"
    )
    page.append("<label for=\"svtype-filter\"><b>SV type:</b></label>")
    page.append(
        "<select id=\"svtype-filter\" class=\"review-select\" onchange=\"applyReviewFilter()\">"
        "<option value=\"all\">All SV types</option>"
        "</select>"
    )
    page.append("<label for=\"hard-conflict-filter\"><b>Hard conflict:</b></label>")
    page.append(
        "<select id=\"hard-conflict-filter\" class=\"review-select\" onchange=\"applyReviewFilter()\">"
        "<option value=\"all\">All hard-conflict states</option>"
        "<option value=\"false\">No hard conflict</option>"
        "<option value=\"true\">Hard conflict</option>"
        "<option value=\"unknown\">Unknown</option>"
        "</select>"
    )
    page.append("</div>")
    page.append("<div id=\"numeric-filter-grid\" class=\"filter-grid\"></div>")
    page.append("<div id=\"review-summary\" class=\"review-summary\"></div>")
    page.append("<div id=\"review-persist-status\" class=\"review-summary\"></div>")
    page.append("</section>")
    page.append("<section class=\"dash\">")
    page.append("<h2>Interactive Summaries</h2>")
    page.append("<div class=\"dash-grid\">")
    page.append("<div id=\"chart-genome-location\" class=\"chart chart-wide\"></div>")
    page.append("<div id=\"chart-chrom-counts\" class=\"chart\"></div>")
    page.append("<div id=\"chart-svlen\" class=\"chart\"></div>")
    page.append("<div id=\"chart-support\" class=\"chart\"></div>")
    page.append("<div id=\"chart-overlap-majority\" class=\"chart\"></div>")
    page.append("<div id=\"chart-celltype\" class=\"chart\"></div>")
    page.append("</div>")
    page.append("</section>")

    if not rows:
        page.append("<div class=\"empty\">No SVs passed the report filters.</div>")
    else:
        for item in rows:
            sv_id_text = str(item["id"])
            sv_id = html.escape(sv_id_text)
            sv_id_attr = html.escape(sv_id_text, quote=True)
            linked_text = _fmt_text(item.get("linked_celltypes", ""))
            primary_text = _fmt_text(item.get("primary_celltype", ""))
            assigned_code_text = _fmt_text(item.get("assigned_code", ""))
            linked = html.escape(linked_text)
            primary = html.escape(primary_text)
            assigned_code = html.escape(assigned_code_text)
            majority = _fmt_float(item.get("majority_pct"))
            overlap = _fmt_float(item.get("overlap_pct"))
            vaf = _fmt_float(item.get("vaf"))
            n_supporting = "NA" if pd.isna(item.get("n_supporting")) else int(item["n_supporting"])
            n_overlapped = "NA" if pd.isna(item.get("n_overlapped")) else int(item["n_overlapped"])
            status = html.escape(str(item.get("viz_status", "")))
            n_supporting_text = str(n_supporting)
            n_overlapped_text = str(n_overlapped)
            status_text = _fmt_text(item.get("viz_status", ""))
            sv_type_text = _fmt_text(item.get("sv_type", ""))
            sv_type = html.escape(sv_type_text)
            hard_conflict_token = "unknown"
            hard_conflict_display = "NA"
            hard_conflict_raw = item.get("has_hard_conflict", pd.NA)
            try:
                if pd.isna(hard_conflict_raw):
                    hard_conflict_token = "unknown"
                    hard_conflict_display = "NA"
                else:
                    hard_conflict_token = "true" if bool(hard_conflict_raw) else "false"
                    hard_conflict_display = "True" if bool(hard_conflict_raw) else "False"
            except TypeError:
                hard_conflict_token = "unknown"
                hard_conflict_display = "NA"
            sv_len_text = "NA"
            sv_len_abs_text = "NA"
            try:
                sv_len_raw_for_display = item.get("sv_len", pd.NA)
                if not pd.isna(sv_len_raw_for_display):
                    sv_len_int = int(float(sv_len_raw_for_display))
                    sv_len_text = str(sv_len_int)
                    sv_len_abs_text = str(abs(sv_len_int))
            except Exception:
                sv_len_text = "NA"
                sv_len_abs_text = "NA"
            sv_len_display = (
                f"{sv_len_text} bp (abs {sv_len_abs_text} bp)"
                if sv_len_text != "NA"
                else "NA"
            )
            sv_igv_text = "NA"
            try:
                sv_chr_raw = item.get("sv_chr", pd.NA)
                sv_pos_raw = item.get("sv_pos", pd.NA)
                sv_len_raw = item.get("sv_len", pd.NA)
                if (not pd.isna(sv_chr_raw)) and (not pd.isna(sv_pos_raw)):
                    sv_chr_text = str(sv_chr_raw).strip()
                    sv_start_1 = int(float(sv_pos_raw))
                    if pd.isna(sv_len_raw):
                        sv_end_1 = sv_start_1
                    else:
                        sv_end_1 = sv_start_1 + max(1, abs(int(float(sv_len_raw)))) - 1
                    sv_igv_text = f"{sv_chr_text}:{sv_start_1}-{sv_end_1}"
            except Exception:
                sv_igv_text = "NA"

            review_status = _normalize_review_status(item.get("review_status", "undecided"))
            review_class = _review_state_class(review_status)
            review_label = _review_state_label(review_status)
            real_active = " is-active" if review_status == "real" else ""
            not_real_active = " is-active" if review_status == "not_real" else ""
            undecided_active = " is-active" if review_status == "undecided" else ""
            page.append(
                f"<section class=\"sv {review_class}\" "
                f"data-sv-id=\"{sv_id_attr}\" "
                f"data-review-status=\"{html.escape(review_status, quote=True)}\" "
                f"data-primary-celltype=\"{html.escape(primary_text, quote=True)}\" "
                f"data-linked-celltypes=\"{html.escape(linked_text, quote=True)}\" "
                f"data-assigned-code=\"{html.escape(assigned_code_text, quote=True)}\" "
                f"data-majority-pct=\"{html.escape(majority, quote=True)}\" "
                f"data-overlap-pct=\"{html.escape(overlap, quote=True)}\" "
                f"data-vaf=\"{html.escape(vaf, quote=True)}\" "
                f"data-n-supporting=\"{html.escape(n_supporting_text, quote=True)}\" "
                f"data-n-overlapped=\"{html.escape(n_overlapped_text, quote=True)}\" "
                f"data-sv-type=\"{html.escape(sv_type_text, quote=True)}\" "
                f"data-has-hard-conflict=\"{html.escape(hard_conflict_token, quote=True)}\" "
                f"data-sv-len=\"{html.escape(sv_len_text, quote=True)}\" "
                f"data-viz-status=\"{html.escape(status_text, quote=True)}\">"
            )
            page.append(f"<h2>{sv_id}</h2>")
            page.append("<div class=\"review-buttons\">")
            page.append(
                f"<button class=\"review-btn review-real{real_active}\" type=\"button\" data-review-value=\"real\" "
                "onclick=\"setSvReview(this)\">Real</button>"
            )
            page.append(
                f"<button class=\"review-btn review-not-real{not_real_active}\" type=\"button\" data-review-value=\"not_real\" "
                "onclick=\"setSvReview(this)\">Not real</button>"
            )
            page.append(
                f"<button class=\"review-btn review-undecided{undecided_active}\" type=\"button\" data-review-value=\"undecided\" "
                "onclick=\"setSvReview(this)\">Undecided</button>"
            )
            page.append("</div>")
            page.append(f"<div class=\"kv\"><b>Review:</b> <span class=\"review-badge\">{html.escape(review_label)}</span></div>")
            page.append(f"<div class=\"kv\"><b>Primary cell type:</b> {primary}</div>")
            page.append(f"<div class=\"kv\"><b>Linked cell types:</b> {linked}</div>")
            page.append(f"<div class=\"kv\"><b>SV type:</b> {sv_type}</div>")
            page.append(f"<div class=\"kv\"><b>Assigned code:</b> {assigned_code}</div>")
            page.append(f"<div class=\"kv\"><b>SV length:</b> {html.escape(sv_len_display)}</div>")
            page.append(f"<div class=\"kv\"><b>IGV SV locus:</b> <code>{html.escape(sv_igv_text)}</code></div>")
            page.append(
                f"<div class=\"kv\"><b>majority_pct:</b> {majority} | <b>overlap_pct:</b> {overlap} | "
                f"<b>vaf:</b> {vaf}</div>"
            )
            page.append(
                f"<div class=\"kv\"><b>n_supporting:</b> {html.escape(n_supporting_text)} | "
                f"<b>n_overlapped:</b> {html.escape(n_overlapped_text)}</div>"
            )
            page.append(f"<div class=\"kv\"><b>has_hard_conflict:</b> {html.escape(hard_conflict_display)}</div>")
            page.append(f"<div class=\"kv\"><b>viz status:</b> {status}</div>")
            igvviz_status = html.escape(str(item.get("igvviz_status", "")))
            page.append(f"<div class=\"kv\"><b>igvviz status:</b> {igvviz_status}</div>")

            err = str(item.get("viz_error", "")).strip()
            if err:
                page.append(f"<div class=\"kv err\">viz error: {html.escape(err)}</div>")
            igv_err = str(item.get("igvviz_error", "")).strip()
            if igv_err:
                page.append(f"<div class=\"kv err\">igvviz error: {html.escape(igv_err)}</div>")

            fig_rel = str(item.get("viz_figure_rel", "")).strip()
            if fig_rel:
                page.append(f"<div style=\"margin-top:10px\"><img src=\"{html.escape(fig_rel)}\" alt=\"SV plot for {sv_id}\" loading=\"lazy\"></div>")
            igv_manifest_rel = str(item.get("igvviz_manifest_rel", "")).strip()
            if igv_manifest_rel:
                page.append(
                    f"<div class=\"kv\" style=\"margin-top:8px\"><b>igvviz manifest:</b> "
                    f"<a href=\"{html.escape(igv_manifest_rel)}\" target=\"_blank\" rel=\"noopener\">{html.escape(igv_manifest_rel)}</a></div>"
                )
            igv_snapshots = item.get("_igvviz_snapshots", [])
            if isinstance(igv_snapshots, list):
                shown = [x for x in igv_snapshots if isinstance(x, dict) and str(x.get("snapshot_rel", "")).strip()]
                if shown:
                    page.append("<div class=\"igv-review-layout\">")
                    page.append("<div class=\"igv-side-actions\">")
                    page.append("<div class=\"kv\"><b>Review controls (IGV):</b></div>")
                    page.append("<div class=\"review-buttons\">")
                    page.append(
                        f"<button class=\"review-btn review-real{real_active}\" type=\"button\" data-review-value=\"real\" "
                        "onclick=\"setSvReview(this)\">Real</button>"
                    )
                    page.append(
                        f"<button class=\"review-btn review-not-real{not_real_active}\" type=\"button\" data-review-value=\"not_real\" "
                        "onclick=\"setSvReview(this)\">Not real</button>"
                    )
                    page.append(
                        f"<button class=\"review-btn review-undecided{undecided_active}\" type=\"button\" data-review-value=\"undecided\" "
                        "onclick=\"setSvReview(this)\">Undecided</button>"
                    )
                    page.append("</div>")
                    page.append("</div>")
                    page.append("<div class=\"igv-grid\">")
                    for snap in shown:
                        snap_rel = html.escape(str(snap.get("snapshot_rel", "")).strip())
                        bam_label = html.escape(str(snap.get("bam_label", "IGV")).strip() or "IGV")
                        page.append("<div class=\"igv-card\">")
                        page.append(f"<div class=\"kv\"><b>{bam_label}</b></div>")
                        page.append(f"<img src=\"{snap_rel}\" alt=\"IGV snapshot {bam_label} for {sv_id}\" loading=\"lazy\">")
                        page.append("</div>")
                    page.append("</div>")
                    page.append("</div>")

            viz_command = str(item.get("viz_command", "")).strip()
            if viz_command:
                escaped_command = html.escape(viz_command, quote=True)
                page.append(f"<div class=\"cmd\"><code>{html.escape(viz_command)}</code></div>")
                page.append(
                    f"<button class=\"copy\" type=\"button\" data-cmd=\"{escaped_command}\" "
                    "onclick=\"copyVizCommand(this)\">Copy viz command</button>"
                )
            igvviz_command = str(item.get("igvviz_command", "")).strip()
            if igvviz_command:
                escaped_igv_command = html.escape(igvviz_command, quote=True)
                page.append(f"<div class=\"cmd\"><code>{html.escape(igvviz_command)}</code></div>")
                page.append(
                    f"<button class=\"copy\" type=\"button\" data-cmd=\"{escaped_igv_command}\" "
                    "onclick=\"copyVizCommand(this)\">Copy igvviz command</button>"
                )
            page.append("</section>")

    page.append("</div>")
    page.append(
        "<script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>"
    )
    script_text = (
        """
const dashboardData=__DASHBOARD_DATA__;
const dashboardFilters=__DASHBOARD_FILTERS__;
const REVIEW_STORAGE_KEY=__REVIEW_STORAGE_KEY__;
const REVIEW_STATUS_LABELS={real:'real',not_real:'not real',undecided:'undecided'};
const NUMERIC_FILTER_CONFIGS=[
  {key:'n_supporting',label:'n_supporting',digits:0,step:1,source:'n_supporting'},
  {key:'overlap_pct',label:'overlap_pct',digits:3,step:0.001,source:'overlap_pct',fixed:[0,1]},
  {key:'majority_pct',label:'majority_pct',digits:3,step:0.001,source:'majority_pct',fixed:[0,1]},
  {key:'sv_len_abs',label:'|sv_len|',digits:0,step:1,source:'sv_len'},
  {key:'vaf',label:'vaf',digits:3,step:0.001,source:'vaf',fixed:[0,1]}
];
let reviewState={};
let numericFilterMeta={};
let numericFilterState={};

function chrSortKey(chr){
  const t=String(chr||'').replace(/^chr/i,'');
  if(/^\\d+$/.test(t)) return [0, parseInt(t,10), t];
  if(t==='X') return [1,23,t];
  if(t==='Y') return [1,24,t];
  if(t==='M'||t==='MT') return [1,25,t];
  return [2, Number.MAX_SAFE_INTEGER, t];
}

function sortChromLabels(labels){
  return labels.slice().sort((a,b)=>{
    const ka=chrSortKey(a);
    const kb=chrSortKey(b);
    if(ka[0]!==kb[0]) return ka[0]-kb[0];
    if(ka[1]!==kb[1]) return ka[1]-kb[1];
    return String(ka[2]).localeCompare(String(kb[2]));
  });
}

function getSvCards(){
  return Array.from(document.querySelectorAll('.sv[data-sv-id]'));
}

function normalizeReviewStatus(value){
  const v=String(value||'').toLowerCase();
  if(v==='real') return 'real';
  if(v==='not_real'||v==='not-real'||v==='not real') return 'not_real';
  return 'undecided';
}

function reviewLabel(status){
  return REVIEW_STATUS_LABELS[status]||'undecided';
}

function escapeHtml(text){
  const span=document.createElement('span');
  span.textContent=String(text==null?'':text);
  return span.innerHTML;
}

function normalizeCelltypeToken(value){
  return String(value==null?'':value).trim().toLowerCase();
}

function splitCelltypeList(text){
  return String(text==null?'':text)
    .split(/[;,|]+/)
    .map(part=>String(part||'').trim())
    .filter(part=>part!==''&&part.toLowerCase()!=='na');
}

function asFiniteNumber(value){
  const num=Number(value);
  return Number.isFinite(num)?num:null;
}

function roundNumericValue(value, digits){
  const num=asFiniteNumber(value);
  if(num==null) return null;
  const factor=Math.pow(10, Number(digits||0));
  return Math.round(num*factor)/factor;
}

function formatNumericValue(value, digits){
  const num=asFiniteNumber(value);
  if(num==null) return 'NA';
  const rounded=roundNumericValue(num, digits);
  return Number(digits||0)===0 ? String(Math.round(rounded)) : rounded.toFixed(Number(digits||0));
}

function normalizeBooleanToken(value){
  const text=String(value==null?'':value).trim().toLowerCase();
  if(text==='true'||text==='1'||text==='yes') return 'true';
  if(text==='false'||text==='0'||text==='no') return 'false';
  return 'unknown';
}

function numericValueByKey(source, key){
  if(key==='sv_len_abs'){
    const value=asFiniteNumber(source&&source.sv_len);
    return value==null?null:Math.abs(value);
  }
  return asFiniteNumber(source&&source[key]);
}

function datasetNumericValue(card, attrName, absolute){
  const raw=card.getAttribute(attrName);
  const value=asFiniteNumber(raw);
  if(value==null) return null;
  return absolute?Math.abs(value):value;
}

function numericValueForCard(card, key){
  if(key==='n_supporting') return datasetNumericValue(card, 'data-n-supporting', false);
  if(key==='overlap_pct') return datasetNumericValue(card, 'data-overlap-pct', false);
  if(key==='majority_pct') return datasetNumericValue(card, 'data-majority-pct', false);
  if(key==='sv_len_abs') return datasetNumericValue(card, 'data-sv-len', true);
  if(key==='vaf') return datasetNumericValue(card, 'data-vaf', false);
  return null;
}

function buildNumericFilterMeta(){
  const out={};
  NUMERIC_FILTER_CONFIGS.forEach(cfg=>{
    const values=(Array.isArray(dashboardData)?dashboardData:[])
      .map(row=>numericValueByKey(row, cfg.key))
      .filter(value=>value!=null);
    if(values.length===0){
      out[cfg.key]={...cfg,enabled:false,min:null,max:null,scaleMin:null,scaleMax:null};
      return;
    }
    let scaleMin=Math.min(...values);
    let scaleMax=Math.max(...values);
    if(Array.isArray(cfg.fixed)&&cfg.fixed.length===2){
      scaleMin=Number(cfg.fixed[0]);
      scaleMax=Number(cfg.fixed[1]);
    }
    scaleMin=roundNumericValue(scaleMin, cfg.digits);
    scaleMax=roundNumericValue(scaleMax, cfg.digits);
    if(scaleMax<scaleMin) scaleMax=scaleMin;
    out[cfg.key]={
      ...cfg,
      enabled:true,
      min:scaleMin,
      max:scaleMax,
      scaleMin:scaleMin,
      scaleMax:scaleMax
    };
  });
  return out;
}

function selectedReviewFilter(){
  const sel=document.getElementById('review-filter');
  return sel?String(sel.value||'all'):'all';
}

function selectedCelltypeFilter(){
  const sel=document.getElementById('celltype-filter');
  return sel?String(sel.value||'all'):'all';
}

function selectedSvtypeFilter(){
  const sel=document.getElementById('svtype-filter');
  return sel?String(sel.value||'all'):'all';
}

function selectedHardConflictFilter(){
  const sel=document.getElementById('hard-conflict-filter');
  return sel?String(sel.value||'all'):'all';
}

function cardCelltypeTokens(card){
  const out=new Set();
  splitCelltypeList(card.getAttribute('data-primary-celltype')||'')
    .forEach(v=>out.add(normalizeCelltypeToken(v)));
  splitCelltypeList(card.getAttribute('data-linked-celltypes')||'')
    .forEach(v=>out.add(normalizeCelltypeToken(v)));
  return out;
}

function rowCelltypeTokens(row){
  const out=new Set();
  splitCelltypeList(row&&row.primary_celltype)
    .forEach(v=>out.add(normalizeCelltypeToken(v)));
  splitCelltypeList(row&&row.linked_celltypes)
    .forEach(v=>out.add(normalizeCelltypeToken(v)));
  return out;
}

function cardMatchesCelltype(card, wantedToken){
  const wanted=normalizeCelltypeToken(wantedToken);
  if(wanted===''||wanted==='all') return true;
  return cardCelltypeTokens(card).has(wanted);
}

function rowMatchesCelltype(row, wantedToken){
  const wanted=normalizeCelltypeToken(wantedToken);
  if(wanted===''||wanted==='all') return true;
  return rowCelltypeTokens(row).has(wanted);
}

function cardMatchesSvType(card, wantedValue){
  const wanted=String(wantedValue||'all').trim().toUpperCase();
  if(wanted===''||wanted==='ALL') return true;
  const actual=String(card.getAttribute('data-sv-type')||'').trim().toUpperCase();
  return actual===wanted;
}

function rowMatchesSvType(row, wantedValue){
  const wanted=String(wantedValue||'all').trim().toUpperCase();
  if(wanted===''||wanted==='ALL') return true;
  const actual=String((row&&row.sv_type)||'').trim().toUpperCase();
  return actual===wanted;
}

function cardMatchesHardConflict(card, wantedValue){
  const wanted=normalizeBooleanToken(wantedValue);
  if(wanted==='unknown' && String(wantedValue||'all')==='all') return true;
  if(String(wantedValue||'all')==='all') return true;
  const actual=normalizeBooleanToken(card.getAttribute('data-has-hard-conflict'));
  return actual===wanted;
}

function rowMatchesHardConflict(row, wantedValue){
  if(String(wantedValue||'all')==='all') return true;
  const actual=normalizeBooleanToken(row&&row.has_hard_conflict);
  const wanted=normalizeBooleanToken(wantedValue);
  return actual===wanted;
}

function numericFiltersMatchValue(value, bounds, meta){
  if(!bounds||!meta) return true;
  if(value==null){
    return bounds.min===meta.scaleMin && bounds.max===meta.scaleMax;
  }
  if(bounds.min!=null && value<bounds.min) return false;
  if(bounds.max!=null && value>bounds.max) return false;
  return true;
}

function cardMatchesNumericFilters(card){
  return NUMERIC_FILTER_CONFIGS.every(cfg=>{
    const meta=numericFilterMeta[cfg.key];
    if(!meta||!meta.enabled) return true;
    const bounds=numericFilterState[cfg.key];
    return numericFiltersMatchValue(numericValueForCard(card, cfg.key), bounds, meta);
  });
}

function rowMatchesNumericFilters(row){
  return NUMERIC_FILTER_CONFIGS.every(cfg=>{
    const meta=numericFilterMeta[cfg.key];
    if(!meta||!meta.enabled) return true;
    const bounds=numericFilterState[cfg.key];
    return numericFiltersMatchValue(numericValueByKey(row, cfg.key), bounds, meta);
  });
}

function filteredDashboardData(filterState){
  const wantedReview=String((filterState&&filterState.review)||'all');
  const wantedCelltype=String((filterState&&filterState.celltype)||'all');
  const wantedSvType=String((filterState&&filterState.svtype)||'all');
  const wantedHardConflict=String((filterState&&filterState.hardConflict)||'all');
  const base=Array.isArray(dashboardData)?dashboardData:[];
  return base.filter(r=>{
    const id=String((r&&r.id!=null)?r.id:'');
    const status=normalizeReviewStatus(reviewState[id]);
    if(wantedReview!=='all'&&status!==wantedReview) return false;
    if(!rowMatchesCelltype(r, wantedCelltype)) return false;
    if(!rowMatchesSvType(r, wantedSvType)) return false;
    if(!rowMatchesHardConflict(r, wantedHardConflict)) return false;
    if(!rowMatchesNumericFilters(r)) return false;
    return true;
  });
}

function numericFilterSummaryParts(){
  return NUMERIC_FILTER_CONFIGS
    .map(cfg=>{
      const meta=numericFilterMeta[cfg.key];
      const bounds=numericFilterState[cfg.key];
      if(!meta||!meta.enabled||!bounds) return '';
      const atDefault=bounds.min===meta.scaleMin && bounds.max===meta.scaleMax;
      if(atDefault) return '';
      return `${cfg.label}: ${formatNumericValue(bounds.min, cfg.digits)}-${formatNumericValue(bounds.max, cfg.digits)}`;
    })
    .filter(Boolean);
}

function summaryScopeLabel(filterState){
  const parts=[];
  const wantedReview=String((filterState&&filterState.review)||'all');
  const wantedCelltype=String((filterState&&filterState.celltype)||'all');
  const wantedSvType=String((filterState&&filterState.svtype)||'all');
  const wantedHardConflict=String((filterState&&filterState.hardConflict)||'all');
  if(wantedReview==='real') parts.push('real');
  else if(wantedReview==='not_real') parts.push('not real');
  else if(wantedReview==='undecided') parts.push('undecided');
  if(wantedCelltype!=='all'){
    let cellLabel=wantedCelltype;
    const sel=document.getElementById('celltype-filter');
    if(sel&&sel.selectedOptions&&sel.selectedOptions.length>0){
      cellLabel=String(sel.selectedOptions[0].textContent||wantedCelltype);
    }
    parts.push(`cell type: ${cellLabel}`);
  }
  if(wantedSvType!=='all') parts.push(`SV type: ${wantedSvType}`);
  if(wantedHardConflict!=='all'){
    parts.push(
      wantedHardConflict==='true'
        ?'hard conflict'
        :(wantedHardConflict==='false'?'no hard conflict':'hard conflict: unknown')
    );
  }
  parts.push(...numericFilterSummaryParts());
  return parts.length===0?'all':parts.join(' | ');
}

function currentDashboardFilterState(){
  return {
    review:selectedReviewFilter(),
    celltype:selectedCelltypeFilter(),
    svtype:selectedSvtypeFilter(),
    hardConflict:selectedHardConflictFilter()
  };
}

function renderSummaries(filterState){
  if(typeof Plotly==='undefined'){
    const msg='Plotly failed to load; interactive plots unavailable.';
    ['chart-genome-location','chart-chrom-counts','chart-svlen','chart-support','chart-overlap-majority','chart-celltype']
      .forEach(id=>{const el=document.getElementById(id); if(el){el.textContent=msg;}});
    return;
  }
  const data=filteredDashboardData(filterState);
  const scope=summaryScopeLabel(filterState);
  const titleSuffix=(scope==='all')?'':` [${scope}]`;
  if(data.length===0){
    const msg=(scope==='all')
      ?'No selected SVs for summary plots.'
      :`No SVs in '${scope}' for summary plots.`;
    ['chart-genome-location','chart-chrom-counts','chart-svlen','chart-support','chart-overlap-majority','chart-celltype']
      .forEach(id=>{const el=document.getElementById(id); if(el){el.textContent=msg;}});
    return;
  }

  const posRows=data.filter(r=>r.sv_chr!=null&&r.sv_pos!=null&&Number.isFinite(Number(r.sv_pos)));
  if(posRows.length>0){
    const yMb=posRows.map(r=>Number(r.sv_pos)/1e6);
    const txt=posRows.map(r=>`${r.id||'NA'}<br>${r.sv_chr}:${r.sv_pos}`);
    Plotly.newPlot(
      'chart-genome-location',
      [{type:'scatter',mode:'markers',x:posRows.map(r=>String(r.sv_chr)),y:yMb,text:txt,hovertemplate:'%{text}<br>Position(Mb): %{y:.3f}<extra></extra>',marker:{size:8,color:'#1f77b4',opacity:0.75}}],
      {title:'Genome-wide SV Locations (selected)'+titleSuffix,xaxis:{title:'Chromosome',categoryorder:'array',categoryarray:sortChromLabels(posRows.map(r=>String(r.sv_chr)))},yaxis:{title:'SV position (Mb)'}},
      {responsive:true,displaylogo:false}
    );
  } else {
    document.getElementById('chart-genome-location').textContent='No sv_chr/sv_pos data.';
  }

  const chrCounts={};
  data.forEach(r=>{
    const c=(r.sv_chr==null||String(r.sv_chr).trim()==='')?'NA':String(r.sv_chr);
    chrCounts[c]=(chrCounts[c]||0)+1;
  });
  const chrLabels=sortChromLabels(Object.keys(chrCounts));
  Plotly.newPlot(
    'chart-chrom-counts',
    [{type:'bar',x:chrLabels,y:chrLabels.map(c=>chrCounts[c]),marker:{color:'#2ca02c'}}],
    {title:'SV Count by Chromosome'+titleSuffix,xaxis:{title:'Chromosome'},yaxis:{title:'SV count'}},
    {responsive:true,displaylogo:false}
  );

  const lenVals=data.map(r=>Math.abs(Number(r.sv_len))).filter(v=>Number.isFinite(v)&&v>0).map(v=>Math.log10(v));
  if(lenVals.length>0){
    Plotly.newPlot(
      'chart-svlen',
      [{type:'histogram',x:lenVals,marker:{color:'#9467bd'}}],
      {title:'SV Length Distribution'+titleSuffix,xaxis:{title:'log10(|sv_len| bp)'},yaxis:{title:'Count'}},
      {responsive:true,displaylogo:false}
    );
  } else {
    document.getElementById('chart-svlen').textContent='No sv_len data.';
  }

  const nSup=data.map(r=>Number(r.n_supporting)).filter(v=>Number.isFinite(v));
  const nOvl=data.map(r=>Number(r.n_overlapped)).filter(v=>Number.isFinite(v));
  if(nSup.length+nOvl.length>0){
    Plotly.newPlot(
      'chart-support',
      [{type:'histogram',x:nSup,name:'n_supporting',opacity:0.65,marker:{color:'#ff7f0e'}},{type:'histogram',x:nOvl,name:'n_overlapped',opacity:0.65,marker:{color:'#17becf'}}],
      {title:'Read Support Distribution'+titleSuffix,xaxis:{title:'Read count'},yaxis:{title:'SV count'},barmode:'overlay'},
      {responsive:true,displaylogo:false}
    );
  } else {
    document.getElementById('chart-support').textContent='No read support data.';
  }

  const omRows=data.filter(r=>Number.isFinite(Number(r.overlap_pct))&&Number.isFinite(Number(r.majority_pct)));
  if(omRows.length>0){
    Plotly.newPlot(
      'chart-overlap-majority',
      [{type:'scatter',mode:'markers',x:omRows.map(r=>Number(r.overlap_pct)),y:omRows.map(r=>Number(r.majority_pct)),text:omRows.map(r=>String(r.id||'')),hovertemplate:'%{text}<br>overlap=%{x:.3f}<br>majority=%{y:.3f}<extra></extra>',marker:{size:9,color:omRows.map(r=>Number.isFinite(Number(r.n_supporting))?Number(r.n_supporting):0),colorscale:'Viridis',showscale:true,colorbar:{title:'n_supporting'}}}],
      {title:'Agreement vs Overlap'+titleSuffix,xaxis:{title:'overlap_pct',range:[0,1]},yaxis:{title:'majority_pct',range:[0,1]},shapes:[{type:'line',x0:Number(dashboardFilters.min_overlap_pct||0),x1:Number(dashboardFilters.min_overlap_pct||0),y0:0,y1:1,line:{dash:'dot',color:'#444'}},{type:'line',y0:Number(dashboardFilters.min_majority_pct||0),y1:Number(dashboardFilters.min_majority_pct||0),x0:0,x1:1,line:{dash:'dot',color:'#444'}}]},
      {responsive:true,displaylogo:false}
    );
  } else {
    document.getElementById('chart-overlap-majority').textContent='No overlap/majority data.';
  }

  const ctCounts={};
  data.forEach(r=>{
    const ct=(r.primary_celltype==null||String(r.primary_celltype).trim()==='')?'NA':String(r.primary_celltype);
    ctCounts[ct]=(ctCounts[ct]||0)+1;
  });
  const ctPairs=Object.entries(ctCounts).sort((a,b)=>b[1]-a[1]).slice(0,15);
  Plotly.newPlot(
    'chart-celltype',
    [{type:'bar',orientation:'h',x:ctPairs.map(p=>p[1]).reverse(),y:ctPairs.map(p=>p[0]).reverse(),marker:{color:'#8c564b'}}],
    {title:'Top Primary Cell Types'+titleSuffix,xaxis:{title:'SV count'},yaxis:{title:'Cell type'}},
    {responsive:true,displaylogo:false}
  );
}

function syncNumericFilterControl(key){
  const meta=numericFilterMeta[key];
  const bounds=numericFilterState[key];
  if(!meta||!meta.enabled||!bounds) return;
  const root=document.querySelector(`.range-filter[data-filter-key="${key}"]`);
  if(!root) return;
  const readout=root.querySelector('.range-readout');
  const minRange=root.querySelector('input[data-role="min-range"]');
  const maxRange=root.querySelector('input[data-role="max-range"]');
  const minInput=root.querySelector('input[data-role="min-input"]');
  const maxInput=root.querySelector('input[data-role="max-input"]');
  if(readout){
    readout.textContent=`${formatNumericValue(bounds.min, meta.digits)} to ${formatNumericValue(bounds.max, meta.digits)}`;
  }
  if(minRange) minRange.value=String(bounds.min);
  if(maxRange) maxRange.value=String(bounds.max);
  if(minInput) minInput.value=formatNumericValue(bounds.min, meta.digits);
  if(maxInput) maxInput.value=formatNumericValue(bounds.max, meta.digits);
}

function setNumericFilterState(key, bound, rawValue){
  const meta=numericFilterMeta[key];
  if(!meta||!meta.enabled) return;
  const fallback=(bound==='min')?meta.scaleMin:meta.scaleMax;
  let value=asFiniteNumber(rawValue);
  if(value==null) value=fallback;
  value=Math.max(meta.scaleMin, Math.min(meta.scaleMax, roundNumericValue(value, meta.digits)));
  const current=numericFilterState[key]||{min:meta.scaleMin,max:meta.scaleMax};
  let minValue=current.min;
  let maxValue=current.max;
  if(bound==='min'){
    minValue=value;
    if(minValue>maxValue) maxValue=minValue;
  }else{
    maxValue=value;
    if(maxValue<minValue) minValue=maxValue;
  }
  numericFilterState[key]={min:minValue,max:maxValue};
  syncNumericFilterControl(key);
  applyReviewFilter();
}

function resetNumericFilter(key){
  const meta=numericFilterMeta[key];
  if(!meta||!meta.enabled) return;
  numericFilterState[key]={min:meta.scaleMin,max:meta.scaleMax};
  syncNumericFilterControl(key);
  applyReviewFilter();
}

function buildNumericFilterControls(){
  const grid=document.getElementById('numeric-filter-grid');
  if(!grid) return;
  grid.innerHTML='';
  NUMERIC_FILTER_CONFIGS.forEach(cfg=>{
    const meta=numericFilterMeta[cfg.key];
    const box=document.createElement('div');
    box.className='range-filter'+((!meta||!meta.enabled)?' is-disabled':'');
    box.setAttribute('data-filter-key', cfg.key);
    if(!meta||!meta.enabled){
      box.innerHTML=`<div class="range-header"><span class="range-label">${escapeHtml(cfg.label)}</span><span class="range-readout">No data</span></div>`;
      grid.appendChild(box);
      return;
    }
    box.innerHTML=`
      <div class="range-header">
        <span class="range-label">${escapeHtml(cfg.label)}</span>
        <span class="range-readout"></span>
      </div>
      <div class="range-sliders">
        <label>Min
          <input type="range" data-role="min-range" min="${meta.scaleMin}" max="${meta.scaleMax}" step="${meta.step}" value="${meta.scaleMin}">
        </label>
        <label>Max
          <input type="range" data-role="max-range" min="${meta.scaleMin}" max="${meta.scaleMax}" step="${meta.step}" value="${meta.scaleMax}">
        </label>
      </div>
      <div class="range-values">
        <label>Min value
          <input type="number" data-role="min-input" min="${meta.scaleMin}" max="${meta.scaleMax}" step="${meta.step}" value="${formatNumericValue(meta.scaleMin, meta.digits)}">
        </label>
        <label>Max value
          <input type="number" data-role="max-input" min="${meta.scaleMin}" max="${meta.scaleMax}" step="${meta.step}" value="${formatNumericValue(meta.scaleMax, meta.digits)}">
        </label>
        <button type="button" class="range-reset" onclick="resetNumericFilter('${cfg.key}')">Reset</button>
      </div>`;
    const minRange=box.querySelector('input[data-role="min-range"]');
    const maxRange=box.querySelector('input[data-role="max-range"]');
    const minInput=box.querySelector('input[data-role="min-input"]');
    const maxInput=box.querySelector('input[data-role="max-input"]');
    minRange.addEventListener('input', event=>setNumericFilterState(cfg.key, 'min', event.target.value));
    maxRange.addEventListener('input', event=>setNumericFilterState(cfg.key, 'max', event.target.value));
    minInput.addEventListener('change', event=>setNumericFilterState(cfg.key, 'min', event.target.value));
    maxInput.addEventListener('change', event=>setNumericFilterState(cfg.key, 'max', event.target.value));
    grid.appendChild(box);
    syncNumericFilterControl(cfg.key);
  });
}

function resetDashboardFilters(){
  ['review-filter','celltype-filter','svtype-filter','hard-conflict-filter'].forEach(id=>{
    const el=document.getElementById(id);
    if(el) el.value='all';
  });
  NUMERIC_FILTER_CONFIGS.forEach(cfg=>{
    const meta=numericFilterMeta[cfg.key];
    if(!meta||!meta.enabled) return;
    numericFilterState[cfg.key]={min:meta.scaleMin,max:meta.scaleMax};
    syncNumericFilterControl(cfg.key);
  });
  applyReviewFilter();
}

function loadDefaultReviewState(){
  const out={};
  getSvCards().forEach(card=>{
    const svId=card.getAttribute('data-sv-id')||'';
    out[svId]=normalizeReviewStatus(card.getAttribute('data-review-status'));
  });
  return out;
}

function hasLocalStorage(){
  try{
    if(typeof window==='undefined'||!window.localStorage) return false;
    const probe='__sniffcell_localstorage_probe__';
    window.localStorage.setItem(probe,'1');
    window.localStorage.removeItem(probe);
    return true;
  }catch(err){
    return false;
  }
}

function loadStoredReviewState(){
  if(!hasLocalStorage()){
    return {state:{},loaded:0,available:false,error:''};
  }
  let parsed={};
  try{
    const raw=window.localStorage.getItem(REVIEW_STORAGE_KEY);
    if(raw) parsed=JSON.parse(raw);
  }catch(err){
    return {state:{},loaded:0,available:true,error:String(err||'unknown error')};
  }
  if(!parsed||typeof parsed!=='object'||Array.isArray(parsed)){
    return {state:{},loaded:0,available:true,error:''};
  }
  const out={};
  let loaded=0;
  Object.entries(parsed).forEach(([svId,status])=>{
    const id=String(svId||'').trim();
    if(!id) return;
    out[id]=normalizeReviewStatus(status);
    loaded+=1;
  });
  return {state:out,loaded:loaded,available:true,error:''};
}

function setPersistStatus(message, isError){
  const el=document.getElementById('review-persist-status');
  if(!el) return;
  el.textContent=String(message||'');
  el.style.color=isError?'#9b1c1c':'#405160';
}

function persistReviewState(){
  if(!hasLocalStorage()) return false;
  try{
    window.localStorage.setItem(REVIEW_STORAGE_KEY, JSON.stringify(reviewState));
    return true;
  }catch(err){
    setPersistStatus(`Could not save to browser localStorage (${String(err||'unknown error')}).`, true);
    return false;
  }
}

function applyReviewToCard(card, status){
  card.setAttribute('data-review-status', status);
  card.classList.toggle('review-state-real', status==='real');
  card.classList.toggle('review-state-not-real', status==='not_real');
  card.classList.toggle('review-state-undecided', status==='undecided');
  const badge=card.querySelector('.review-badge');
  if(badge) badge.textContent=reviewLabel(status);
  card.querySelectorAll('.review-btn').forEach(btn=>{
    const active=btn.getAttribute('data-review-value')===status;
    btn.classList.toggle('is-active', active);
  });
}

function populateCelltypeFilterOptions(){
  const sel=document.getElementById('celltype-filter');
  if(!sel) return;
  const oldValue=String(sel.value||'all');
  const labelByToken=new Map();
  getSvCards().forEach(card=>{
    splitCelltypeList(card.getAttribute('data-primary-celltype')||'').forEach(label=>{
      const token=normalizeCelltypeToken(label);
      if(token&&!labelByToken.has(token)) labelByToken.set(token, label);
    });
    splitCelltypeList(card.getAttribute('data-linked-celltypes')||'').forEach(label=>{
      const token=normalizeCelltypeToken(label);
      if(token&&!labelByToken.has(token)) labelByToken.set(token, label);
    });
  });
  sel.innerHTML='';
  const allOpt=document.createElement('option');
  allOpt.value='all';
  allOpt.textContent='All assigned cell types';
  sel.appendChild(allOpt);
  Array.from(labelByToken.entries())
    .sort((a,b)=>String(a[1]).localeCompare(String(b[1])))
    .forEach(([token,label])=>{
      const opt=document.createElement('option');
      opt.value=token;
      opt.textContent=label;
      sel.appendChild(opt);
    });
  if(Array.from(sel.options).some(opt=>opt.value===oldValue)) sel.value=oldValue;
  else sel.value='all';
}

function populateSvtypeFilterOptions(){
  const sel=document.getElementById('svtype-filter');
  if(!sel) return;
  const oldValue=String(sel.value||'all');
  const values=new Set();
  getSvCards().forEach(card=>{
    const svType=String(card.getAttribute('data-sv-type')||'').trim().toUpperCase();
    if(svType&&svType!=='NA') values.add(svType);
  });
  sel.innerHTML='';
  const allOpt=document.createElement('option');
  allOpt.value='all';
  allOpt.textContent='All SV types';
  sel.appendChild(allOpt);
  Array.from(values).sort((a,b)=>a.localeCompare(b)).forEach(value=>{
    const opt=document.createElement('option');
    opt.value=value;
    opt.textContent=value;
    sel.appendChild(opt);
  });
  if(Array.from(sel.options).some(opt=>opt.value===oldValue)) sel.value=oldValue;
  else sel.value='all';
}

function syncReviewSummary(filterState){
  const cards=getSvCards();
  const counts={real:0,not_real:0,undecided:0};
  let visible=0;
  cards.forEach(card=>{
    const status=normalizeReviewStatus(card.getAttribute('data-review-status'));
    counts[status]=(counts[status]||0)+1;
    if(card.style.display!=='none') visible+=1;
  });
  const el=document.getElementById('review-summary');
  if(el){
    const scope=summaryScopeLabel(filterState);
    const scopeText=scope==='all' ? 'filters: none' : `filters: ${scope}`;
    el.textContent=`real: ${counts.real} | not real: ${counts.not_real} | undecided: ${counts.undecided} | visible: ${visible}/${cards.length} | ${scopeText}`;
  }
}

function applyReviewFilter(){
  const filterState=currentDashboardFilterState();
  getSvCards().forEach(card=>{
    const status=normalizeReviewStatus(card.getAttribute('data-review-status'));
    const showReview=(filterState.review==='all'||filterState.review===status);
    const showCelltype=cardMatchesCelltype(card, filterState.celltype);
    const showSvType=cardMatchesSvType(card, filterState.svtype);
    const showHardConflict=cardMatchesHardConflict(card, filterState.hardConflict);
    const showNumeric=cardMatchesNumericFilters(card);
    card.style.display=(showReview&&showCelltype&&showSvType&&showHardConflict&&showNumeric)?'':'none';
  });
  syncReviewSummary(filterState);
  renderSummaries(filterState);
}

function setSvReview(btn){
  const card=btn.closest('.sv[data-sv-id]');
  if(!card) return;
  const svId=card.getAttribute('data-sv-id')||'';
  const status=normalizeReviewStatus(btn.getAttribute('data-review-value'));
  reviewState[svId]=status;
  applyReviewToCard(card, status);
  if(persistReviewState()){
    setPersistStatus('Saved review labels to browser localStorage.', false);
  }
  applyReviewFilter();
}

function initReviewControls(){
  const defaults=loadDefaultReviewState();
  const stored=loadStoredReviewState();
  const validIds=new Set(Object.keys(defaults));
  reviewState={...defaults};

  let applied=0;
  Object.entries(stored.state).forEach(([svId,status])=>{
    if(!validIds.has(svId)) return;
    reviewState[svId]=normalizeReviewStatus(status);
    applied+=1;
  });

  if(!stored.available){
    setPersistStatus('Browser localStorage is unavailable; review labels are session-only.', true);
  } else if(stored.error){
    setPersistStatus(`Could not load browser localStorage state (${stored.error}).`, true);
  } else if(applied>0){
    setPersistStatus(`Loaded ${applied} review labels from browser localStorage.`, false);
  } else {
    setPersistStatus('No saved browser review labels found for this report.', false);
  }

  getSvCards().forEach(card=>{
    const svId=card.getAttribute('data-sv-id')||'';
    const status=normalizeReviewStatus(reviewState[svId]);
    reviewState[svId]=status;
    applyReviewToCard(card, status);
  });
  numericFilterMeta=buildNumericFilterMeta();
  NUMERIC_FILTER_CONFIGS.forEach(cfg=>{
    const meta=numericFilterMeta[cfg.key];
    if(meta&&meta.enabled){
      numericFilterState[cfg.key]={min:meta.scaleMin,max:meta.scaleMax};
    }
  });
  populateCelltypeFilterOptions();
  populateSvtypeFilterOptions();
  buildNumericFilterControls();
  applyReviewFilter();
  if(stored.available){
    persistReviewState();
  }
}

initReviewControls();

function copyVizCommand(btn){
  const cmd=btn.getAttribute('data-cmd')||'';
  const done=()=>{
    const old=btn.textContent;
    btn.textContent='Copied';
    setTimeout(()=>{btn.textContent=old;},1200);
  };
  if(navigator.clipboard&&window.isSecureContext){
    navigator.clipboard.writeText(cmd).then(done).catch(()=>{});
    return;
  }
  const ta=document.createElement('textarea');
  ta.value=cmd;
  document.body.appendChild(ta);
  ta.select();
  try{document.execCommand('copy');done();}catch(e){}
  document.body.removeChild(ta);
}
    """
        .replace("__DASHBOARD_DATA__", json.dumps(dashboard_records))
        .replace("__DASHBOARD_FILTERS__", json.dumps(filters))
        .replace("__REVIEW_STORAGE_KEY__", json.dumps(review_storage_key))
    )
    page.append(f"<script>{script_text}</script>")
    page.append("</body>")
    page.append("</html>")
    return "\n".join(page)


def report_main(args) -> None:
    logger = logging.getLogger("sniffcell.report")
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    anno_output = Path(args.anno_output)
    if not anno_output.exists():
        raise FileNotFoundError(f"anno_output does not exist: {anno_output}")

    manifest_path = anno_output / "anno_run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Could not find {manifest_path}. "
            "sniffcell report needs an anno output folder generated by sniffcell anno."
        )

    sv_assignment_path = anno_output / "sv_assignment.tsv"
    if not sv_assignment_path.exists():
        raise FileNotFoundError(f"Could not find sv_assignment.tsv: {sv_assignment_path}")

    if int(args.max_reads) <= 0:
        raise ValueError("max_reads must be > 0")
    if int(args.window) < 0:
        raise ValueError("window must be >= 0")
    if int(args.figure_threads) <= 0:
        raise ValueError("figure_threads must be > 0")
    if int(getattr(args, "igv_snapshot_width", 3600)) <= 0:
        raise ValueError("igv_snapshot_width must be > 0")
    if int(getattr(args, "igv_snapshot_height", 1600)) <= 0:
        raise ValueError("igv_snapshot_height must be > 0")

    figure_profile = str(getattr(args, "figure_profile", "full")).strip().lower()
    if figure_profile not in {"fast", "full"}:
        raise ValueError("figure_profile must be one of: fast, full")
    requested_dpi = int(getattr(args, "figure_dpi", 160))
    if requested_dpi <= 0:
        raise ValueError("figure_dpi must be > 0")

    requested_max_reads = int(args.max_reads)
    effective_max_reads = requested_max_reads
    effective_dpi = requested_dpi
    skip_methylation_overlay = (figure_profile == "fast")
    shared_threads = int(args.figure_threads)
    with_igvviz = bool(getattr(args, "with_igvviz", False))
    igv_cmd = str(getattr(args, "igv_cmd", "igv.sh"))
    igv_visibility_window = int(args.window)
    igv_phase_tag = "HP"
    igv_support_tag = "SC"
    igv_snapshot_format = str(getattr(args, "igv_snapshot_format", "png"))
    igv_snapshot_width = int(getattr(args, "igv_snapshot_width", 3600))
    igv_snapshot_height = int(getattr(args, "igv_snapshot_height", 1600))
    reuse_existing_igvviz = bool(getattr(args, "reuse_existing_igvviz", False))
    igv_bams = igvviz_module._split_bam_args(getattr(args, "igv_bams", None))
    with_igvreport = bool(getattr(args, "with_igvreport", False))

    # Fast profile intentionally trades panel detail for speed.
    if figure_profile == "fast":
        effective_max_reads = min(requested_max_reads, 120)
        effective_dpi = min(requested_dpi, 160)

    report_dir, figure_dir, html_path, archive_path = _resolve_report_paths(anno_output, args.output)
    report_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    igvviz_root = report_dir / "igvviz"
    igvviz_root.mkdir(parents=True, exist_ok=True)
    igvreport_root = report_dir / "igvreport"
    review_storage_key = f"sniffcell_report_review::{str(anno_output.resolve())}"
    manifest_payload = _load_json_dict(manifest_path)

    sv_df = _load_sv_assignment(sv_assignment_path)
    sv_df = _backfill_sv_fields_from_manifest_vcf(sv_df, manifest_payload, logger)
    selected = _select_high_confidence_svs(
        sv_df,
        min_overlap_pct=float(args.min_overlap_pct),
        min_majority_pct=float(args.min_majority_pct),
        include_unassigned=bool(args.include_unassigned),
        allow_hard_conflict=bool(args.allow_hard_conflict),
        max_sv=int(args.max_sv),
    )
    logger.info(
        "Selected %d/%d SVs for report (min_overlap_pct=%.3f min_majority_pct=%.3f include_unassigned=%s allow_hard_conflict=%s max_sv=%d)",
        len(selected),
        len(sv_df),
        float(args.min_overlap_pct),
        float(args.min_majority_pct),
        bool(args.include_unassigned),
        bool(args.allow_hard_conflict),
        int(args.max_sv),
    )
    logger.info(
        "Figure profile: %s (threads=%d dpi=%d max_reads=%d skip_methylation_overlay=%s)",
        figure_profile,
        shared_threads,
        int(effective_dpi),
        int(effective_max_reads),
        bool(skip_methylation_overlay),
    )
    logger.info("Review persistence uses browser localStorage key: %s", review_storage_key)
    if requested_max_reads != effective_max_reads:
        logger.info(
            "Fast profile capped max_reads from %d to %d for faster rendering.",
            requested_max_reads,
            effective_max_reads,
        )
    if requested_dpi != effective_dpi:
        logger.info(
            "Fast profile capped figure_dpi from %d to %d for faster rendering.",
            requested_dpi,
            effective_dpi,
        )
    logger.info(
        "IGV rendering: with_igvviz=%s threads=%d igv_cmd=%s igv_bams=%s visibility_window=%d phase_tag=%s support_tag=%s include_non_supporting=%s keep_intermediates=%s snapshot_format=%s size=%dx%d",
        with_igvviz,
        shared_threads,
        igv_cmd,
        "|".join(igv_bams) if igv_bams else "<anno-manifest-bam>",
        igv_visibility_window,
        igv_phase_tag,
        igv_support_tag,
        True,
        True,
        igv_snapshot_format,
        igv_snapshot_width,
        igv_snapshot_height,
    )
    logger.info(
        "Alternate IGV report: with_igvreport=%s flanking=%d",
        with_igvreport,
        int(args.window),
    )
    if bool(args.with_figures) and len(selected) >= 25:
        logger.warning(
            "Report will render %d SV panels; this can take a long time. "
            "Use --max_sv to limit candidates, --reuse_existing_viz to skip existing plots, "
            "or increase --figure_threads.",
            len(selected),
        )
    if with_igvviz and len(selected) >= 10:
        logger.warning(
            "Report will run igvviz for %d SVs across %d BAM(s); this can take a long time. "
            "Use --max_sv to limit candidates or increase --figure_threads.",
            len(selected),
            max(1, len(igv_bams)),
        )
    if (not bool(args.with_figures)) and len(selected) > 0:
        logger.info(
            "Figure rendering is disabled (default figure-less mode). "
            "Use --with_figures to render panels, or copy per-SV viz commands from the report HTML."
        )
    if (not with_igvviz) and len(selected) > 0:
        logger.info(
            "IGV rendering is disabled. Use --with_igvviz to generate igvviz screenshots per selected SV."
        )
    if (not with_igvreport) and len(selected) > 0:
        logger.info(
            "Alternate IGV.js report is disabled. Use --with_igvreport to generate an igv-reports HTML page for the selected SVs."
        )

    rows_for_report: list[dict[str, object]] = []
    slug_counts: dict[str, int] = {}
    render_jobs: list[tuple[int, str, Path]] = []
    igvviz_jobs: list[tuple[int, str, Path, Path]] = []

    for row_idx, row in enumerate(selected.to_dict(orient="records")):
        sv_id = str(row["id"])
        base_slug = _safe_slug(sv_id)
        slug_counts[base_slug] = slug_counts.get(base_slug, 0) + 1
        slug = base_slug if slug_counts[base_slug] == 1 else f"{base_slug}_{slug_counts[base_slug]}"
        figure_path = figure_dir / f"{slug}.viz.{args.format}"
        igvviz_dir = igvviz_root / slug
        sv_slug_for_igv = _safe_slug(sv_id)
        igvviz_manifest_path = igvviz_dir / f"{sv_slug_for_igv}.igvviz.manifest.json"
        viz_command = _build_viz_cli_command(
            anno_output=anno_output,
            sv_id=sv_id,
            output_path=figure_path,
            window=int(args.window),
            max_reads=int(effective_max_reads),
            fmt=str(args.format),
            dpi=int(effective_dpi),
            exact_window=True,
            skip_methylation_overlay=bool(skip_methylation_overlay),
        )
        igvviz_command = _build_igvviz_cli_command(
            anno_output=anno_output,
            sv_id=sv_id,
            output_dir=igvviz_dir,
            window=int(args.window),
            igv_bams=igv_bams,
            igv_cmd=igv_cmd,
            snapshot_format=igv_snapshot_format,
            snapshot_width=igv_snapshot_width,
            snapshot_height=igv_snapshot_height,
        )

        report_row = dict(row)
        report_row["review_status"] = _normalize_review_status(report_row.get("review_status", "undecided"))
        report_row["viz_status"] = "not_rendered"
        report_row["viz_error"] = ""
        report_row["viz_figure"] = str(figure_path)
        report_row["viz_figure_rel"] = ""
        report_row["viz_command"] = viz_command
        report_row["igvviz_status"] = "not_rendered"
        report_row["igvviz_error"] = ""
        report_row["igvviz_dir"] = str(igvviz_dir)
        report_row["igvviz_dir_rel"] = ""
        report_row["igvviz_manifest"] = str(igvviz_manifest_path)
        report_row["igvviz_manifest_rel"] = ""
        report_row["igvviz_snapshots_rel"] = ""
        report_row["igvviz_snapshot_bams"] = ""
        report_row["_igvviz_snapshots"] = []
        report_row["igvviz_command"] = igvviz_command
        rows_for_report.append(report_row)
        render_jobs.append((row_idx, sv_id, figure_path))
        igvviz_jobs.append((row_idx, sv_id, igvviz_dir, igvviz_manifest_path))

    if bool(args.with_figures) and render_jobs:
        n_threads = int(args.figure_threads)
        progress = None
        if tqdm is not None:
            progress = tqdm(total=len(render_jobs), desc="Rendering SV figures", unit="sv")
        if n_threads <= 1 or len(render_jobs) <= 1:
            try:
                for row_idx, sv_id, figure_path in render_jobs:
                    status, error = _render_one_viz_panel(
                        anno_output=anno_output,
                        sv_id=sv_id,
                        figure_path=figure_path,
                        window=int(args.window),
                        max_reads=int(effective_max_reads),
                        fmt=str(args.format),
                        dpi=int(effective_dpi),
                        exact_window=True,
                        skip_methylation_overlay=bool(skip_methylation_overlay),
                        reuse_existing_viz=bool(args.reuse_existing_viz),
                    )
                    rows_for_report[row_idx]["viz_status"] = status
                    rows_for_report[row_idx]["viz_error"] = error
                    if progress is not None:
                        progress.update(1)
            finally:
                if progress is not None:
                    progress.close()
        else:
            try:
                with ThreadPoolExecutor(max_workers=n_threads) as executor:
                    future_map = {
                        executor.submit(
                            _render_one_viz_panel,
                            anno_output=anno_output,
                            sv_id=sv_id,
                            figure_path=figure_path,
                            window=int(args.window),
                            max_reads=int(effective_max_reads),
                            fmt=str(args.format),
                            dpi=int(effective_dpi),
                            exact_window=True,
                            skip_methylation_overlay=bool(skip_methylation_overlay),
                            reuse_existing_viz=bool(args.reuse_existing_viz),
                        ): (row_idx, sv_id)
                        for row_idx, sv_id, figure_path in render_jobs
                    }
                    for future in as_completed(future_map):
                        row_idx, sv_id = future_map[future]
                        try:
                            status, error = future.result()
                        except Exception as exc:
                            status, error = "failed", str(exc)
                            logger.exception("Unexpected viz worker failure for SV %s", sv_id)
                        rows_for_report[row_idx]["viz_status"] = status
                        rows_for_report[row_idx]["viz_error"] = error
                        if str(status).startswith("failed"):
                            logger.error("viz failed for SV %s: %s", sv_id, error)
                        if progress is not None:
                            progress.update(1)
            finally:
                if progress is not None:
                    progress.close()

    if with_igvviz and igvviz_jobs:
        n_threads = shared_threads
        progress = None
        if tqdm is not None:
            progress = tqdm(total=len(igvviz_jobs), desc="Rendering IGV screenshots", unit="sv")
        if n_threads <= 1 or len(igvviz_jobs) <= 1:
            try:
                for row_idx, sv_id, igvviz_dir, igvviz_manifest_path in igvviz_jobs:
                    status, error = _render_one_igvviz_bundle(
                        anno_output=anno_output,
                        sv_id=sv_id,
                        output_dir=igvviz_dir,
                        expected_manifest_path=igvviz_manifest_path,
                        window=int(args.window),
                        igv_bams=igv_bams,
                        igv_cmd=igv_cmd,
                        snapshot_format=igv_snapshot_format,
                        snapshot_width=igv_snapshot_width,
                        snapshot_height=igv_snapshot_height,
                        reuse_existing_igvviz=bool(reuse_existing_igvviz),
                    )
                    rows_for_report[row_idx]["igvviz_status"] = status
                    rows_for_report[row_idx]["igvviz_error"] = error
                    if progress is not None:
                        progress.update(1)
            finally:
                if progress is not None:
                    progress.close()
        else:
            try:
                with ThreadPoolExecutor(max_workers=n_threads) as executor:
                    future_map = {
                        executor.submit(
                            _render_one_igvviz_bundle,
                            anno_output=anno_output,
                            sv_id=sv_id,
                            output_dir=igvviz_dir,
                            expected_manifest_path=igvviz_manifest_path,
                            window=int(args.window),
                            igv_bams=igv_bams,
                            igv_cmd=igv_cmd,
                            snapshot_format=igv_snapshot_format,
                            snapshot_width=igv_snapshot_width,
                            snapshot_height=igv_snapshot_height,
                            reuse_existing_igvviz=bool(reuse_existing_igvviz),
                        ): (row_idx, sv_id)
                        for row_idx, sv_id, igvviz_dir, igvviz_manifest_path in igvviz_jobs
                    }
                    for future in as_completed(future_map):
                        row_idx, sv_id = future_map[future]
                        try:
                            status, error = future.result()
                        except Exception as exc:
                            status, error = "failed", str(exc)
                            logger.exception("Unexpected igvviz worker failure for SV %s", sv_id)
                        rows_for_report[row_idx]["igvviz_status"] = status
                        rows_for_report[row_idx]["igvviz_error"] = error
                        if str(status).startswith("failed"):
                            logger.error("igvviz failed for SV %s: %s", sv_id, error)
                        if progress is not None:
                            progress.update(1)
            finally:
                if progress is not None:
                    progress.close()

    igvreport_html_path = igvreport_root / "index.html"
    igvreport_manifest_path = igvreport_root / "igvreport_manifest.json"
    igvreport_result: dict[str, object] = {
        "enabled": bool(with_igvreport),
        "status": "not_rendered",
        "error": "",
        "command": "",
        "html": str(igvreport_html_path),
        "html_rel": "",
        "manifest": str(igvreport_manifest_path),
        "manifest_rel": "",
    }
    if with_igvreport and len(selected) > 0:
        igvreport_result.update(
            igvreport_module.render_igvreport_bundle(
                anno_output=anno_output,
                selected_df=selected.copy(),
                output_dir=igvreport_root,
                native_report_html=html_path,
                igv_bams=igv_bams,
                window=int(args.window),
            )
        )
    elif (not with_igvreport) and igvreport_html_path.exists() and igvreport_manifest_path.exists():
        igvreport_result["status"] = "existing"

    igvreport_payload = _load_json_dict(igvreport_manifest_path)
    if (not str(igvreport_result.get("command", "")).strip()) and igvreport_payload:
        igvreport_result["command"] = str(igvreport_payload.get("igvreport_command", ""))
    if (not str(igvreport_result.get("error", "")).strip()) and igvreport_payload:
        igvreport_result["error"] = str(igvreport_payload.get("error", ""))
    if (not str(igvreport_result.get("status", "")).strip()) and igvreport_payload:
        igvreport_result["status"] = str(igvreport_payload.get("status", ""))

    rendered_count = 0
    failed_count = 0
    igvviz_rendered_count = 0
    igvviz_failed_count = 0
    igvreport_rendered_count = 0
    igvreport_failed_count = 0
    for report_row in rows_for_report:
        fig_path = Path(str(report_row["viz_figure"]))
        if (not bool(args.with_figures)) and fig_path.exists():
            report_row["viz_status"] = "existing"
        if fig_path.exists():
            rendered_count += 1
            report_row["viz_figure_rel"] = fig_path.relative_to(html_path.parent).as_posix()
        if str(report_row["viz_status"]).startswith("failed"):
            failed_count += 1

        igv_manifest_path = Path(str(report_row["igvviz_manifest"]))
        igv_dir_path = Path(str(report_row["igvviz_dir"]))
        if (not with_igvviz) and igv_manifest_path.exists():
            report_row["igvviz_status"] = "existing"
        if igv_dir_path.exists():
            report_row["igvviz_dir_rel"] = igv_dir_path.relative_to(html_path.parent).as_posix()
        if igv_manifest_path.exists():
            igvviz_rendered_count += 1
            report_row["igvviz_manifest_rel"] = igv_manifest_path.relative_to(html_path.parent).as_posix()
            snap_rows = _load_igvviz_snapshot_rows(igv_manifest_path, html_path.parent)
            report_row["_igvviz_snapshots"] = snap_rows
            rels = [str(x.get("snapshot_rel", "")).strip() for x in snap_rows if str(x.get("snapshot_rel", "")).strip()]
            bam_labels = [str(x.get("bam_label", "")).strip() for x in snap_rows if str(x.get("snapshot_rel", "")).strip()]
            report_row["igvviz_snapshots_rel"] = ";".join(rels)
            report_row["igvviz_snapshot_bams"] = ";".join(bam_labels)
        if str(report_row["igvviz_status"]).startswith("failed"):
            igvviz_failed_count += 1

    if igvreport_html_path.exists():
        igvreport_rendered_count = 1
        igvreport_result["html_rel"] = igvreport_html_path.relative_to(html_path.parent).as_posix()
    if igvreport_manifest_path.exists():
        igvreport_result["manifest_rel"] = igvreport_manifest_path.relative_to(html_path.parent).as_posix()
    if str(igvreport_result.get("status", "")).startswith("failed"):
        igvreport_failed_count = 1

    selected_report_df = pd.DataFrame(rows_for_report)
    if "review_status" not in selected_report_df.columns:
        selected_report_df["review_status"] = pd.Series(dtype="string")
    selected_report_df["review_status"] = selected_report_df["review_status"].map(_normalize_review_status).astype("string")
    for col in (
        "viz_status",
        "viz_error",
        "viz_figure",
        "viz_figure_rel",
        "viz_command",
        "igvviz_status",
        "igvviz_error",
        "igvviz_dir",
        "igvviz_dir_rel",
        "igvviz_manifest",
        "igvviz_manifest_rel",
        "igvviz_snapshots_rel",
        "igvviz_snapshot_bams",
        "igvviz_command",
    ):
        if col not in selected_report_df.columns:
            selected_report_df[col] = pd.Series(dtype="string")

    selected_tsv_path = report_dir / "high_confidence_sv.tsv"
    selected_report_df.to_csv(selected_tsv_path, sep="\t", index=False)

    failed_tsv_path = report_dir / "failed_viz.tsv"
    failed_only = selected_report_df[selected_report_df["viz_status"].astype(str).str.startswith("failed")].copy()
    failed_only.to_csv(failed_tsv_path, sep="\t", index=False)
    failed_igvviz_tsv_path = report_dir / "failed_igvviz.tsv"
    failed_igv_only = selected_report_df[selected_report_df["igvviz_status"].astype(str).str.startswith("failed")].copy()
    failed_igv_only.to_csv(failed_igvviz_tsv_path, sep="\t", index=False)

    filters = {
        "min_overlap_pct": float(args.min_overlap_pct),
        "min_majority_pct": float(args.min_majority_pct),
        "include_unassigned": bool(args.include_unassigned),
        "allow_hard_conflict": bool(args.allow_hard_conflict),
        "max_sv": int(args.max_sv),
    }
    viz_cfg = {
        "with_figures": bool(args.with_figures),
        "figure_threads": int(shared_threads),
        "figure_profile": str(figure_profile),
        "figure_dpi": int(effective_dpi),
        "skip_methylation_overlay": bool(skip_methylation_overlay),
        "exact_window": True,
        "window": int(args.window),
        "max_reads": int(effective_max_reads),
        "format": str(args.format),
        "export_tables": True,
        "reuse_existing_viz": bool(args.reuse_existing_viz),
    }
    igvviz_cfg = {
        "with_igvviz": bool(with_igvviz),
        "igv_bams": list(igv_bams),
        "igv_threads": int(shared_threads),
        "igv_cmd": str(igv_cmd),
        "igv_visibility_window": int(igv_visibility_window),
        "igv_phase_tag": str(igv_phase_tag),
        "igv_support_tag": str(igv_support_tag),
        "igv_include_non_supporting": True,
        "igv_snapshot_format": str(igv_snapshot_format),
        "igv_snapshot_width": int(igv_snapshot_width),
        "igv_snapshot_height": int(igv_snapshot_height),
        "reuse_existing_igvviz": bool(reuse_existing_igvviz),
        "igv_keep_intermediates": True,
        "igv_batch_only": False,
    }
    igvreport_cfg = {
        "with_igvreport": bool(with_igvreport),
        "igvreport_flanking": int(args.window),
    }
    dashboard_records = _build_dashboard_records(selected_report_df)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_text = _build_report_html(
        generated_at=generated_at,
        anno_output=anno_output,
        sv_assignment_path=sv_assignment_path,
        review_storage_key=review_storage_key,
        filters=filters,
        viz=viz_cfg,
        igvreport=igvreport_result,
        total_sv=int(len(sv_df)),
        selected_count=int(len(selected)),
        rendered_count=int(rendered_count),
        failed_count=int(failed_count),
        rows=rows_for_report,
        dashboard_records=dashboard_records,
    )
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(html_text, encoding="utf-8")

    manifest_payload = {
        "command": "report",
        "generated_at": generated_at,
        "inputs": {
            "anno_output": str(anno_output.resolve()),
            "sv_assignment": str(sv_assignment_path.resolve()),
        },
        "filters": filters,
        "viz": viz_cfg,
        "igvviz": igvviz_cfg,
        "igvreport": {
            **igvreport_cfg,
            "status": str(igvreport_result.get("status", "")),
            "error": str(igvreport_result.get("error", "")),
            "command": str(igvreport_result.get("command", "")),
        },
        "counts": {
            "sv_total": int(len(sv_df)),
            "sv_selected": int(len(selected)),
            "viz_rendered_or_reused": int(rendered_count),
            "viz_failed": int(failed_count),
            "igvviz_rendered_or_reused": int(igvviz_rendered_count),
            "igvviz_failed": int(igvviz_failed_count),
            "igvreport_rendered_or_reused": int(igvreport_rendered_count),
            "igvreport_failed": int(igvreport_failed_count),
        },
        "outputs": {
            "report_dir": str(report_dir.resolve()),
            "html": str(html_path.resolve()),
            "figures_dir": str(figure_dir.resolve()),
            "igvviz_dir": str(igvviz_root.resolve()),
            "igvreport_dir": str(igvreport_root.resolve()),
            "igvreport_html": str(igvreport_html_path.resolve()),
            "igvreport_manifest": str(igvreport_manifest_path.resolve()),
            "high_confidence_sv_tsv": str(selected_tsv_path.resolve()),
            "review_storage_key": review_storage_key,
            "failed_viz_tsv": str(failed_tsv_path.resolve()),
            "failed_igvviz_tsv": str(failed_igvviz_tsv_path.resolve()),
            "report_archive": (str(archive_path.resolve()) if archive_path is not None else ""),
        },
    }
    manifest_out = report_dir / "report_manifest.json"
    manifest_out.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True), encoding="utf-8")

    archive_out: Path | None = None
    if archive_path is not None:
        archive_out = _write_report_archive(report_dir, archive_path)

    logger.info("Wrote HTML report: %s", html_path)
    logger.info("Wrote selected SV table: %s", selected_tsv_path)
    logger.info("Review state localStorage key: %s", review_storage_key)
    logger.info("Wrote failed viz table: %s", failed_tsv_path)
    logger.info("Wrote failed igvviz table: %s", failed_igvviz_tsv_path)
    if str(igvreport_result.get("status", "")).strip() not in {"", "not_rendered"}:
        logger.info(
            "IGV alternate report status=%s html=%s manifest=%s",
            igvreport_result.get("status", ""),
            igvreport_html_path,
            igvreport_manifest_path,
        )
    logger.info("Wrote report manifest: %s", manifest_out)
    if archive_out is not None:
        logger.info("Wrote gzipped report archive: %s", archive_out)
