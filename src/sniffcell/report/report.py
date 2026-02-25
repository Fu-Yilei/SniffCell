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

from sniffcell.viz import viz as viz_module
from sniffcell.viz import igvviz as igvviz_module

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


def _build_dashboard_records(selected_report_df: pd.DataFrame) -> list[dict[str, object]]:
    if selected_report_df.empty:
        return []

    cols = [
        "id",
        "sv_chr",
        "sv_pos",
        "sv_len",
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

    if "has_hard_conflict" not in out.columns:
        out["has_hard_conflict"] = pd.Series(pd.array([pd.NA] * len(out), dtype="boolean"))
    else:
        out["has_hard_conflict"] = out["has_hard_conflict"].map(_parse_bool).astype("boolean")

    out["id"] = out["id"].astype("string")
    return out


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
    filters: dict[str, object],
    viz: dict[str, object],
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
        ".review-controls h2{margin:0 0 10px 0;font-size:20px;}"
        ".review-row{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-bottom:8px;}"
        ".review-select{border:1px solid #c8d1da;border-radius:6px;padding:6px 10px;background:#fff;color:#14212c;font-size:14px;}"
        ".review-export-btn{border:0;border-radius:6px;padding:7px 12px;background:#0a7f5a;color:#fff;cursor:pointer;font-size:13px;}"
        ".review-export-btn:hover{background:#076d4d;}"
        ".review-summary{font-size:13px;color:#405160;margin:4px 0 10px 0;}"
        ".review-export-wrap{overflow:auto;max-height:320px;border:1px solid #d8dee4;border-radius:8px;background:#fff;}"
        ".review-table{width:100%;border-collapse:collapse;font-size:13px;}"
        ".review-table th,.review-table td{padding:6px 8px;border-bottom:1px solid #e9eef2;text-align:left;white-space:nowrap;}"
        ".review-table th{position:sticky;top:0;background:#f8fafc;z-index:1;}"
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
        ".igv-grid{display:flex;flex-direction:column;gap:10px;margin-top:10px;}"
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
        f"exact_window={viz['exact_window']}"
        "</div>"
    )
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
    page.append(
        "<button id=\"export-review\" class=\"review-export-btn\" type=\"button\" "
        "onclick=\"exportReviewTable()\">Export review table</button>"
    )
    page.append("</div>")
    page.append("<div id=\"review-summary\" class=\"review-summary\"></div>")
    page.append("<div id=\"review-export-container\" class=\"review-export-wrap\"></div>")
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
            n_supporting = "NA" if pd.isna(item.get("n_supporting")) else int(item["n_supporting"])
            n_overlapped = "NA" if pd.isna(item.get("n_overlapped")) else int(item["n_overlapped"])
            status = html.escape(str(item.get("viz_status", "")))
            n_supporting_text = str(n_supporting)
            n_overlapped_text = str(n_overlapped)
            status_text = _fmt_text(item.get("viz_status", ""))
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

            page.append(
                "<section class=\"sv review-state-undecided\" "
                f"data-sv-id=\"{sv_id_attr}\" "
                "data-review-status=\"undecided\" "
                f"data-primary-celltype=\"{html.escape(primary_text, quote=True)}\" "
                f"data-linked-celltypes=\"{html.escape(linked_text, quote=True)}\" "
                f"data-assigned-code=\"{html.escape(assigned_code_text, quote=True)}\" "
                f"data-majority-pct=\"{html.escape(majority, quote=True)}\" "
                f"data-overlap-pct=\"{html.escape(overlap, quote=True)}\" "
                f"data-n-supporting=\"{html.escape(n_supporting_text, quote=True)}\" "
                f"data-n-overlapped=\"{html.escape(n_overlapped_text, quote=True)}\" "
                f"data-sv-len=\"{html.escape(sv_len_text, quote=True)}\" "
                f"data-viz-status=\"{html.escape(status_text, quote=True)}\">"
            )
            page.append(f"<h2>{sv_id}</h2>")
            page.append("<div class=\"review-buttons\">")
            page.append(
                "<button class=\"review-btn review-real\" type=\"button\" data-review-value=\"real\" "
                "onclick=\"setSvReview(this)\">Real</button>"
            )
            page.append(
                "<button class=\"review-btn review-not-real\" type=\"button\" data-review-value=\"not_real\" "
                "onclick=\"setSvReview(this)\">Not real</button>"
            )
            page.append(
                "<button class=\"review-btn review-undecided is-active\" type=\"button\" data-review-value=\"undecided\" "
                "onclick=\"setSvReview(this)\">Undecided</button>"
            )
            page.append("</div>")
            page.append("<div class=\"kv\"><b>Review:</b> <span class=\"review-badge\">undecided</span></div>")
            page.append(f"<div class=\"kv\"><b>Primary cell type:</b> {primary}</div>")
            page.append(f"<div class=\"kv\"><b>Linked cell types:</b> {linked}</div>")
            page.append(f"<div class=\"kv\"><b>Assigned code:</b> {assigned_code}</div>")
            page.append(f"<div class=\"kv\"><b>SV length:</b> {html.escape(sv_len_display)}</div>")
            page.append(f"<div class=\"kv\"><b>IGV SV locus:</b> <code>{html.escape(sv_igv_text)}</code></div>")
            page.append(f"<div class=\"kv\"><b>majority_pct:</b> {majority} | <b>overlap_pct:</b> {overlap}</div>")
            page.append(
                f"<div class=\"kv\"><b>n_supporting:</b> {html.escape(n_supporting_text)} | "
                f"<b>n_overlapped:</b> {html.escape(n_overlapped_text)}</div>"
            )
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
                    page.append("<div class=\"igv-grid\">")
                    for snap in shown:
                        snap_rel = html.escape(str(snap.get("snapshot_rel", "")).strip())
                        bam_label = html.escape(str(snap.get("bam_label", "IGV")).strip() or "IGV")
                        page.append("<div class=\"igv-card\">")
                        page.append(f"<div class=\"kv\"><b>{bam_label}</b></div>")
                        page.append(f"<img src=\"{snap_rel}\" alt=\"IGV snapshot {bam_label} for {sv_id}\" loading=\"lazy\">")
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
    page.append(
        "<script>"
        f"const dashboardData={json.dumps(dashboard_records)};"
        f"const dashboardFilters={json.dumps(filters)};"
        "function chrSortKey(chr){"
        "const t=String(chr||'').replace(/^chr/i,'');"
        "if(/^\\d+$/.test(t)) return [0, parseInt(t,10), t];"
        "if(t==='X') return [1,23,t];"
        "if(t==='Y') return [1,24,t];"
        "if(t==='M'||t==='MT') return [1,25,t];"
        "return [2, Number.MAX_SAFE_INTEGER, t];"
        "}"
        "function sortChromLabels(labels){"
        "return labels.slice().sort((a,b)=>{const ka=chrSortKey(a), kb=chrSortKey(b);"
        "if(ka[0]!==kb[0]) return ka[0]-kb[0]; if(ka[1]!==kb[1]) return ka[1]-kb[1];"
        "return String(ka[2]).localeCompare(String(kb[2]));});"
        "}"
        "function selectedReviewFilter(){const sel=document.getElementById('review-filter');return sel?String(sel.value||'all'):'all';}"
        "function filteredDashboardData(filterValue){"
        "const wanted=String(filterValue||'all');"
        "const base=Array.isArray(dashboardData)?dashboardData:[];"
        "if(wanted==='all'){return base;}"
        "return base.filter(r=>{"
        "const id=String((r&&r.id!=null)?r.id:'');"
        "const status=normalizeReviewStatus(reviewState[id]);"
        "return status===wanted;"
        "});"
        "}"
        "function summaryScopeLabel(filterValue){"
        "const wanted=String(filterValue||'all');"
        "if(wanted==='real') return 'real';"
        "if(wanted==='not_real') return 'not real';"
        "if(wanted==='undecided') return 'undecided';"
        "return 'all';"
        "}"
        "function renderSummaries(filterValue){"
        "if(typeof Plotly==='undefined'){"
        "const msg='Plotly failed to load; interactive plots unavailable.';"
        "['chart-genome-location','chart-chrom-counts','chart-svlen','chart-support','chart-overlap-majority','chart-celltype']"
        ".forEach(id=>{const el=document.getElementById(id); if(el){el.textContent=msg;}}); return;}"
        "const data=filteredDashboardData(filterValue);"
        "const scope=summaryScopeLabel(filterValue);"
        "const titleSuffix=(scope==='all')?'':` [${scope}]`;"
        "if(data.length===0){"
        "const msg=(scope==='all')?'No selected SVs for summary plots.':`No SVs in '${scope}' for summary plots.`;"
        "['chart-genome-location','chart-chrom-counts','chart-svlen','chart-support','chart-overlap-majority','chart-celltype']"
        ".forEach(id=>{const el=document.getElementById(id); if(el){el.textContent=msg;}}); return;}"
        "const posRows=data.filter(r=>r.sv_chr!=null&&r.sv_pos!=null&&Number.isFinite(Number(r.sv_pos)));"
        "if(posRows.length>0){"
        "const yMb=posRows.map(r=>Number(r.sv_pos)/1e6);"
        "const txt=posRows.map(r=>`${r.id||'NA'}<br>${r.sv_chr}:${r.sv_pos}`);"
        "Plotly.newPlot('chart-genome-location',[{type:'scatter',mode:'markers',x:posRows.map(r=>String(r.sv_chr)),y:yMb,text:txt,hovertemplate:'%{text}<br>Position(Mb): %{y:.3f}<extra></extra>',marker:{size:8,color:'#1f77b4',opacity:0.75}}],"
        "{title:'Genome-wide SV Locations (selected)'+titleSuffix,xaxis:{title:'Chromosome',categoryorder:'array',categoryarray:sortChromLabels(posRows.map(r=>String(r.sv_chr)))},yaxis:{title:'SV position (Mb)'}},{responsive:true,displaylogo:false});"
        "}else{document.getElementById('chart-genome-location').textContent='No sv_chr/sv_pos data.';}"
        "const chrCounts={};"
        "data.forEach(r=>{const c=(r.sv_chr==null||String(r.sv_chr).trim()==='')?'NA':String(r.sv_chr); chrCounts[c]=(chrCounts[c]||0)+1;});"
        "const chrLabels=sortChromLabels(Object.keys(chrCounts));"
        "Plotly.newPlot('chart-chrom-counts',[{type:'bar',x:chrLabels,y:chrLabels.map(c=>chrCounts[c]),marker:{color:'#2ca02c'}}],{title:'SV Count by Chromosome'+titleSuffix,xaxis:{title:'Chromosome'},yaxis:{title:'SV count'}},{responsive:true,displaylogo:false});"
        "const lenVals=data.map(r=>Math.abs(Number(r.sv_len))).filter(v=>Number.isFinite(v)&&v>0).map(v=>Math.log10(v));"
        "if(lenVals.length>0){"
        "Plotly.newPlot('chart-svlen',[{type:'histogram',x:lenVals,marker:{color:'#9467bd'}}],{title:'SV Length Distribution'+titleSuffix,xaxis:{title:'log10(|sv_len| bp)'},yaxis:{title:'Count'}},{responsive:true,displaylogo:false});"
        "}else{document.getElementById('chart-svlen').textContent='No sv_len data.';}"
        "const nSup=data.map(r=>Number(r.n_supporting)).filter(v=>Number.isFinite(v));"
        "const nOvl=data.map(r=>Number(r.n_overlapped)).filter(v=>Number.isFinite(v));"
        "if(nSup.length+nOvl.length>0){"
        "Plotly.newPlot('chart-support',[{type:'histogram',x:nSup,name:'n_supporting',opacity:0.65,marker:{color:'#ff7f0e'}},{type:'histogram',x:nOvl,name:'n_overlapped',opacity:0.65,marker:{color:'#17becf'}}],"
        "{title:'Read Support Distribution'+titleSuffix,xaxis:{title:'Read count'},yaxis:{title:'SV count'},barmode:'overlay'},{responsive:true,displaylogo:false});"
        "}else{document.getElementById('chart-support').textContent='No read support data.';}"
        "const omRows=data.filter(r=>Number.isFinite(Number(r.overlap_pct))&&Number.isFinite(Number(r.majority_pct)));"
        "if(omRows.length>0){"
        "Plotly.newPlot('chart-overlap-majority',[{type:'scatter',mode:'markers',x:omRows.map(r=>Number(r.overlap_pct)),y:omRows.map(r=>Number(r.majority_pct)),text:omRows.map(r=>String(r.id||'')),"
        "hovertemplate:'%{text}<br>overlap=%{x:.3f}<br>majority=%{y:.3f}<extra></extra>',marker:{size:9,color:omRows.map(r=>Number.isFinite(Number(r.n_supporting))?Number(r.n_supporting):0),colorscale:'Viridis',showscale:true,colorbar:{title:'n_supporting'}}}],"
        "{title:'Agreement vs Overlap'+titleSuffix,xaxis:{title:'overlap_pct',range:[0,1]},yaxis:{title:'majority_pct',range:[0,1]},"
        "shapes:[{type:'line',x0:Number(dashboardFilters.min_overlap_pct||0),x1:Number(dashboardFilters.min_overlap_pct||0),y0:0,y1:1,line:{dash:'dot',color:'#444'}},"
        "{type:'line',y0:Number(dashboardFilters.min_majority_pct||0),y1:Number(dashboardFilters.min_majority_pct||0),x0:0,x1:1,line:{dash:'dot',color:'#444'}}]},"
        "{responsive:true,displaylogo:false});"
        "}else{document.getElementById('chart-overlap-majority').textContent='No overlap/majority data.';}"
        "const ctCounts={};"
        "data.forEach(r=>{const ct=(r.primary_celltype==null||String(r.primary_celltype).trim()==='')?'NA':String(r.primary_celltype); ctCounts[ct]=(ctCounts[ct]||0)+1;});"
        "const ctPairs=Object.entries(ctCounts).sort((a,b)=>b[1]-a[1]).slice(0,15);"
        "Plotly.newPlot('chart-celltype',[{type:'bar',orientation:'h',x:ctPairs.map(p=>p[1]).reverse(),y:ctPairs.map(p=>p[0]).reverse(),marker:{color:'#8c564b'}}],"
        "{title:'Top Primary Cell Types'+titleSuffix,xaxis:{title:'SV count'},yaxis:{title:'Cell type'}},{responsive:true,displaylogo:false});"
        "}"
        "const REVIEW_STORAGE_KEY='sniffcell_review::'+(window.location.pathname||'report');"
        "const REVIEW_STATUS_LABELS={real:'real',not_real:'not real',undecided:'undecided'};"
        "let reviewState={};"
        "function getSvCards(){return Array.from(document.querySelectorAll('.sv[data-sv-id]'));}"
        "function normalizeReviewStatus(value){const v=String(value||'').toLowerCase();"
        "if(v==='real') return 'real'; if(v==='not_real'||v==='not-real'||v==='not real') return 'not_real'; return 'undecided';}"
        "function reviewLabel(status){return REVIEW_STATUS_LABELS[status]||'undecided';}"
        "function escapeHtml(text){const span=document.createElement('span');span.textContent=String(text==null?'':text);return span.innerHTML;}"
        "function loadReviewState(){"
        "try{const raw=window.localStorage.getItem(REVIEW_STORAGE_KEY);"
        "if(!raw){return {};}"
        "const parsed=JSON.parse(raw);"
        "if(parsed&&typeof parsed==='object'){return parsed;}"
        "}catch(e){}"
        "return {};"
        "}"
        "function saveReviewState(){"
        "try{window.localStorage.setItem(REVIEW_STORAGE_KEY, JSON.stringify(reviewState));}catch(e){}"
        "}"
        "function applyReviewToCard(card,status){"
        "card.setAttribute('data-review-status', status);"
        "card.classList.toggle('review-state-real', status==='real');"
        "card.classList.toggle('review-state-not-real', status==='not_real');"
        "card.classList.toggle('review-state-undecided', status==='undecided');"
        "const badge=card.querySelector('.review-badge');"
        "if(badge){badge.textContent=reviewLabel(status);}"
        "card.querySelectorAll('.review-btn').forEach(btn=>{"
        "const active=btn.getAttribute('data-review-value')===status;"
        "btn.classList.toggle('is-active', active);"
        "});"
        "}"
        "function syncReviewSummary(){"
        "const cards=getSvCards();"
        "const counts={real:0,not_real:0,undecided:0};"
        "let visible=0;"
        "cards.forEach(card=>{"
        "const status=normalizeReviewStatus(card.getAttribute('data-review-status'));"
        "counts[status]=(counts[status]||0)+1;"
        "if(card.style.display!=='none'){visible+=1;}"
        "});"
        "const el=document.getElementById('review-summary');"
        "if(el){"
        "el.textContent=`real: ${counts.real} | not real: ${counts.not_real} | undecided: ${counts.undecided} | visible: ${visible}/${cards.length}`;"
        "}"
        "}"
        "function applyReviewFilter(){"
        "const sel=document.getElementById('review-filter');"
        "const wanted=sel?String(sel.value||'all'):'all';"
        "getSvCards().forEach(card=>{"
        "const status=normalizeReviewStatus(card.getAttribute('data-review-status'));"
        "const show=(wanted==='all'||wanted===status);"
        "card.style.display=show?'':'none';"
        "});"
        "syncReviewSummary();"
        "renderSummaries(wanted);"
        "}"
        "function setSvReview(btn){"
        "const card=btn.closest('.sv[data-sv-id]');"
        "if(!card){return;}"
        "const svId=card.getAttribute('data-sv-id')||'';"
        "const status=normalizeReviewStatus(btn.getAttribute('data-review-value'));"
        "reviewState[svId]=status;"
        "applyReviewToCard(card, status);"
        "saveReviewState();"
        "applyReviewFilter();"
        "}"
        "function collectReviewRows(filterValue){"
        "const rows=[];"
        "getSvCards().forEach(card=>{"
        "const status=normalizeReviewStatus(card.getAttribute('data-review-status'));"
        "if(filterValue!=='all'&&status!==filterValue){return;}"
        "rows.push({"
        "id:card.getAttribute('data-sv-id')||'',"
        "sv_len:card.getAttribute('data-sv-len')||'NA',"
        "review_status:reviewLabel(status),"
        "primary_celltype:card.getAttribute('data-primary-celltype')||'NA',"
        "linked_celltypes:card.getAttribute('data-linked-celltypes')||'NA',"
        "assigned_code:card.getAttribute('data-assigned-code')||'NA',"
        "majority_pct:card.getAttribute('data-majority-pct')||'NA',"
        "overlap_pct:card.getAttribute('data-overlap-pct')||'NA',"
        "n_supporting:card.getAttribute('data-n-supporting')||'NA',"
        "n_overlapped:card.getAttribute('data-n-overlapped')||'NA',"
        "viz_status:card.getAttribute('data-viz-status')||'NA'"
        "});"
        "});"
        "return rows;"
        "}"
        "function renderExportTable(rows){"
        "const wrap=document.getElementById('review-export-container');"
        "if(!wrap){return;}"
        "if(rows.length===0){"
        "wrap.innerHTML='<div class=\"kv\">No SVs match the selected review filter.</div>';"
        "return;"
        "}"
        "const cols=['id','sv_len','review_status','primary_celltype','linked_celltypes','assigned_code','majority_pct','overlap_pct','n_supporting','n_overlapped','viz_status'];"
        "const labels={id:'SV ID',sv_len:'sv_len',review_status:'review_status',primary_celltype:'primary_celltype',linked_celltypes:'linked_celltypes',assigned_code:'assigned_code',majority_pct:'majority_pct',overlap_pct:'overlap_pct',n_supporting:'n_supporting',n_overlapped:'n_overlapped',viz_status:'viz_status'};"
        "const chunks=['<table class=\"review-table\"><thead><tr>'];"
        "cols.forEach(col=>{chunks.push(`<th>${escapeHtml(labels[col]||col)}</th>`);});"
        "chunks.push('</tr></thead><tbody>');"
        "rows.forEach(row=>{"
        "chunks.push('<tr>');"
        "cols.forEach(col=>{chunks.push(`<td>${escapeHtml(row[col])}</td>`);});"
        "chunks.push('</tr>');"
        "});"
        "chunks.push('</tbody></table>');"
        "wrap.innerHTML=chunks.join('');"
        "}"
        "function downloadReviewRows(rows, filterValue){"
        "const cols=['id','sv_len','review_status','primary_celltype','linked_celltypes','assigned_code','majority_pct','overlap_pct','n_supporting','n_overlapped','viz_status'];"
        "const header=['SV ID','sv_len','review_status','primary_celltype','linked_celltypes','assigned_code','majority_pct','overlap_pct','n_supporting','n_overlapped','viz_status'];"
        "const lines=[header.join('\\t')];"
        "rows.forEach(row=>{"
        "lines.push(cols.map(col=>String(row[col]==null?'':row[col]).replace(/[\\t\\n\\r]+/g,' ')).join('\\t'));"
        "});"
        "const blob=new Blob([lines.join('\\n')+'\\n'],{type:'text/tab-separated-values;charset=utf-8'});"
        "const url=URL.createObjectURL(blob);"
        "const a=document.createElement('a');"
        "const stamp=new Date().toISOString().replace(/[:.]/g,'-');"
        "a.href=url;"
        "a.download=`sniffcell_sv_review_${filterValue||'all'}_${stamp}.tsv`;"
        "document.body.appendChild(a);"
        "a.click();"
        "document.body.removeChild(a);"
        "URL.revokeObjectURL(url);"
        "}"
        "function exportReviewTable(){"
        "const sel=document.getElementById('review-filter');"
        "const filterValue=sel?String(sel.value||'all'):'all';"
        "const rows=collectReviewRows(filterValue);"
        "renderExportTable(rows);"
        "downloadReviewRows(rows, filterValue);"
        "}"
        "function initReviewControls(){"
        "reviewState=loadReviewState();"
        "getSvCards().forEach(card=>{"
        "const svId=card.getAttribute('data-sv-id')||'';"
        "const status=normalizeReviewStatus(reviewState[svId]);"
        "reviewState[svId]=status;"
        "applyReviewToCard(card, status);"
        "});"
        "saveReviewState();"
        "applyReviewFilter();"
        "renderExportTable(collectReviewRows('all'));"
        "}"
        "initReviewControls();"
        "function copyVizCommand(btn){"
        "const cmd=btn.getAttribute('data-cmd')||'';"
        "const done=()=>{const old=btn.textContent;btn.textContent='Copied';setTimeout(()=>{btn.textContent=old;},1200);};"
        "if(navigator.clipboard&&window.isSecureContext){navigator.clipboard.writeText(cmd).then(done).catch(()=>{});return;}"
        "const ta=document.createElement('textarea');ta.value=cmd;document.body.appendChild(ta);ta.select();"
        "try{document.execCommand('copy');done();}catch(e){}"
        "document.body.removeChild(ta);"
        "}"
        "</script>"
    )
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

    # Fast profile intentionally trades panel detail for speed.
    if figure_profile == "fast":
        effective_max_reads = min(requested_max_reads, 120)
        effective_dpi = min(requested_dpi, 160)

    report_dir, figure_dir, html_path, archive_path = _resolve_report_paths(anno_output, args.output)
    report_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    igvviz_root = report_dir / "igvviz"
    igvviz_root.mkdir(parents=True, exist_ok=True)

    sv_df = _load_sv_assignment(sv_assignment_path)
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

    rendered_count = 0
    failed_count = 0
    igvviz_rendered_count = 0
    igvviz_failed_count = 0
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

    selected_report_df = pd.DataFrame(rows_for_report)
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
    dashboard_records = _build_dashboard_records(selected_report_df)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_text = _build_report_html(
        generated_at=generated_at,
        anno_output=anno_output,
        sv_assignment_path=sv_assignment_path,
        filters=filters,
        viz=viz_cfg,
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
        "counts": {
            "sv_total": int(len(sv_df)),
            "sv_selected": int(len(selected)),
            "viz_rendered_or_reused": int(rendered_count),
            "viz_failed": int(failed_count),
            "igvviz_rendered_or_reused": int(igvviz_rendered_count),
            "igvviz_failed": int(igvviz_failed_count),
        },
        "outputs": {
            "report_dir": str(report_dir.resolve()),
            "html": str(html_path.resolve()),
            "figures_dir": str(figure_dir.resolve()),
            "igvviz_dir": str(igvviz_root.resolve()),
            "high_confidence_sv_tsv": str(selected_tsv_path.resolve()),
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
    logger.info("Wrote failed viz table: %s", failed_tsv_path)
    logger.info("Wrote failed igvviz table: %s", failed_igvviz_tsv_path)
    logger.info("Wrote report manifest: %s", manifest_out)
    if archive_out is not None:
        logger.info("Wrote gzipped report archive: %s", archive_out)
