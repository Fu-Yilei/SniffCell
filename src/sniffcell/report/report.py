from __future__ import annotations

import html
import json
import logging
import re
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from sniffcell.viz import viz as viz_module


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


def _resolve_report_paths(anno_output: Path, output: str | None) -> tuple[Path, Path, Path]:
    if output is None:
        report_dir = anno_output / "report"
        html_path = report_dir / "index.html"
    else:
        out = Path(output)
        if out.suffix.lower() == ".html":
            html_path = out
            report_dir = out.parent if str(out.parent) else Path(".")
        else:
            report_dir = out
            html_path = report_dir / "index.html"
    figure_dir = report_dir / "figures"
    return report_dir, figure_dir, html_path


def _build_report_html(
    *,
    generated_at: str,
    anno_output: Path,
    sv_assignment_path: Path,
    filters: dict[str, object],
    total_sv: int,
    selected_count: int,
    rendered_count: int,
    failed_count: int,
    rows: list[dict[str, object]],
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
        "img{max-width:100%;height:auto;border:1px solid #d8dee4;border-radius:8px;background:#fff;}"
        ".empty{background:#ffffff;border-radius:10px;padding:18px;font-size:15px;}"
        "code{background:#edf2f7;padding:1px 4px;border-radius:4px;}"
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
        f"max_sv={filters['max_sv']}"
        "</div>"
    )

    if not rows:
        page.append("<div class=\"empty\">No SVs passed the report filters.</div>")
    else:
        for item in rows:
            sv_id = html.escape(str(item["id"]))
            linked = html.escape(str(item.get("linked_celltypes", "")))
            primary = html.escape(str(item.get("primary_celltype", "")))
            assigned_code = html.escape(str(item.get("assigned_code", "")))
            majority = _fmt_float(item.get("majority_pct"))
            overlap = _fmt_float(item.get("overlap_pct"))
            n_supporting = "NA" if pd.isna(item.get("n_supporting")) else int(item["n_supporting"])
            n_overlapped = "NA" if pd.isna(item.get("n_overlapped")) else int(item["n_overlapped"])
            status = html.escape(str(item.get("viz_status", "")))

            page.append("<section class=\"sv\">")
            page.append(f"<h2>{sv_id}</h2>")
            page.append(f"<div class=\"kv\"><b>Primary cell type:</b> {primary or 'NA'}</div>")
            page.append(f"<div class=\"kv\"><b>Linked cell types:</b> {linked or 'NA'}</div>")
            page.append(f"<div class=\"kv\"><b>Assigned code:</b> {assigned_code or 'NA'}</div>")
            page.append(f"<div class=\"kv\"><b>majority_pct:</b> {majority} | <b>overlap_pct:</b> {overlap}</div>")
            page.append(f"<div class=\"kv\"><b>n_supporting:</b> {n_supporting} | <b>n_overlapped:</b> {n_overlapped}</div>")
            page.append(f"<div class=\"kv\"><b>viz status:</b> {status}</div>")

            err = str(item.get("viz_error", "")).strip()
            if err:
                page.append(f"<div class=\"kv err\">viz error: {html.escape(err)}</div>")

            fig_rel = str(item.get("viz_figure_rel", "")).strip()
            if fig_rel:
                page.append(f"<div style=\"margin-top:10px\"><img src=\"{html.escape(fig_rel)}\" alt=\"SV plot for {sv_id}\" loading=\"lazy\"></div>")
            page.append("</section>")

    page.append("</div>")
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

    report_dir, figure_dir, html_path = _resolve_report_paths(anno_output, args.output)
    report_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

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
    if len(selected) >= 25:
        logger.warning(
            "Report will render %d SV panels; this can take a long time. "
            "Use --max_sv to limit candidates or --reuse_existing_viz to skip existing plots.",
            len(selected),
        )

    rows_for_report: list[dict[str, object]] = []
    slug_counts: dict[str, int] = {}
    rendered_count = 0
    failed_count = 0

    for row in selected.to_dict(orient="records"):
        sv_id = str(row["id"])
        base_slug = _safe_slug(sv_id)
        slug_counts[base_slug] = slug_counts.get(base_slug, 0) + 1
        slug = base_slug if slug_counts[base_slug] == 1 else f"{base_slug}_{slug_counts[base_slug]}"
        figure_path = figure_dir / f"{slug}.viz.{args.format}"

        viz_status = "pending"
        viz_error = ""

        try:
            if bool(args.reuse_existing_viz) and figure_path.exists():
                viz_status = "reused"
            else:
                viz_args = SimpleNamespace(
                    anno_output=str(anno_output),
                    sv_id=sv_id,
                    input=None,
                    vcf=None,
                    reference=None,
                    bed=None,
                    read_assignment=None,
                    kanpig_read_names=None,
                    window=int(args.window),
                    max_reads=int(args.max_reads),
                    format=args.format,
                    export_tables=bool(args.export_tables),
                    output=str(figure_path),
                )
                viz_module.viz_main(viz_args)
                viz_status = "rendered"

            if figure_path.exists():
                rendered_count += 1
            else:
                failed_count += 1
                viz_status = "failed_no_output"
                viz_error = "viz completed but output figure was not found."
        except Exception as exc:
            failed_count += 1
            viz_status = "failed"
            viz_error = str(exc)
            logger.exception("viz failed for SV %s", sv_id)

        figure_rel = ""
        if figure_path.exists():
            figure_rel = figure_path.relative_to(html_path.parent).as_posix()

        report_row = dict(row)
        report_row["viz_status"] = viz_status
        report_row["viz_error"] = viz_error
        report_row["viz_figure"] = str(figure_path)
        report_row["viz_figure_rel"] = figure_rel
        rows_for_report.append(report_row)

    selected_report_df = pd.DataFrame(rows_for_report)
    for col in ("viz_status", "viz_error", "viz_figure", "viz_figure_rel"):
        if col not in selected_report_df.columns:
            selected_report_df[col] = pd.Series(dtype="string")
    selected_tsv_path = report_dir / "high_confidence_sv.tsv"
    selected_report_df.to_csv(selected_tsv_path, sep="\t", index=False)

    failed_tsv_path = report_dir / "failed_viz.tsv"
    failed_only = selected_report_df[selected_report_df["viz_status"].astype(str).str.startswith("failed")].copy()
    failed_only.to_csv(failed_tsv_path, sep="\t", index=False)

    filters = {
        "min_overlap_pct": float(args.min_overlap_pct),
        "min_majority_pct": float(args.min_majority_pct),
        "include_unassigned": bool(args.include_unassigned),
        "allow_hard_conflict": bool(args.allow_hard_conflict),
        "max_sv": int(args.max_sv),
    }
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_text = _build_report_html(
        generated_at=generated_at,
        anno_output=anno_output,
        sv_assignment_path=sv_assignment_path,
        filters=filters,
        total_sv=int(len(sv_df)),
        selected_count=int(len(selected)),
        rendered_count=int(rendered_count),
        failed_count=int(failed_count),
        rows=rows_for_report,
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
        "viz": {
            "window": int(args.window),
            "max_reads": int(args.max_reads),
            "format": str(args.format),
            "export_tables": bool(args.export_tables),
            "reuse_existing_viz": bool(args.reuse_existing_viz),
        },
        "counts": {
            "sv_total": int(len(sv_df)),
            "sv_selected": int(len(selected)),
            "viz_rendered_or_reused": int(rendered_count),
            "viz_failed": int(failed_count),
        },
        "outputs": {
            "report_dir": str(report_dir.resolve()),
            "html": str(html_path.resolve()),
            "figures_dir": str(figure_dir.resolve()),
            "high_confidence_sv_tsv": str(selected_tsv_path.resolve()),
            "failed_viz_tsv": str(failed_tsv_path.resolve()),
        },
    }
    manifest_out = report_dir / "report_manifest.json"
    manifest_out.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True), encoding="utf-8")

    logger.info("Wrote HTML report: %s", html_path)
    logger.info("Wrote selected SV table: %s", selected_tsv_path)
    logger.info("Wrote failed viz table: %s", failed_tsv_path)
    logger.info("Wrote report manifest: %s", manifest_out)
