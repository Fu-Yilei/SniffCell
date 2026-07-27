from __future__ import annotations

import base64
from datetime import datetime
import gzip
import html
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import pandas as pd

from sniffcell.viz import igvviz as igvviz_module
from sniffcell.viz import viz as viz_module

_INFO_COLUMNS = [
    "site_locus",
    "primary_celltype",
    "linked_celltypes",
    "assigned_code",
    "majority_pct",
    "overlap_pct",
    "n_supporting",
    "n_overlapped",
]

_BASEMOD_COLOR_MODE = "basemod2"
_SUPPORT_READ_HIGHLIGHT_COLOR = "#f59e0b"
_SUPPORTED_TRACK_SUFFIXES = (
    ".bam",
    ".cram",
    ".vcf",
    ".vcf.gz",
    ".bed",
    ".bed.gz",
    ".bedgraph",
    ".bedgraph.gz",
    ".bb",
    ".bw",
    ".bigbed",
    ".bigwig",
    ".gff",
    ".gff.gz",
    ".gff3",
    ".gff3.gz",
    ".gtf",
    ".gtf.gz",
    ".bedpe",
    ".bedpe.gz",
)

_SESSION_DICTIONARY_RE = re.compile(
    r"(const sessionDictionary = )(.+?)(\n\s*let igvBrowser)",
    flags=re.DOTALL,
)


def _fmt_site_scalar(value: object) -> str:
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    text = str(value).strip()
    return text


def _resolve_igvreport_runner() -> list[str]:
    create_report_path = shutil.which("create_report")
    if create_report_path:
        return [create_report_path]

    if importlib.util.find_spec("igv_reports.report") is not None:
        return [sys.executable, "-m", "igv_reports.report"]

    raise RuntimeError(
        "igv-reports is not installed. Install it with `pip install igv-reports`."
    )


def _resolve_igvreport_inputs(anno_output: Path, igv_bams: list[str] | None) -> dict[str, object]:
    manifest = viz_module._load_anno_manifest(str(anno_output))
    manifest_inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}

    bam_paths = igvviz_module._split_bam_args(igv_bams)
    if not bam_paths:
        manifest_bam = manifest_inputs.get("bam")
        if manifest_bam:
            bam_paths = [str(manifest_bam)]

    reference_path = manifest_inputs.get("reference")
    if not reference_path:
        raise ValueError(
            "igvreport needs a reference FASTA from anno_run_manifest.json. "
            "Run sniffcell anno with a reference, or update the manifest."
        )

    vcf_path = manifest_inputs.get("vcf")
    bed_path = manifest_inputs.get("bed")
    if (not bam_paths) and (not vcf_path) and (not bed_path):
        raise ValueError(
            "igvreport needs at least one track file (BAM/VCF/BED) from --igv_bams or anno_run_manifest.json."
        )

    return {
        "bam_paths": [str(Path(x).expanduser()) for x in bam_paths],
        "reference_path": str(reference_path),
        "vcf_path": (str(vcf_path) if vcf_path else None),
        "bed_path": (str(bed_path) if bed_path else None),
    }


def _sv_interval_from_row(row: pd.Series, vcf_path: str | None) -> tuple[str, int, int]:
    sv_id = str(row.get("id", "")).strip()
    if vcf_path:
        try:
            payload = viz_module._get_sv_payload(vcf_path, sv_id)
            return str(payload["chrom"]), int(payload["start"]), int(payload["end"])
        except Exception:
            pass

    chrom = _fmt_site_scalar(row.get("sv_chr", ""))
    if not chrom:
        raise ValueError(f"Could not resolve locus for SV '{sv_id}'. Missing sv_chr and VCF lookup failed.")

    sv_pos_raw = row.get("sv_pos", pd.NA)
    if pd.isna(sv_pos_raw):
        raise ValueError(f"Could not resolve locus for SV '{sv_id}'. Missing sv_pos and VCF lookup failed.")

    start_1 = int(float(sv_pos_raw))
    start_0 = max(0, start_1 - 1)
    sv_len_raw = row.get("sv_len", pd.NA)
    if pd.isna(sv_len_raw):
        end_0 = start_0 + 1
    else:
        end_0 = start_0 + max(1, abs(int(float(sv_len_raw))))
    return chrom, start_0, end_0


def _build_sites_table(selected_df: pd.DataFrame, vcf_path: str | None) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for item in selected_df.to_dict(orient="records"):
        row = pd.Series(item)
        chrom, start_0, end_0 = _sv_interval_from_row(row, vcf_path)
        sv_id = str(item.get("id", "")).strip()
        sv_locus = f"{chrom}:{start_0 + 1}-{end_0}"
        rows.append(
            {
                "chrom": chrom,
                "start": int(start_0),
                "end": int(end_0),
                "id": sv_id,
                "site_locus": sv_locus,
                "sv_locus": sv_locus,
                "variant_class": _fmt_site_scalar(item.get("variant_class", "")),
                "primary_celltype": _fmt_site_scalar(item.get("primary_celltype", "")),
                "linked_celltypes": _fmt_site_scalar(item.get("linked_celltypes", "")),
                "assigned_code": _fmt_site_scalar(item.get("assigned_code", "")),
                "majority_pct": _fmt_site_scalar(item.get("majority_pct", "")),
                "overlap_pct": _fmt_site_scalar(item.get("overlap_pct", "")),
                "n_supporting": _fmt_site_scalar(item.get("n_supporting", "")),
                "n_overlapped": _fmt_site_scalar(item.get("n_overlapped", "")),
            }
        )
    return pd.DataFrame(rows)


def _write_header_html(
    *,
    header_path: Path,
    native_report_html: Path,
    anno_output: Path,
    selected_count: int,
    track_count: int,
) -> None:
    rel_native = os.path.relpath(str(native_report_html.resolve()), start=str(header_path.parent.resolve()))
    header_html = (
        "<section style=\"font-family:Helvetica,Arial,sans-serif;padding:16px 20px 8px 20px;"
        "background:#f3f5f7;border-bottom:1px solid #d8dee4;margin-bottom:10px;\">"
        "<div style=\"max-width:1200px;margin:0 auto;\">"
        "<div style=\"font-size:24px;font-weight:700;color:#0f2233;margin-bottom:6px;\">"
        "SniffCell Alternate IGV Report"
        "</div>"
        "<div style=\"font-size:14px;color:#44525f;line-height:1.5;\">"
        f"Selected variants: {selected_count} | Tracks: {track_count} | "
        f"anno_output: <code>{html.escape(str(anno_output))}</code> | "
        f"<a href=\"{html.escape(rel_native)}\">Open native SniffCell report</a>"
        "</div>"
        "<div style=\"font-size:13px;color:#5b6b77;line-height:1.5;margin-top:6px;\">"
        "Alignment tracks default to DNA methylation two-color mode; "
        "Variant-supporting reads are highlighted when read names are available."
        "</div>"
        "</div>"
        "</section>"
    )
    header_path.write_text(header_html, encoding="utf-8")


def _is_supported_track_path(path_text: str | None) -> bool:
    text = str(path_text or "").strip().lower()
    if not text:
        return False
    return any(text.endswith(suffix) for suffix in _SUPPORTED_TRACK_SUFFIXES)


def _decode_igv_data_uri_json(data_uri: str) -> dict[str, object]:
    text = str(data_uri).strip()
    if "," not in text:
        raise ValueError("Unsupported IGV data URI format.")
    prefix, payload = text.split(",", 1)
    raw = base64.b64decode(payload)
    if "application/gzip" in prefix:
        raw = gzip.decompress(raw)
    return json.loads(raw.decode("utf-8"))


def _encode_igv_data_uri_json(payload: dict[str, object]) -> str:
    session_text = json.dumps(payload, separators=(",", ":"))
    compressed = gzip.compress(session_text.encode("utf-8"))
    return "data:application/gzip;base64," + base64.b64encode(compressed).decode("ascii")


def _normalize_supporting_read_names(value: object) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for name in viz_module._parse_support_read_names(value):
        text = str(name).strip()
        if (not text) or (text in seen):
            continue
        seen.add(text)
        out.append(text)
    return out


def _build_supporting_reads_by_session(selected_df: pd.DataFrame, vcf_path: str | None) -> dict[str, list[str]]:
    support_map: dict[str, list[str]] = {}
    for idx, item in enumerate(selected_df.to_dict(orient="records")):
        support_names = _normalize_supporting_read_names(item.get("supporting_reads", ""))
        group_a_names = _normalize_supporting_read_names(item.get("group_a_read_names", ""))
        group_b_names = _normalize_supporting_read_names(item.get("group_b_read_names", ""))
        if group_a_names or group_b_names:
            support_names = _normalize_supporting_read_names([*support_names, *group_a_names, *group_b_names])
        sv_id = str(item.get("id", "")).strip()
        if (not support_names) and vcf_path and sv_id:
            try:
                payload = viz_module._get_sv_payload(vcf_path, sv_id)
                support_names = sorted(str(x).strip() for x in payload.get("supporting_reads", set()) if str(x).strip())
            except Exception:
                support_names = []
        support_map[str(idx)] = support_names
    return support_map


def _apply_alignment_track_defaults(session_payload: dict[str, object]) -> dict[str, object]:
    out = dict(session_payload)
    raw_tracks = out.get("tracks", [])
    if not isinstance(raw_tracks, list):
        return out

    tracks: list[object] = []
    for track in raw_tracks:
        if not isinstance(track, dict):
            tracks.append(track)
            continue
        next_track = dict(track)
        if str(next_track.get("type", "")).strip().lower() == "alignment":
            next_track["colorBy"] = _BASEMOD_COLOR_MODE
        tracks.append(next_track)
    out["tracks"] = tracks
    return out


def _inject_igvreport_custom_js(html_text: str, supporting_reads_by_session: dict[str, list[str]]) -> str:
    helper_block = (
        f"\nlet igvBrowser\n"
        "let sniffcellActiveSessionId = \"0\"\n"
        "let sniffcellSupportHighlightEnabled = true\n"
        f"const sniffcellSupportingReadsBySession = {json.dumps(supporting_reads_by_session, sort_keys=True)}\n"
        f"const sniffcellSupportingReadHighlightColor = {json.dumps(_SUPPORT_READ_HIGHLIGHT_COLOR)}\n\n"
        "function sniffcellGetSupportingReads(uniqueId) {\n"
        "    const sessionKey = String(uniqueId)\n"
        "    const readNames = sniffcellSupportingReadsBySession[sessionKey]\n"
        "    return Array.isArray(readNames) ? readNames : []\n"
        "}\n\n"
        "function sniffcellEnsureSupportToolbar() {\n"
        "    if (document.getElementById(\"sniffcell-support-toolbar\")) {\n"
        "        return\n"
        "    }\n"
        "    const host = document.getElementById(\"igvContainer\") || document.getElementById(\"container\") || document.body\n"
        "    const toolbar = document.createElement(\"div\")\n"
        "    toolbar.id = \"sniffcell-support-toolbar\"\n"
        "    toolbar.style.cssText = \"display:flex;flex-wrap:wrap;gap:10px;align-items:center;padding:10px 12px;margin:0 0 10px 0;background:#f7f8fa;border:1px solid #d8dee4;border-radius:8px;font-family:Helvetica,Arial,sans-serif;\"\n"
        "\n"
        "    const title = document.createElement(\"strong\")\n"
        "    title.textContent = \"Supporting reads\"\n"
        "    title.style.color = \"#0f2233\"\n"
        "    toolbar.appendChild(title)\n"
        "\n"
        "    const highlightBtn = document.createElement(\"button\")\n"
        "    highlightBtn.id = \"sniffcell-highlight-support-btn\"\n"
        "    highlightBtn.type = \"button\"\n"
        "    highlightBtn.textContent = \"Highlight supporting reads\"\n"
        "    highlightBtn.style.cssText = \"padding:6px 10px;border:1px solid #cbd5e1;border-radius:6px;background:#ffffff;cursor:pointer;\"\n"
        "    highlightBtn.onclick = function () {\n"
        "        sniffcellApplySupportingReadHighlights(sniffcellActiveSessionId, true)\n"
        "    }\n"
        "    toolbar.appendChild(highlightBtn)\n"
        "\n"
        "    const clearBtn = document.createElement(\"button\")\n"
        "    clearBtn.id = \"sniffcell-clear-support-btn\"\n"
        "    clearBtn.type = \"button\"\n"
        "    clearBtn.textContent = \"Clear support highlight\"\n"
        "    clearBtn.style.cssText = \"padding:6px 10px;border:1px solid #cbd5e1;border-radius:6px;background:#ffffff;cursor:pointer;\"\n"
        "    clearBtn.onclick = function () {\n"
        "        sniffcellApplySupportingReadHighlights(sniffcellActiveSessionId, false)\n"
        "    }\n"
        "    toolbar.appendChild(clearBtn)\n"
        "\n"
        "    const status = document.createElement(\"span\")\n"
        "    status.id = \"sniffcell-support-status\"\n"
        "    status.style.cssText = \"color:#44525f;font-size:13px;\"\n"
        "    toolbar.appendChild(status)\n"
        "\n"
        "    if (host.firstChild) {\n"
        "        host.insertBefore(toolbar, host.firstChild)\n"
        "    } else {\n"
        "        host.appendChild(toolbar)\n"
        "    }\n"
        "}\n\n"
        "function sniffcellUpdateSupportToolbar(uniqueId) {\n"
        "    const readNames = sniffcellGetSupportingReads(uniqueId)\n"
        "    const highlightBtn = document.getElementById(\"sniffcell-highlight-support-btn\")\n"
        "    const clearBtn = document.getElementById(\"sniffcell-clear-support-btn\")\n"
        "    const status = document.getElementById(\"sniffcell-support-status\")\n"
        "    if (highlightBtn) {\n"
        "        highlightBtn.disabled = readNames.length === 0\n"
        "    }\n"
        "    if (clearBtn) {\n"
        "        clearBtn.disabled = readNames.length === 0\n"
        "    }\n"
        "    if (status) {\n"
        "        if (readNames.length === 0) {\n"
        "            status.textContent = \"No supporting read names available for this SV.\"\n"
        "        } else if (sniffcellSupportHighlightEnabled) {\n"
        "            status.textContent = \"Supporting reads: \" + String(readNames.length) + \". Highlighting is supported; true read selection/filtering is not exposed by IGV.js here.\"\n"
        "        } else {\n"
        "            status.textContent = \"Supporting reads: \" + String(readNames.length) + \". Highlight is currently off.\"\n"
        "        }\n"
        "    }\n"
        "}\n\n"
        "function sniffcellFindAlignmentTracks() {\n"
        "    if (!igvBrowser) {\n"
        "        return []\n"
        "    }\n"
        "    try {\n"
        "        if (typeof igvBrowser.findTracks === \"function\") {\n"
        "            const tracks = igvBrowser.findTracks(\"type\", \"alignment\")\n"
        "            if (Array.isArray(tracks)) {\n"
        "                return tracks\n"
        "            }\n"
        "        }\n"
        "    } catch (err) {\n"
        "    }\n"
        "    try {\n"
        "        if (Array.isArray(igvBrowser.trackViews)) {\n"
        "            return igvBrowser.trackViews\n"
        "                .map(function (trackView) { return trackView ? trackView.track : null })\n"
        "                .filter(function (track) { return track && track.type === \"alignment\" })\n"
        "        }\n"
        "    } catch (err) {\n"
        "    }\n"
        "    return []\n"
        "}\n\n"
        "function sniffcellApplySupportingReadHighlights(uniqueId, enabled) {\n"
        "    const sessionKey = String(uniqueId)\n"
        "    sniffcellActiveSessionId = sessionKey\n"
        "    if (typeof enabled === \"boolean\") {\n"
        "        sniffcellSupportHighlightEnabled = enabled\n"
        "    }\n"
        "    const readNames = sniffcellSupportHighlightEnabled ? sniffcellGetSupportingReads(sessionKey) : []\n"
        "    let attempts = 0\n\n"
        "    function applyOnce() {\n"
        "        const tracks = sniffcellFindAlignmentTracks()\n"
        "        if (!tracks.length) {\n"
        "            return false\n"
        "        }\n"
        "        tracks.forEach(function (track) {\n"
        "            if (track && typeof track.setHighlightedReads === \"function\") {\n"
        "                track.setHighlightedReads(readNames, sniffcellSupportingReadHighlightColor)\n"
        "            }\n"
        "        })\n"
        "        try {\n"
        "            if (igvBrowser && typeof igvBrowser.repaintViews === \"function\") {\n"
        "                igvBrowser.repaintViews()\n"
        "            }\n"
        "        } catch (err) {\n"
        "        }\n"
        "        sniffcellUpdateSupportToolbar(sessionKey)\n"
        "        return true\n"
        "    }\n\n"
        "    function retry() {\n"
        "        attempts += 1\n"
        "        if (applyOnce() || attempts >= 25) {\n"
        "            sniffcellUpdateSupportToolbar(sessionKey)\n"
        "            return\n"
        "        }\n"
        "        window.setTimeout(retry, 200)\n"
        "    }\n\n"
        "    retry()\n"
        "}\n"
    )
    if "let igvBrowser" not in html_text:
        raise ValueError("Could not find igvBrowser declaration in igvreport HTML.")
    out = html_text.replace("let igvBrowser", helper_block, 1)

    create_browser_old = "igvBrowser = b\n                initTable()"
    create_browser_new = "igvBrowser = b\n                sniffcellEnsureSupportToolbar()\n                sniffcellApplySupportingReadHighlights(\"0\")\n                initTable()"
    if create_browser_old not in out:
        raise ValueError("Could not find IGV browser initialization block in igvreport HTML.")
    out = out.replace(create_browser_old, create_browser_new, 1)

    load_session_old = (
        "                igvBrowser.loadSession({\n"
        "                    url: session\n"
        "                })"
    )
    load_session_new = (
        "                igvBrowser.loadSession({\n"
        "                    url: session\n"
        "                }).then(function () {\n"
        "                    sniffcellApplySupportingReadHighlights(uniqueId, sniffcellSupportHighlightEnabled)\n"
        "                }).catch(function () {\n"
        "                    sniffcellApplySupportingReadHighlights(uniqueId, sniffcellSupportHighlightEnabled)\n"
        "                })"
    )
    if load_session_old not in out:
        raise ValueError("Could not find session load block in igvreport HTML.")
    return out.replace(load_session_old, load_session_new, 1)


def _customize_igvreport_html(
    html_path: Path,
    *,
    supporting_reads_by_session: dict[str, list[str]],
) -> None:
    html_text = html_path.read_text(encoding="utf-8")
    match = _SESSION_DICTIONARY_RE.search(html_text)
    if match is None:
        raise ValueError("Could not locate sessionDictionary in igvreport HTML.")

    session_dict = json.loads(match.group(2))
    if not isinstance(session_dict, dict):
        raise ValueError("igvreport sessionDictionary is not a JSON object.")

    customized_sessions: dict[str, str] = {}
    for session_id, session_uri in session_dict.items():
        if not isinstance(session_id, str):
            session_id = str(session_id)
        if not isinstance(session_uri, str):
            raise ValueError(f"Unexpected session URI payload for session {session_id!r}.")
        session_payload = _decode_igv_data_uri_json(session_uri)
        customized_sessions[session_id] = _encode_igv_data_uri_json(
            _apply_alignment_track_defaults(session_payload)
        )

    out = html_text[: match.start(2)] + json.dumps(customized_sessions, sort_keys=True) + html_text[match.end(2) :]
    out = _inject_igvreport_custom_js(out, supporting_reads_by_session)
    html_path.write_text(out, encoding="utf-8")


def build_igvreport_cli_command(
    *,
    runner_prefix: list[str],
    sites_path: Path,
    output_html: Path,
    reference_path: str,
    tracks: list[str],
    info_columns: list[str],
    flanking: int,
    header_path: Path,
) -> list[str]:
    cmd = [
        *runner_prefix,
        str(sites_path),
        "--fasta",
        str(reference_path),
        "--output",
        str(output_html),
        "--sequence",
        "1",
        "--begin",
        "2",
        "--end",
        "3",
        "--zero_based",
        "true",
        "--flanking",
        str(int(flanking)),
        "--title",
        "SniffCell Alternate IGV Report",
        "--header",
        str(header_path),
        "--info-columns",
        *[str(x) for x in info_columns],
    ]
    if tracks:
        cmd.extend(["--tracks", *[str(x) for x in tracks]])
    return cmd


def render_igvreport_bundle(
    *,
    anno_output: Path,
    selected_df: pd.DataFrame,
    output_dir: Path,
    native_report_html: Path,
    igv_bams: list[str] | None,
    window: int,
) -> dict[str, object]:
    if selected_df.empty:
        raise ValueError("igvreport needs at least one selected variant.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    html_path = output_dir / "index.html"
    manifest_path = output_dir / "igvreport_manifest.json"
    sites_path = output_dir / "selected_sites.tsv"
    legacy_sites_path = output_dir / "selected_sv_sites.tsv"
    header_path = output_dir / "sniffcell_igvreport_header.html"

    base_result: dict[str, object] = {
        "status": "not_rendered",
        "error": "",
        "command": "",
        "html_path": str(html_path.resolve()),
        "manifest_path": str(manifest_path.resolve()),
        "sites_path": str(sites_path.resolve()),
        "legacy_sites_path": str(legacy_sites_path.resolve()),
        "header_path": str(header_path.resolve()),
    }

    stdout_text = ""
    stderr_text = ""
    status = "failed"
    error = ""
    command_text = ""
    inputs_payload: dict[str, object] = {}
    tracks: list[str] = []
    postprocess_status = "not_applied"
    postprocess_error = ""
    supporting_reads_available_sessions = 0

    try:
        resolved = _resolve_igvreport_inputs(anno_output, igv_bams)
        bam_paths = list(resolved["bam_paths"])
        reference_path = str(resolved["reference_path"])
        vcf_path = resolved["vcf_path"]
        bed_path = resolved["bed_path"]
        supporting_reads_by_session = _build_supporting_reads_by_session(selected_df, vcf_path)
        supporting_reads_available_sessions = int(sum(1 for reads in supporting_reads_by_session.values() if reads))

        tracks = [*bam_paths]
        if vcf_path:
            tracks.append(str(vcf_path))
        if bed_path and _is_supported_track_path(bed_path):
            tracks.append(str(bed_path))

        _write_header_html(
            header_path=header_path,
            native_report_html=native_report_html,
            anno_output=anno_output,
            selected_count=int(len(selected_df)),
            track_count=int(len(tracks)),
        )

        sites_df = _build_sites_table(selected_df, vcf_path)

        runner_prefix = _resolve_igvreport_runner()
        command = build_igvreport_cli_command(
            runner_prefix=runner_prefix,
            sites_path=sites_path,
            output_html=html_path,
            reference_path=reference_path,
            tracks=tracks,
            info_columns=_INFO_COLUMNS,
            flanking=int(window),
            header_path=header_path,
        )
        command_text = " ".join(shlex.quote(part) for part in command)

        sites_df.to_csv(sites_path, sep="\t", index=False)
        sites_df.to_csv(legacy_sites_path, sep="\t", index=False)
        proc = subprocess.run(command, check=False, capture_output=True, text=True)
        stdout_text = str(proc.stdout or "")
        stderr_text = str(proc.stderr or "")
        if proc.returncode != 0:
            status = "failed"
            error = (
                f"igvreport command failed (exit={proc.returncode}). "
                f"Install igv-reports. stderr: {stderr_text.strip() or 'none'}"
            )
        elif html_path.exists():
            status = "rendered"
            error = ""
            try:
                _customize_igvreport_html(
                    html_path,
                    supporting_reads_by_session=supporting_reads_by_session,
                )
                postprocess_status = "applied"
            except Exception as exc:
                postprocess_status = "failed"
                postprocess_error = str(exc)
        else:
            status = "failed_no_output"
            error = "igvreport completed but output HTML was not found."

        inputs_payload = {
            "anno_output": str(anno_output.resolve()),
            "reference": str(Path(reference_path).expanduser().resolve()),
            "bams": [str(Path(x).expanduser().resolve()) for x in bam_paths],
            "vcf": (str(Path(vcf_path).expanduser().resolve()) if vcf_path else ""),
            "bed": (str(Path(bed_path).expanduser().resolve()) if bed_path else ""),
        }
    except FileNotFoundError as exc:
        status = "failed"
        error = (
            f"igvreport runner was not found: {exc}. "
            "Install igv-reports (`pip install igv-reports`)."
        )
    except Exception as exc:
        status = "failed"
        error = str(exc)

    manifest_payload = {
        "command": "igvreport",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "status": status,
        "error": error,
        "igvreport_command": command_text,
        "inputs": inputs_payload,
        "outputs": {
            "output_dir": str(output_dir.resolve()),
            "html": str(html_path.resolve()),
            "manifest": str(manifest_path.resolve()),
            "sites_tsv": str(sites_path.resolve()),
            "selected_sites_tsv": str(sites_path.resolve()),
            "selected_sv_sites_tsv": str(legacy_sites_path.resolve()),
            "header_html": str(header_path.resolve()),
        },
        "runtime": {
            "stdout": stdout_text,
            "stderr": stderr_text,
        },
        "tracks": tracks,
        "info_columns": list(_INFO_COLUMNS),
        "flanking": int(window),
        "session_customization": {
            "alignment_color_by": _BASEMOD_COLOR_MODE,
            "supporting_read_highlight_color": _SUPPORT_READ_HIGHLIGHT_COLOR,
            "supporting_reads_available_sessions": supporting_reads_available_sessions,
            "status": postprocess_status,
            "error": postprocess_error,
        },
        "native_report_html": str(native_report_html.resolve()),
        "n_selected": int(len(selected_df)),
    }
    manifest_path.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True), encoding="utf-8")

    return {
        **base_result,
        "status": status,
        "error": error,
        "command": command_text,
    }
