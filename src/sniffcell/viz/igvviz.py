from __future__ import annotations

import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pysam

from sniffcell.viz.viz import (
    _get_sv_payload,
    _load_anno_manifest,
    _load_kanpig_supporting_reads,
    _norm_chr,
)

_BASEMOD_BINARY_THRESHOLD = 0.70
_BASEMOD_BINARY_ON_COLOR = "220,38,38"
_BASEMOD_BINARY_OFF_COLOR = "65,105,225"
_HG38_GENE_TRACK = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz"
_BASEMOD_MOD_COLOR_KEYS = (
    "BASEMOD.M_COLOR",
    "BASEMOD.H_COLOR",
    "BASEMOD.F_COLOR",
    "BASEMOD.C_COLOR",
    "BASEMOD.G_COLOR",
    "BASEMOD.E_COLOR",
    "BASEMOD.B_COLOR",
    "BASEMOD.A_COLOR",
    "BASEMOD.O_COLOR",
    "BASEMOD.17082_COLOR",
    "BASEMOD.17596_COLOR",
    "BASEMOD.21839_COLOR",
    "BASEMOD.OTHER_COLOR",
)
_BASEMOD_NONE_COLOR_KEYS = (
    "BASEMOD.NONE_A_COLOR",
    "BASEMOD.NONE_C_COLOR",
    "BASEMOD.NONE_T_COLOR",
    "BASEMOD.NONE_G_COLOR",
    "BASEMOD.NONE_N_COLOR",
)


def _split_bam_args(values: list[str] | None) -> list[str]:
    if not values:
        return []
    out: list[str] = []
    seen: set[str] = set()
    for raw in values:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        chunks = [text]
        if "," in text:
            chunks = [x.strip() for x in text.split(",") if x.strip()]
        for chunk in chunks:
            key = str(Path(chunk).expanduser())
            if key in seen:
                continue
            seen.add(key)
            out.append(key)
    return out


def _sanitize_token(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text).strip())
    cleaned = cleaned.strip("._-")
    return cleaned or "item"


def _igv_quote(value: str | Path) -> str:
    text = str(value)
    if any(ch.isspace() for ch in text):
        return f"\"{text}\""
    return text


def _resolve_igvviz_runtime_inputs(args, logger: logging.Logger) -> dict:
    anno_output = getattr(args, "anno_output", None)
    manifest = {}
    if anno_output:
        manifest = _load_anno_manifest(anno_output)
        logger.debug("Loaded anno run manifest from: %s", Path(anno_output) / "anno_run_manifest.json")

    manifest_inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}
    manifest_runtime = manifest.get("runtime", {}) if isinstance(manifest, dict) else {}

    bam_paths = _split_bam_args(getattr(args, "input", None))
    if not bam_paths:
        manifest_bam = manifest_inputs.get("bam")
        if manifest_bam:
            bam_paths = [str(manifest_bam)]

    vcf_path = args.vcf or manifest_inputs.get("vcf")
    reference_path = args.reference or manifest_inputs.get("reference")
    bed_path = args.bed or manifest_inputs.get("bed")

    if not bam_paths:
        raise ValueError(
            "igvviz needs one or more BAMs. Provide -i/--input, or set --anno_output with anno_run_manifest.json."
        )
    if not vcf_path:
        raise ValueError(
            "igvviz needs a VCF. Provide -v/--vcf, or set --anno_output with anno_run_manifest.json."
        )

    output_dir_arg = getattr(args, "output", None)
    if output_dir_arg:
        output_dir = Path(output_dir_arg)
    elif anno_output:
        output_dir = Path(anno_output) / "igvviz"
    else:
        output_dir = Path.cwd() / "igvviz"

    effective_window = int(args.window)
    if anno_output and int(args.window) == 5000 and "window" in manifest_runtime:
        try:
            effective_window = int(manifest_runtime["window"])
        except (TypeError, ValueError):
            pass

    return {
        "bam_paths": [str(Path(x).expanduser()) for x in bam_paths],
        "vcf_path": str(vcf_path),
        "reference_path": (str(reference_path) if reference_path else None),
        "bed_path": (str(bed_path) if bed_path else None),
        "window": effective_window,
        "output_dir": output_dir,
        "kanpig_read_names": getattr(args, "kanpig_read_names", None),
    }


def _candidate_chrom_names(chrom: str) -> list[str]:
    text = str(chrom).strip()
    if not text:
        return []
    norm = _norm_chr(text)
    out = [text]
    if norm:
        out.append(norm)
        out.append(f"chr{norm}")
    dedup = []
    seen = set()
    for item in out:
        if item and item not in seen:
            seen.add(item)
            dedup.append(item)
    return dedup


def _resolve_chrom_for_bam(bam_path: str, chrom: str) -> str:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = set(bam.references)
    for cand in _candidate_chrom_names(chrom):
        if cand in refs:
            return cand
    raise ValueError(f"Chromosome '{chrom}' (or chr-normalized aliases) not found in BAM header: {bam_path}")


def _load_ctdmr_markers(
    bed_path: str | None,
    sv_chrom: str,
    region_start: int,
    region_end: int,
    logger: logging.Logger,
) -> pd.DataFrame:
    cols = ["start", "end", "label", "best_dir"]
    if bed_path is None:
        return pd.DataFrame(columns=cols)

    path = Path(bed_path)
    if not path.exists():
        logger.warning("ctDMR file was not found and will be skipped: %s", path)
        return pd.DataFrame(columns=cols)

    dmrs = pd.DataFrame()
    try:
        dmrs = pd.read_csv(path, sep="\t")
        if (not dmrs.empty) and str(dmrs.columns[0]).startswith("#"):
            dmrs = dmrs.rename(columns={dmrs.columns[0]: str(dmrs.columns[0]).lstrip("#")})
        required = {"chr", "start", "end"}
        if not required.issubset(set(dmrs.columns)):
            dmrs = pd.DataFrame()
    except Exception:
        dmrs = pd.DataFrame()

    # Fallback for BED-like files without headers (for example *.igv.bed).
    if dmrs.empty:
        try:
            raw = pd.read_csv(path, sep="\t", header=None, comment="#")
            if raw.shape[1] >= 3:
                dmrs = pd.DataFrame(
                    {
                        "chr": raw.iloc[:, 0].astype(str),
                        "start": raw.iloc[:, 1],
                        "end": raw.iloc[:, 2],
                        "name": (raw.iloc[:, 3].astype(str) if raw.shape[1] >= 4 else "ctDMR"),
                    }
                )
        except Exception:
            logger.warning("Failed to parse ctDMR file; skipping marker track: %s", path)
            return pd.DataFrame(columns=cols)

    if dmrs.empty:
        return pd.DataFrame(columns=cols)

    dmrs = dmrs.copy()
    dmrs["chr_norm"] = dmrs["chr"].map(_norm_chr)
    dmrs["start"] = pd.to_numeric(dmrs["start"], errors="coerce")
    dmrs["end"] = pd.to_numeric(dmrs["end"], errors="coerce")
    dmrs = dmrs.dropna(subset=["chr_norm", "start", "end"])
    if dmrs.empty:
        return pd.DataFrame(columns=cols)
    dmrs["start"] = dmrs["start"].astype(int)
    dmrs["end"] = dmrs["end"].astype(int)
    dmrs = dmrs[dmrs["end"] > dmrs["start"]]
    dmrs = dmrs[
        (dmrs["chr_norm"] == _norm_chr(sv_chrom))
        & (dmrs["start"] < int(region_end))
        & (dmrs["end"] > int(region_start))
    ].copy()
    if dmrs.empty:
        return pd.DataFrame(columns=cols)

    label = pd.Series(["ctDMR"] * len(dmrs), index=dmrs.index, dtype="object")
    if "best_group_leaves" in dmrs.columns:
        best_group_leaves = dmrs["best_group_leaves"].astype(str).str.strip()
        label = best_group_leaves.where(best_group_leaves.ne(""), label)
    if "best_group" in dmrs.columns:
        best_group = dmrs["best_group"].astype(str).str.strip()
        label = best_group.where(best_group.ne(""), label)
    if "name" in dmrs.columns:
        name_col = dmrs["name"].astype(str).str.strip()
        label = name_col.where(name_col.ne(""), label)

    if "best_dir" not in dmrs.columns:
        dmrs["best_dir"] = ""
    dmrs["label"] = label
    return dmrs[["start", "end", "label", "best_dir"]].reset_index(drop=True)


def _write_ctdmr_track_bed(
    markers: pd.DataFrame,
    chrom_name: str,
    output_bed: Path,
) -> None:
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    if markers.empty:
        output_bed.write_text("", encoding="utf-8")
        return

    rows: list[tuple] = []
    for row in markers.itertuples(index=False):
        best_dir = str(row.best_dir).strip().lower()
        if best_dir == "hyper":
            rgb = "215,48,39"
        elif best_dir == "hypo":
            rgb = "69,117,180"
        else:
            rgb = "153,153,153"
        rows.append(
            (
                chrom_name,
                int(row.start),
                int(row.end),
                str(row.label) if str(row.label).strip() else "ctDMR",
                0,
                ".",
                int(row.start),
                int(row.end),
                rgb,
            )
        )

    bed_df = pd.DataFrame(
        rows,
        columns=["chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"],
    )
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)


def _tag_region_bam(
    bam_path: str,
    chrom_name: str,
    region_start: int,
    region_end: int,
    supporting_reads: set[str],
    support_tag: str,
    phase_tag: str,
    support_phase_group_tag: str,
    output_bam: Path,
) -> dict[str, int]:
    def _phase_bucket(read_obj) -> str:
        try:
            if read_obj.has_tag(str(phase_tag)):
                phase_value = str(read_obj.get_tag(str(phase_tag))).strip()
                if phase_value:
                    return _sanitize_token(phase_value)
        except Exception:
            pass
        return "NA"

    output_bam.parent.mkdir(parents=True, exist_ok=True)
    n_reads = 0
    n_support = 0
    with pysam.AlignmentFile(bam_path, "rb") as in_bam:
        with pysam.AlignmentFile(str(output_bam), "wb", template=in_bam) as out_bam:
            for read in in_bam.fetch(chrom_name, int(region_start), int(region_end)):
                if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                    continue
                if int(read.reference_end) <= int(region_start) or int(read.reference_start) >= int(region_end):
                    continue
                is_support = str(read.query_name) in supporting_reads
                n_reads += 1
                if is_support:
                    n_support += 1
                support_value = "A_SUPPORT" if is_support else "Z_OTHER"
                phase_value = _phase_bucket(read)
                support_phase_value = f"{support_value}_HP{phase_value}"
                read.set_tag(str(support_tag), support_value, value_type="Z")
                read.set_tag(str(support_phase_group_tag), support_phase_value, value_type="Z")
                out_bam.write(read)
    pysam.index(str(output_bam))
    return {
        "n_reads": int(n_reads),
        "n_supporting_reads": int(n_support),
        "n_non_supporting_reads": int(max(0, n_reads - n_support)),
    }


def _igv_locus(chrom: str, start0: int, end0: int) -> str:
    return f"{chrom}:{int(start0) + 1}-{int(end0)}"


def _infer_gene_track(reference_path: str | None) -> str | None:
    if not reference_path:
        return None
    ref_text = str(reference_path).strip().lower()
    if not ref_text:
        return None
    # Auto-enable compact gene labels for hg38 references.
    if any(token in ref_text for token in ("hg38", "grch38", "gca_000001405.15", "gcf_000001405.26")):
        return _HG38_GENE_TRACK
    return None


def _build_igv_batch_lines(
    *,
    jobs: list[dict],
    snapshot_dir: Path,
    reference_path: str | None,
    gene_track_path: str | None,
    visibility_window: int,
    phase_tag: str,
    support_phase_group_tag: str,
    snapshot_width: int,
    snapshot_height: int,
    hide_methylation: bool = False,
) -> list[str]:
    def _append_binary_basemod_preferences(out_lines: list[str]) -> None:
        out_lines.append(f"preference BASEMOD.THRESHOLD {_BASEMOD_BINARY_THRESHOLD}")
        out_lines.append("preference SAM.SHOW_GROUP_SEPARATOR true")
        for key in _BASEMOD_MOD_COLOR_KEYS:
            out_lines.append(f"preference {key} {_BASEMOD_BINARY_ON_COLOR}")
        for key in _BASEMOD_NONE_COLOR_KEYS:
            out_lines.append(f"preference {key} {_BASEMOD_BINARY_OFF_COLOR}")

    lines: list[str] = ["new"]
    if reference_path:
        lines.append(f"genome {_igv_quote(reference_path)}")
    lines.append(f"snapshotDirectory {_igv_quote(snapshot_dir.resolve())}")
    lines.append("maxPanelHeight 2000")
    lines.append("setSleepInterval 500")
    if int(visibility_window) > 0:
        lines.append(f"preference SAM.MAX_VISIBLE_RANGE {int(visibility_window)}")
    if not hide_methylation:
        _append_binary_basemod_preferences(lines)

    for job in jobs:
        lines.append("new")
        if reference_path:
            lines.append(f"genome {_igv_quote(reference_path)}")
        lines.append(f"snapshotDirectory {_igv_quote(snapshot_dir.resolve())}")
        lines.append("maxPanelHeight 2000")
        lines.append("setSleepInterval 500")
        if int(visibility_window) > 0:
            lines.append(f"preference SAM.MAX_VISIBLE_RANGE {int(visibility_window)}")
        if not hide_methylation:
            _append_binary_basemod_preferences(lines)
        if gene_track_path:
            lines.append(f"load {_igv_quote(gene_track_path)}")
        lines.append(f"load {_igv_quote(job['tagged_bam'])}")
        ctdmr_track = str(job.get("ctdmr_track", "")).strip()
        if ctdmr_track:
            lines.append(f"load {_igv_quote(ctdmr_track)}")
        lines.append(f"goto {job['locus']}")
        lines.append("expand")
        lines.append(f"group TAG {support_phase_group_tag}")
        lines.append(f"sort TAG {phase_tag}")
        if not hide_methylation:
            # Two-color base-mod rendering: methylated marks are colored, non-methylated marks are de-emphasized.
            lines.append("colorBy BASE_MODIFICATION_2COLOR")
        lines.append(f"snapshot {job['snapshot_name']}")

    lines.append("exit")
    return lines


def _prepare_igv_subprocess_env(work_dir: Path) -> dict[str, str]:
    env_root = work_dir / "igv_runtime_env"
    home_dir = env_root / "home"
    igv_user_dir = home_dir / "igv"
    xdg_runtime_dir = env_root / "xdg_runtime"
    igv_user_dir.mkdir(parents=True, exist_ok=True)
    xdg_runtime_dir.mkdir(parents=True, exist_ok=True)
    prefs_path = igv_user_dir / "prefs.properties"
    prefs_path.write_text(
        "PORT_ENABLED=false\n"
        "PORT_NUMBER=0\n"
        "CIRC_VIEW_ENABLED=false\n",
        encoding="utf-8",
    )
    env = os.environ.copy()
    env["HOME"] = str(home_dir)
    env["XDG_RUNTIME_DIR"] = str(xdg_runtime_dir)
    env.setdefault("DBUS_SESSION_BUS_ADDRESS", "")
    return env


def _run_igv_batch(
    igv_cmd: str,
    batch_path: Path,
    *,
    snapshot_width: int,
    snapshot_height: int,
    env: dict[str, str] | None = None,
) -> tuple[str, str]:
    cmd = shlex.split(str(igv_cmd).strip())
    if not cmd:
        raise ValueError("--igv_cmd cannot be empty")
    # IGV requires a display; auto-wrap in xvfb-run on headless nodes when possible.
    if (not os.environ.get("DISPLAY")) and (shutil.which("xvfb-run") is not None):
        first = str(cmd[0]).strip().lower()
        if ("xvfb-run" not in first) and (Path(first).name != "xvfb-run"):
            screen_w = max(1920, int(snapshot_width) + 200)
            screen_h = max(1080, int(snapshot_height) + 200)
            cmd = ["xvfb-run", "-a", "-s", f"-screen 0 {screen_w}x{screen_h}x24", *cmd]
    cmd += ["-b", str(batch_path.resolve())]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        stderr_text = str(proc.stderr or "")
        harmless_xvfb_teardown = (
            "/usr/bin/xvfb-run: line 186: kill:" in stderr_text
            and "No such process" in stderr_text
        )
        if harmless_xvfb_teardown:
            return proc.stdout, proc.stderr
        raise RuntimeError(
            f"IGV batch execution failed (exit={proc.returncode}). "
            f"cmd={' '.join(cmd)}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc.stdout, proc.stderr


def igvviz_main(args) -> None:
    logger = logging.getLogger("sniffcell.igvviz")
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    resolved = _resolve_igvviz_runtime_inputs(args, logger)
    window = int(resolved["window"])
    if window < 0:
        raise ValueError("window must be >= 0")
    raw_visibility_window = getattr(args, "visibility_window", None)
    if raw_visibility_window is None:
        visibility_window = int(window)
    else:
        visibility_window = int(raw_visibility_window)
    if visibility_window < 0:
        raise ValueError("visibility_window must be >= 0")
    snapshot_width = int(getattr(args, "snapshot_width", 3600))
    snapshot_height = int(getattr(args, "snapshot_height", 1600))
    if snapshot_width <= 0:
        raise ValueError("snapshot_width must be > 0")
    if snapshot_height <= 0:
        raise ValueError("snapshot_height must be > 0")
    hide_methylation = bool(getattr(args, "hide_methylation", False))
    support_phase_group_tag = "SG"
    if str(args.phase_tag).strip().upper() == support_phase_group_tag:
        support_phase_group_tag = "SX"
    if str(args.support_tag).strip().upper() == support_phase_group_tag:
        support_phase_group_tag = "SY"
    gene_track_path = _infer_gene_track(resolved["reference_path"])

    sv = _get_sv_payload(resolved["vcf_path"], args.sv_id)
    override_support = _load_kanpig_supporting_reads(resolved["kanpig_read_names"], args.sv_id)
    if override_support:
        sv["supporting_reads"] = override_support
    supporting_reads = set(str(x) for x in sv.get("supporting_reads", set()) if str(x).strip())

    region_start = max(0, int(sv["start"]) - window)
    region_end = int(sv["end"]) + window

    output_dir = Path(resolved["output_dir"]).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    keep_intermediates = bool(getattr(args, "keep_intermediates", False) or getattr(args, "batch_only", False))
    sv_slug = _sanitize_token(str(args.sv_id))

    temp_ctx: tempfile.TemporaryDirectory[str] | None = None
    if keep_intermediates:
        work_dir = output_dir / f"{sv_slug}.igvviz_work"
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        temp_ctx = tempfile.TemporaryDirectory(prefix=f"{sv_slug}.igvviz.", dir=str(output_dir))
        work_dir = Path(temp_ctx.name)

    markers = _load_ctdmr_markers(
        resolved["bed_path"],
        sv_chrom=str(sv["chrom"]),
        region_start=region_start,
        region_end=region_end,
        logger=logger,
    )

    jobs: list[dict] = []
    summary_rows: list[dict] = []
    try:
        igv_subprocess_env = _prepare_igv_subprocess_env(work_dir)
        for idx, bam in enumerate(resolved["bam_paths"], start=1):
            bam_path = str(Path(bam).expanduser().resolve())
            bam_slug = _sanitize_token(Path(bam_path).stem)
            chrom_for_bam = _resolve_chrom_for_bam(bam_path, str(sv["chrom"]))
            locus = _igv_locus(chrom_for_bam, region_start, region_end)

            tagged_bam = work_dir / f"{idx:02d}.{bam_slug}.{sv_slug}.tagged.region.bam"
            counts = _tag_region_bam(
                bam_path=bam_path,
                chrom_name=chrom_for_bam,
                region_start=region_start,
                region_end=region_end,
                supporting_reads=supporting_reads,
                support_tag=str(args.support_tag),
                phase_tag=str(args.phase_tag),
                support_phase_group_tag=str(support_phase_group_tag),
                output_bam=tagged_bam,
            )

            ctdmr_track_path: Path | None = None
            if not markers.empty:
                ctdmr_track_path = work_dir / f"{idx:02d}.{bam_slug}.{sv_slug}.ctdmr.bed"
                _write_ctdmr_track_bed(markers, chrom_for_bam, ctdmr_track_path)

            snapshot_name = f"{sv_slug}.{idx:02d}.{bam_slug}.igv.{args.snapshot_format}"
            jobs.append(
                {
                    "bam": bam_path,
                    "tagged_bam": str(tagged_bam.resolve()),
                    "ctdmr_track": (str(ctdmr_track_path.resolve()) if ctdmr_track_path is not None else ""),
                    "locus": locus,
                    "snapshot_name": snapshot_name,
                }
            )
            summary_rows.append(
                {
                    "bam": bam_path,
                    "tagged_bam": str(tagged_bam.resolve()),
                    "ctdmr_track": (str(ctdmr_track_path.resolve()) if ctdmr_track_path is not None else ""),
                    "locus": locus,
                    "n_reads_in_window": int(counts["n_reads"]),
                    "n_supporting_reads_in_window": int(counts["n_supporting_reads"]),
                    "n_non_supporting_reads_in_window": int(counts["n_non_supporting_reads"]),
                    "snapshot": str((output_dir / snapshot_name).resolve()),
                }
            )

        batch_lines = _build_igv_batch_lines(
            jobs=jobs,
            snapshot_dir=output_dir,
            reference_path=resolved["reference_path"],
            gene_track_path=gene_track_path,
            visibility_window=visibility_window,
            phase_tag=str(args.phase_tag),
            support_phase_group_tag=str(support_phase_group_tag),
            snapshot_width=snapshot_width,
            snapshot_height=snapshot_height,
            hide_methylation=hide_methylation,
        )
        batch_path = output_dir / f"{sv_slug}.igvviz.batch.txt"
        batch_path.write_text("\n".join(batch_lines) + "\n", encoding="utf-8")

        summary_path = output_dir / f"{sv_slug}.igvviz.summary.tsv"
        pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

        manifest_path = output_dir / f"{sv_slug}.igvviz.manifest.json"
        manifest_payload = {
            "sv_id": str(args.sv_id),
            "sv_locus": _igv_locus(str(sv["chrom"]), int(sv["start"]), int(sv["end"])),
            "window_locus": _igv_locus(str(sv["chrom"]), region_start, region_end),
            "window": int(window),
            "snapshot_width": int(snapshot_width),
            "snapshot_height": int(snapshot_height),
            "include_non_supporting": True,
            "n_supporting_reads_listed": int(len(supporting_reads)),
            "support_phase_group_tag": str(support_phase_group_tag),
            "hide_methylation": hide_methylation,
            "basemod_scheme": (None if hide_methylation else "binary_methylated_vs_none"),
            "basemod_threshold": (None if hide_methylation else float(_BASEMOD_BINARY_THRESHOLD)),
            "basemod_color_on": (None if hide_methylation else str(_BASEMOD_BINARY_ON_COLOR)),
            "basemod_color_off": (None if hide_methylation else str(_BASEMOD_BINARY_OFF_COLOR)),
            "gene_track": str(gene_track_path or ""),
            "igv_cmd": str(args.igv_cmd),
            "batch_file": str(batch_path.resolve()),
            "summary_tsv": str(summary_path.resolve()),
            "snapshot_dir": str(output_dir.resolve()),
            "jobs": summary_rows,
        }
        manifest_path.write_text(json.dumps(manifest_payload, indent=2), encoding="utf-8")

        logger.info("Wrote IGV batch script: %s", batch_path)
        logger.info("Wrote IGV summary TSV: %s", summary_path)
        logger.info("Wrote IGV manifest: %s", manifest_path)

        if getattr(args, "batch_only", False):
            logger.info("batch_only is set; skipping IGV execution.")
            return

        stdout_text, stderr_text = _run_igv_batch(
            str(args.igv_cmd),
            batch_path,
            snapshot_width=int(snapshot_width),
            snapshot_height=int(snapshot_height),
            env=igv_subprocess_env,
        )
        if stdout_text.strip():
            logger.debug("IGV stdout:\n%s", stdout_text)
        if stderr_text.strip():
            logger.debug("IGV stderr:\n%s", stderr_text)

        missing = [row["snapshot"] for row in summary_rows if not Path(str(row["snapshot"])).exists()]
        if missing:
            logger.warning(
                "IGV finished but %d snapshot files were not found. Check IGV logs and batch script.",
                len(missing),
            )
        else:
            logger.info("Wrote %d IGV snapshot(s) to %s", len(summary_rows), output_dir)
    finally:
        if temp_ctx is not None:
            temp_ctx.cleanup()
