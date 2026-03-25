from __future__ import annotations

import json
import logging
import threading
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from textwrap import shorten

import numpy as np
import pandas as pd
import pysam

from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam
from sniffcell.anno.variant_assignment import (
    _build_group_leaf_sets,
    _decode_linked_celltypes_from_row,
    _resolve_hierarchy_labels,
)

_THREAD_LOCAL_IO = threading.local()


def _path_cache_key(path: str | None) -> str | None:
    if path is None:
        return None
    return str(Path(path).expanduser().resolve())


def _get_thread_bam_handle(bam_path: str) -> pysam.AlignmentFile:
    key = _path_cache_key(bam_path)
    if key is None:
        raise ValueError("bam_path is required")
    cache = getattr(_THREAD_LOCAL_IO, "bam_handles", None)
    if cache is None:
        cache = {}
        _THREAD_LOCAL_IO.bam_handles = cache
    handle = cache.get(key)
    if handle is None:
        handle = pysam.AlignmentFile(key, "rb")
        cache[key] = handle
    return handle


def _get_thread_fasta_handle(reference_path: str | None) -> pysam.FastaFile | None:
    key = _path_cache_key(reference_path)
    if key is None:
        return None
    cache = getattr(_THREAD_LOCAL_IO, "fasta_handles", None)
    if cache is None:
        cache = {}
        _THREAD_LOCAL_IO.fasta_handles = cache
    handle = cache.get(key)
    if handle is None:
        handle = pysam.FastaFile(key)
        cache[key] = handle
    return handle


def _norm_chr(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    return text[3:] if text.lower().startswith("chr") else text


def _parse_support_read_names(value: object) -> list[str]:
    null_tokens = {"", "NA", "NAN", "NONE", "."}

    if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
        out: list[str] = []
        for item in value:
            text = str(item).strip()
            if text and text.upper() not in null_tokens:
                out.append(text)
        return out

    if isinstance(value, str):
        text = value.strip()
        if not text or text.upper() in null_tokens:
            return []
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = json.loads(text)
            except Exception:
                parsed = None
            if isinstance(parsed, list):
                return [str(item).strip() for item in parsed if str(item).strip() and str(item).strip().upper() not in null_tokens]
        for delim in (",", "|", ";"):
            if delim in text:
                return [t.strip() for t in text.split(delim) if t.strip() and t.strip().upper() not in null_tokens]
        return [text]

    try:
        if pd.isna(value):
            return []
    except TypeError:
        pass

    text = str(value).strip()
    if not text or text.upper() in null_tokens:
        return []
    return [text]


def _parse_pipe_values(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [token.strip() for token in text.split("|") if token.strip()]


def _load_kanpig_supporting_reads(path: str | None, sv_id: str) -> set[str]:
    if path is None:
        return set()
    mapping = pd.read_csv(path, sep="\t", header=None, names=["sv_id", "read_name"], dtype=str)
    subset = mapping.loc[mapping["sv_id"] == sv_id, "read_name"].dropna().astype(str)
    return {x.strip() for x in subset if x.strip()}


def _first_scalar(value: object) -> object:
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return value[0]
    return value


def _safe_info_get(info, key: str, default=None):
    try:
        value = info.get(key, default)
    except (KeyError, ValueError):
        return default
    if value is None:
        return default
    return value


def _build_sv_payload(record: pysam.VariantRecord) -> dict:
    sv_start = int(record.start)
    sv_end = int(record.stop)
    if sv_end <= sv_start:
        svlen_val = _first_scalar(record.info.get("SVLEN", 1))
        try:
            svlen_abs = abs(int(svlen_val))
        except (TypeError, ValueError):
            svlen_abs = 1
        sv_end = sv_start + max(1, svlen_abs)

    svtype = str(_safe_info_get(record.info, "SVTYPE", "SV"))
    support_names = set(_parse_support_read_names(_safe_info_get(record.info, "RNAMES", [])))
    return {
        "id": str(record.id),
        "chrom": str(record.chrom),
        "start": sv_start,
        "end": sv_end,
        "svtype": svtype,
        "variant_class": "SV",
        "svlen": _first_scalar(_safe_info_get(record.info, "SVLEN", pd.NA)),
        "supporting_reads": support_names,
    }


def _build_variant_payload_from_table_row(row: pd.Series) -> dict:
    start = int(pd.to_numeric(row.get("start"), errors="coerce"))
    end = int(pd.to_numeric(row.get("end"), errors="coerce"))
    if end <= start:
        end = start + 1
    group_a_reads = _parse_support_read_names(row.get("group_a_read_names", ""))
    group_b_reads = _parse_support_read_names(row.get("group_b_read_names", ""))
    support_names = set(group_a_reads) | set(group_b_reads)
    svlen = pd.to_numeric(row.get("change_size_bp", pd.NA), errors="coerce")
    return {
        "id": str(row.get("variant_id", row.get("id", ""))),
        "chrom": str(row.get("chrom", row.get("chr", ""))),
        "start": start,
        "end": end,
        "svtype": str(row.get("variant_subtype", row.get("sv_type", "SV"))),
        "variant_class": str(row.get("variant_class", "SV")),
        "svlen": (_first_scalar(svlen) if pd.notna(svlen) else pd.NA),
        "supporting_reads": support_names,
    }


@lru_cache(maxsize=8)
def _load_sv_payload_index(vcf_path: str) -> dict[str, dict]:
    key = _path_cache_key(vcf_path)
    if key is None:
        raise ValueError("vcf_path is required")

    payloads: dict[str, dict] = {}
    try:
        with pysam.VariantFile(key) as vf:
            for record in vf.fetch():
                rec_id = str(record.id)
                if not rec_id:
                    continue
                payloads[rec_id] = _build_sv_payload(record)
        return payloads
    except Exception:
        pass

    table = pd.read_csv(key, sep="\t")
    required = {"chrom", "start", "end"}
    if not required.issubset(set(table.columns)):
        raise ValueError(f"Variant table is missing required columns: {sorted(required - set(table.columns))}")
    if "variant_id" not in table.columns and "id" not in table.columns:
        raise ValueError("Variant table is missing required column 'variant_id' or 'id'")
    for _, row in table.iterrows():
        payload = _build_variant_payload_from_table_row(row)
        rec_id = str(payload["id"]).strip()
        if not rec_id:
            continue
        payloads[rec_id] = payload
    return payloads


def _get_sv_payload(vcf_path: str, sv_id: str) -> dict:
    payloads = _load_sv_payload_index(vcf_path)
    payload = payloads.get(str(sv_id))
    if payload is None:
        raise ValueError(f"SV ID '{sv_id}' was not found in VCF: {vcf_path}")
    out = dict(payload)
    out["supporting_reads"] = set(payload.get("supporting_reads", set()))
    return out


def _resolve_output_path(output: str, fmt: str) -> Path:
    out = Path(output)
    if out.suffix.lower() in {".png", ".pdf"}:
        return out
    return out.with_suffix(f".{fmt}")


@lru_cache(maxsize=8)
def _load_anno_manifest(anno_output: str) -> dict:
    anno_root = _path_cache_key(anno_output) or str(anno_output)
    manifest_path = Path(anno_root) / "anno_run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Could not find anno run manifest: {manifest_path}. "
            "Run `sniffcell anno` again (new versions write this file), or pass -i/-v/-r/-b/-a explicitly."
        )
    with manifest_path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid anno manifest format: {manifest_path}")
    return payload


def _resolve_viz_runtime_inputs(args, logger: logging.Logger) -> dict:
    anno_output = getattr(args, "anno_output", None)
    manifest = {}
    if anno_output:
        manifest = _load_anno_manifest(anno_output)
        logger.debug("Loaded anno run manifest from: %s", Path(anno_output) / "anno_run_manifest.json")

    manifest_inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}
    manifest_runtime = manifest.get("runtime", {}) if isinstance(manifest, dict) else {}
    manifest_outputs = manifest.get("outputs", {}) if isinstance(manifest, dict) else {}

    bam_path = args.input or manifest_inputs.get("bam")
    vcf_path = args.vcf or manifest_inputs.get("vcf") or manifest_inputs.get("variants")
    reference_path = args.reference or manifest_inputs.get("reference")
    bed_path = args.bed or manifest_inputs.get("bed")

    read_assignment = args.read_assignment
    if read_assignment is None:
        read_assignment = manifest_outputs.get("reads_classification")
        if read_assignment is None and anno_output:
            guessed = Path(anno_output) / "reads_classification.tsv"
            if guessed.exists():
                read_assignment = str(guessed)

    sv_assignment = manifest_outputs.get("sv_assignment")
    if sv_assignment is None and anno_output:
        guessed = Path(anno_output) / "sv_assignment.tsv"
        if guessed.exists():
            sv_assignment = str(guessed)

    if not bam_path or not vcf_path:
        raise ValueError(
            "viz needs BAM and VCF. Provide -i/-v, or provide --anno_output pointing to an anno folder with anno_run_manifest.json."
        )

    if args.output:
        output_path = _resolve_output_path(args.output, args.format)
    elif anno_output:
        output_path = Path(anno_output) / f"{args.sv_id}.viz.{args.format}"
    else:
        output_path = Path.cwd() / f"{args.sv_id}.viz.{args.format}"

    effective_window = int(args.window)
    use_manifest_window = (not bool(getattr(args, "exact_window", False)))
    if use_manifest_window and anno_output and int(args.window) == 5000 and "window" in manifest_runtime:
        try:
            effective_window = int(manifest_runtime["window"])
        except (TypeError, ValueError):
            pass

    assignment_window = effective_window
    if "window" in manifest_runtime:
        try:
            assignment_window = int(manifest_runtime["window"])
        except (TypeError, ValueError):
            assignment_window = effective_window

    return {
        "bam_path": str(bam_path),
        "vcf_path": str(vcf_path),
        "reference_path": (str(reference_path) if reference_path else None),
        "bed_path": (str(bed_path) if bed_path else None),
        "read_assignment_path": (str(read_assignment) if read_assignment else None),
        "sv_assignment_path": (str(sv_assignment) if sv_assignment else None),
        "output_path": output_path,
        "window": effective_window,
        "assignment_window": assignment_window,
        "kanpig_read_names": getattr(args, "kanpig_read_names", None),
    }


def _fetch_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    supporting_reads: set[str],
    max_reads: int,
    bam_handle: pysam.AlignmentFile | None = None,
    support_haplotype_only: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    def _normalize_haplotype_tag(value: object) -> int | None:
        if value is None:
            return None
        try:
            hp = int(value)
        except (TypeError, ValueError):
            text = str(value).strip()
            if not text:
                return None
            if text.upper().startswith("HP"):
                text = text[2:]
            try:
                hp = int(text)
            except (TypeError, ValueError):
                return None
        return hp if hp > 0 else None

    def _order_reads_for_plot(reads: pd.DataFrame) -> pd.DataFrame:
        if reads.empty:
            return reads.copy()
        ordered = reads.copy()
        if "haplotype" not in ordered.columns:
            ordered["haplotype"] = pd.Series(pd.NA, index=ordered.index, dtype="Int64")
        ordered["haplotype"] = pd.to_numeric(ordered["haplotype"], errors="coerce").astype("Int64")
        ordered["haplotype_label"] = np.where(
            ordered["haplotype"].notna(),
            "HP" + ordered["haplotype"].astype(str),
            "unphased",
        )
        ordered["_haplotype_missing"] = ordered["haplotype"].isna()
        ordered["_haplotype_sort"] = ordered["haplotype"].fillna(10**9).astype("int64")
        ordered = ordered.sort_values(
            ["_haplotype_missing", "_haplotype_sort", "is_supporting", "start", "end", "read_name"],
            ascending=[True, True, False, True, True, True],
            kind="stable",
            ignore_index=True,
        )
        return ordered.drop(columns=["_haplotype_missing", "_haplotype_sort"])

    rows: list[dict] = []
    bam = bam_handle if bam_handle is not None else _get_thread_bam_handle(bam_path)
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.reference_start is None or read.reference_end is None:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        r_start = max(int(read.reference_start), start)
        r_end = min(int(read.reference_end), end)
        if r_end <= r_start:
            continue

        qname = str(read.query_name)
        rows.append(
            {
                "read_name": qname,
                "start": r_start,
                "end": r_end,
                "is_supporting": qname in supporting_reads,
                "haplotype": _normalize_haplotype_tag(read.get_tag("HP") if read.has_tag("HP") else None),
            }
        )

    all_reads = pd.DataFrame(rows)
    if all_reads.empty:
        all_reads.attrs["applied_support_haplotype"] = None
        return all_reads, all_reads

    all_reads = _order_reads_for_plot(all_reads)
    applied_support_haplotype: int | None = None
    display_reads = all_reads
    if support_haplotype_only:
        support_df = all_reads.loc[all_reads["is_supporting"]].copy()
        phased_support_haplotypes = (
            pd.to_numeric(support_df.get("haplotype"), errors="coerce").dropna().astype(int).unique().tolist()
            if not support_df.empty
            else []
        )
        if len(phased_support_haplotypes) == 1:
            applied_support_haplotype = int(phased_support_haplotypes[0])
            same_haplotype_mask = all_reads["haplotype"].fillna(-1).astype("int64").eq(applied_support_haplotype)
            display_reads = all_reads.loc[same_haplotype_mask | all_reads["is_supporting"]].copy()

    if max_reads <= 0 or len(display_reads) <= max_reads:
        shown = display_reads.copy()
    else:
        supporting_df = display_reads[display_reads["is_supporting"]].copy()
        non_support_df = display_reads[~display_reads["is_supporting"]].copy()
        if len(supporting_df) >= max_reads:
            shown = supporting_df.iloc[:max_reads].copy()
        else:
            n_keep_non = max_reads - len(supporting_df)
            shown = pd.concat([supporting_df, non_support_df.iloc[:n_keep_non]], ignore_index=True)
        shown = _order_reads_for_plot(shown)
    shown.attrs["applied_support_haplotype"] = applied_support_haplotype
    all_reads.attrs["applied_support_haplotype"] = applied_support_haplotype
    return shown, all_reads


def _collect_large_indels_from_cigar(
    *,
    cigartuples: list[tuple[int, int]] | None,
    reference_start: int,
    region_start: int,
    region_end: int,
    min_indel_bp: int,
) -> list[dict[str, int | str]]:
    events: list[dict[str, int | str]] = []
    if not cigartuples or min_indel_bp <= 0:
        return events

    ref_pos = int(reference_start)
    for op, length in cigartuples:
        size = int(length)
        if op in {0, 7, 8}:  # M, =, X
            ref_pos += size
            continue
        if op == 1:  # insertion relative to reference
            if size >= min_indel_bp and region_start <= ref_pos <= region_end:
                events.append(
                    {
                        "event_type": "INS",
                        "start": int(ref_pos),
                        "end": int(ref_pos),
                        "pos": int(ref_pos),
                        "length": size,
                    }
                )
            continue
        if op == 2:  # deletion relative to reference
            del_start = int(ref_pos)
            del_end = int(ref_pos + size)
            if size >= min_indel_bp and del_end > region_start and del_start < region_end:
                events.append(
                    {
                        "event_type": "DEL",
                        "start": max(region_start, del_start),
                        "end": min(region_end, del_end),
                        "pos": int(del_start),
                        "length": size,
                    }
                )
            ref_pos = del_end
            continue
        if op == 3:  # ref skip
            ref_pos += size
            continue
        # op in {4,5,6,9,10}: soft/hard clip, pad, back
    return events


def _fetch_large_indels(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    read_names: set[str],
    min_indel_bp: int,
    bam_handle: pysam.AlignmentFile | None = None,
) -> pd.DataFrame:
    cols = ["read_name", "event_type", "start", "end", "pos", "length"]
    if min_indel_bp <= 0 or not read_names:
        return pd.DataFrame(columns=cols)

    bam = bam_handle if bam_handle is not None else _get_thread_bam_handle(bam_path)
    rows: list[dict[str, int | str]] = []
    pending = set(str(x) for x in read_names if str(x).strip())
    for read in bam.fetch(chrom, start, end):
        if not pending:
            break
        if read.is_unmapped or read.reference_start is None or read.reference_end is None:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        read_name = str(read.query_name)
        if read_name not in pending:
            continue
        for event in _collect_large_indels_from_cigar(
            cigartuples=read.cigartuples,
            reference_start=int(read.reference_start),
            region_start=int(start),
            region_end=int(end),
            min_indel_bp=int(min_indel_bp),
        ):
            rows.append({"read_name": read_name, **event})
        pending.discard(read_name)

    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows, columns=cols)


def _insertion_bar_width_bp(indel_len: int, region_span: int) -> float:
    clipped = max(0, min(int(indel_len), 4000))
    if clipped <= 0:
        return 0.0
    span = max(1, int(region_span))
    min_bp = max(18.0, 0.010 * span)
    max_bp = max(min_bp + 1.0, 0.18 * span)
    # Use a wide linear ramp so 200-600 bp insertions separate clearly, but
    # keep the bar bounded within a reasonable fraction of the shown window.
    scaled = 20.0 + (0.52 * clipped)
    return float(min(max_bp, max(min_bp, scaled)))


@lru_cache(maxsize=4)
def _load_ctdmr_table_cached(bed_path: str) -> pd.DataFrame:
    cols = ["chr", "start", "end", "best_group", "best_group_leaves", "best_dir", "name"]
    dmrs = pd.read_csv(bed_path, sep="\t")
    if dmrs.empty:
        return pd.DataFrame(columns=cols + ["label", "chr_norm"])
    if dmrs.columns[0].startswith("#"):
        dmrs = dmrs.rename(columns={dmrs.columns[0]: dmrs.columns[0].lstrip("#")})
    required = {"chr", "start", "end"}
    if not required.issubset(dmrs.columns):
        missing = sorted(required - set(dmrs.columns))
        raise ValueError(f"ctDMR file missing required columns: {missing}")

    dmrs = dmrs.copy()
    dmrs["start"] = pd.to_numeric(dmrs["start"], errors="coerce")
    dmrs["end"] = pd.to_numeric(dmrs["end"], errors="coerce")
    dmrs = dmrs.dropna(subset=["chr", "start", "end"])
    dmrs["start"] = dmrs["start"].astype(int)
    dmrs["end"] = dmrs["end"].astype(int)
    dmrs = dmrs[dmrs["end"] > dmrs["start"]]
    dmrs["chr_norm"] = dmrs["chr"].map(_norm_chr)

    for col in cols:
        if col not in dmrs.columns:
            dmrs[col] = ""
        dmrs[col] = dmrs[col].fillna("")

    dmrs["label"] = np.where(
        dmrs["best_group_leaves"].astype(str).str.strip().ne(""),
        dmrs["best_group_leaves"].astype(str),
        np.where(
            dmrs["best_group"].astype(str).str.strip().ne(""),
            dmrs["best_group"].astype(str),
            dmrs["name"].astype(str),
        ),
    )
    dmrs = dmrs.sort_values(["chr_norm", "start", "end"], kind="stable", ignore_index=True)
    mean_cols = [c for c in dmrs.columns if isinstance(c, str) and c.startswith("mean_")]
    keep_cols = cols + ["label", "chr_norm"] + [c for c in mean_cols if c not in cols]
    return dmrs[keep_cols]


def _read_ctdmrs(bed_path: str | None, chrom: str, start: int, end: int) -> pd.DataFrame:
    cols = ["chr", "start", "end", "best_group", "best_group_leaves", "best_dir", "name", "label", "chr_norm"]
    if bed_path is None:
        return pd.DataFrame(columns=cols)

    key = _path_cache_key(bed_path)
    if key is None:
        return pd.DataFrame(columns=cols)
    dmrs = _load_ctdmr_table_cached(key)
    if dmrs.empty:
        return dmrs.copy()

    chrom_norm = _norm_chr(chrom)
    out = dmrs[(dmrs["chr_norm"] == chrom_norm) & (dmrs["start"] < end) & (dmrs["end"] > start)]
    if out.empty:
        return pd.DataFrame(columns=list(dmrs.columns))
    return out.copy()


def _summarize_ctdmr_overlap(
    dmrs: pd.DataFrame,
    reads: pd.DataFrame,
    sv_start: int,
    sv_end: int,
) -> pd.DataFrame:
    cols = [
        "chr",
        "start",
        "end",
        "label",
        "best_group",
        "best_group_leaves",
        "best_dir",
        "overlaps_sv_core",
        "sv_core_overlap_bp",
        "supporting_read_overlap_count",
        "non_supporting_read_overlap_count",
        "read_overlap_count",
    ]
    if dmrs.empty:
        return pd.DataFrame(columns=cols)
    if reads.empty:
        out = dmrs[["chr", "start", "end", "label", "best_group", "best_group_leaves", "best_dir"]].copy()
        out["overlaps_sv_core"] = (out["start"] < sv_end) & (out["end"] > sv_start)
        out["sv_core_overlap_bp"] = np.maximum(0, np.minimum(out["end"], sv_end) - np.maximum(out["start"], sv_start))
        out["supporting_read_overlap_count"] = 0
        out["non_supporting_read_overlap_count"] = 0
        out["read_overlap_count"] = 0
        return out[cols]

    read_start = reads["start"].to_numpy(np.int64)
    read_end = reads["end"].to_numpy(np.int64)
    read_support = reads["is_supporting"].to_numpy(bool)

    rows = []
    for row in dmrs.itertuples(index=False):
        ov = (read_start < int(row.end)) & (read_end > int(row.start))
        n_all = int(np.count_nonzero(ov))
        n_support = int(np.count_nonzero(ov & read_support))
        n_non_support = n_all - n_support
        sv_overlap_bp = int(max(0, min(int(row.end), sv_end) - max(int(row.start), sv_start)))
        rows.append(
            {
                "chr": row.chr,
                "start": int(row.start),
                "end": int(row.end),
                "label": str(row.label),
                "best_group": str(row.best_group),
                "best_group_leaves": str(row.best_group_leaves),
                "best_dir": str(row.best_dir),
                "overlaps_sv_core": sv_overlap_bp > 0,
                "sv_core_overlap_bp": sv_overlap_bp,
                "supporting_read_overlap_count": n_support,
                "non_supporting_read_overlap_count": n_non_support,
                "read_overlap_count": n_all,
            }
        )
    return pd.DataFrame(rows, columns=cols)


@lru_cache(maxsize=4)
def _load_read_assignment_table_cached(path: str) -> pd.DataFrame:
    assignment = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        dtype={
            "code": "string",
            "code_order": "string",
            "best_group": "string",
            "best_group_leaves": "string",
            "other_group": "string",
            "other_group_leaves": "string",
        },
    )
    if assignment.empty:
        assignment["read_name"] = pd.Series(dtype="string")
        assignment["chr_norm"] = pd.Series(dtype="string")
        return assignment

    assignment = assignment.copy()
    assignment["read_name"] = assignment.index.astype(str)

    for col in ["chr", "start", "end", "code", "code_order", "best_group", "best_group_leaves", "other_group", "other_group_leaves"]:
        if col not in assignment.columns:
            assignment[col] = pd.NA

    assignment["chr"] = assignment["chr"].astype("string")
    assignment["chr_norm"] = assignment["chr"].map(_norm_chr).astype("string")
    assignment["start"] = pd.to_numeric(assignment["start"], errors="coerce").astype("Int64")
    assignment["end"] = pd.to_numeric(assignment["end"], errors="coerce").astype("Int64")

    if "is_best_group" not in assignment.columns:
        assignment["is_best_group"] = False
    else:
        assignment["is_best_group"] = (
            assignment["is_best_group"]
            .map(lambda x: str(x).strip().lower() in {"1", "true", "t", "yes"})
            .fillna(False)
            .astype(bool)
        )
    return assignment


def _load_read_assignment_table(path: str | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    key = _path_cache_key(path)
    if key is None:
        return pd.DataFrame()
    return _load_read_assignment_table_cached(key)


@lru_cache(maxsize=8)
def _load_sv_assignment_table_cached(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _load_sv_assignment_table(path: str | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    key = _path_cache_key(path)
    if key is None:
        return pd.DataFrame()
    return _load_sv_assignment_table_cached(key)


def _get_sv_assignment_row(path: str | None, sv_id: str) -> pd.Series | None:
    table = _load_sv_assignment_table(path)
    if table.empty or "id" not in table.columns:
        return None
    subset = table[table["id"].astype(str) == str(sv_id)]
    if subset.empty:
        return None
    return subset.iloc[0]


def _decode_read_assignment_rows(evidence: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "read_name",
        "chr",
        "chr_norm",
        "start",
        "end",
        "assigned_celltypes",
        "assigned_schema",
        "code",
        "code_order",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=cols)

    group_leaf_sets = _build_group_leaf_sets(evidence)
    rows = []
    for row in evidence.itertuples(index=False):
        row_s = pd.Series(row._asdict())
        schema, linked = _decode_linked_celltypes_from_row(row_s)
        if not linked:
            continue
        resolved = _resolve_hierarchy_labels(linked, group_leaf_sets)
        if not resolved:
            continue
        rows.append(
            {
                "read_name": str(row_s.get("read_name", "")),
                "chr": str(row_s.get("chr", "")),
                "chr_norm": _norm_chr(row_s.get("chr", "")),
                "start": int(row_s.get("start")),
                "end": int(row_s.get("end")),
                "assigned_celltypes": "|".join(resolved),
                "assigned_schema": str(schema),
                "code": str(row_s.get("code", "")),
                "code_order": str(row_s.get("code_order", "")),
            }
        )
    return pd.DataFrame(rows, columns=cols)


def _collect_supporting_assignment_evidence(
    assignment_df: pd.DataFrame,
    supporting_reads: set[str],
    sv_chrom: str,
    sv_start: int,
    sv_end: int,
    link_window: int,
    clip_start: int | None = None,
    clip_end: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if assignment_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    support_list = sorted(str(x) for x in supporting_reads if str(x).strip())
    if not support_list:
        return pd.DataFrame(), pd.DataFrame()

    support_idx = assignment_df.index.intersection(pd.Index(support_list))
    if support_idx.empty:
        return pd.DataFrame(), pd.DataFrame()

    evidence = assignment_df.loc[support_idx].copy()
    if "chr_norm" not in evidence.columns:
        evidence["chr_norm"] = evidence["chr"].map(_norm_chr)

    if not pd.api.types.is_integer_dtype(evidence["start"].dtype):
        evidence["start"] = pd.to_numeric(evidence["start"], errors="coerce").astype("Int64")
    if not pd.api.types.is_integer_dtype(evidence["end"].dtype):
        evidence["end"] = pd.to_numeric(evidence["end"], errors="coerce").astype("Int64")

    evidence = evidence.dropna(subset=["chr_norm", "start", "end"])
    if evidence.empty:
        return pd.DataFrame(), pd.DataFrame()

    evidence["start"] = evidence["start"].astype(int)
    evidence["end"] = evidence["end"].astype(int)

    sv_chrom_norm = _norm_chr(sv_chrom)
    same_chr = evidence["chr_norm"].eq(sv_chrom_norm)
    overlap_padded = (int(sv_start) < (evidence["end"] + int(link_window))) & (int(sv_end) > (evidence["start"] - int(link_window)))
    overlap_core = (int(sv_start) <= evidence["end"]) & (int(sv_end) >= evidence["start"])
    in_sv_window = same_chr & overlap_padded & (~overlap_core)
    evidence = evidence[in_sv_window].copy()

    if clip_start is not None and clip_end is not None:
        evidence = evidence[(evidence["start"] < int(clip_end)) & (evidence["end"] > int(clip_start))].copy()

    evidence = evidence[evidence["read_name"].isin(support_list)].copy()
    if evidence.empty:
        return evidence, pd.DataFrame()

    decoded_detail = _decode_read_assignment_rows(evidence)
    return evidence, decoded_detail


def _summarize_supporting_read_assignments(
    assignment_df: pd.DataFrame,
    supporting_reads: set[str],
    sv_chrom: str,
    sv_start: int,
    sv_end: int,
    window: int,
    region_start: int,
    region_end: int,
    assignment_available: bool,
    clip_to_region: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_cols = [
        "read_name",
        "is_supporting",
        "assignment_status",
        "is_assigned",
        "assigned_celltypes",
        "assigned_celltype_counts",
        "n_candidate_rows",
        "n_assignment_rows",
        "n_assignment_regions",
    ]
    detail_cols = [
        "read_name",
        "chr",
        "chr_norm",
        "start",
        "end",
        "assigned_celltypes",
        "assigned_schema",
        "code",
        "code_order",
    ]

    support_list = sorted(str(x) for x in supporting_reads if str(x).strip())
    if not support_list:
        return pd.DataFrame(columns=summary_cols), pd.DataFrame(columns=detail_cols)

    if assignment_df.empty:
        status = "unassigned_no_assignment_file" if not assignment_available else "unassigned_no_assignment_rows"
        summary = pd.DataFrame(
            {
                "read_name": support_list,
                "is_supporting": True,
                "assignment_status": status,
                "is_assigned": False,
                "assigned_celltypes": "",
                "assigned_celltype_counts": "",
                "n_candidate_rows": 0,
                "n_assignment_rows": 0,
                "n_assignment_regions": 0,
            }
        )
        return summary, pd.DataFrame(columns=detail_cols)

    evidence, decoded_detail = _collect_supporting_assignment_evidence(
        assignment_df=assignment_df,
        supporting_reads=supporting_reads,
        sv_chrom=sv_chrom,
        sv_start=sv_start,
        sv_end=sv_end,
        link_window=window,
        clip_start=(region_start if clip_to_region else None),
        clip_end=(region_end if clip_to_region else None),
    )
    if evidence.empty:
        status = "unassigned_no_overlap_rows"
        summary = pd.DataFrame(
            {
                "read_name": support_list,
                "is_supporting": True,
                "assignment_status": status,
                "is_assigned": False,
                "assigned_celltypes": "",
                "assigned_celltype_counts": "",
                "n_candidate_rows": 0,
                "n_assignment_rows": 0,
                "n_assignment_regions": 0,
            }
        )
        return summary, pd.DataFrame(columns=detail_cols)

    raw_counts = evidence.groupby("read_name", sort=False).size().to_dict()
    by_read_counts: dict[str, dict[str, int]] = {}
    by_read_majority: dict[str, str] = {}
    decoded_n_by_read: dict[str, int] = {}
    decoded_regions_by_read: dict[str, int] = {}
    if not decoded_detail.empty:
        decoded_n_by_read = (
            decoded_detail.groupby("read_name", sort=False)
            .size()
            .astype(int)
            .to_dict()
        )
        decoded_regions_by_read = (
            decoded_detail[["read_name", "chr_norm", "start", "end"]]
            .drop_duplicates(ignore_index=True)
            .groupby("read_name", sort=False)
            .size()
            .astype(int)
            .to_dict()
        )
        per_read_link_counts = (
            decoded_detail
            .groupby(["read_name", "assigned_celltypes"], sort=False)
            .size()
            .rename("count")
            .reset_index()
            .sort_values(["read_name", "count", "assigned_celltypes"], ascending=[True, False, True], kind="stable")
        )
        winners = per_read_link_counts.drop_duplicates(subset=["read_name"], keep="first")
        by_read_majority = {str(r.read_name): str(r.assigned_celltypes) for r in winners.itertuples(index=False)}

        for read_name, read_df in per_read_link_counts.groupby("read_name", sort=False):
            counter: dict[str, int] = defaultdict(int)
            for token, cnt in zip(read_df["assigned_celltypes"], read_df["count"]):
                for ct in [x.strip() for x in str(token).split("|") if x.strip()]:
                    counter[ct] += int(cnt)
            by_read_counts[str(read_name)] = dict(counter)

    rows = []
    for read_name in support_list:
        raw_n = int(raw_counts.get(read_name, 0))
        ct_counter = by_read_counts.get(read_name, {})
        if ct_counter:
            ranked = sorted(ct_counter.items(), key=lambda x: (-x[1], x[0]))
            majority_link = by_read_majority.get(read_name, "")
            majority_celltypes = [ct for ct in str(majority_link).split("|") if ct.strip()]
            celltypes = majority_celltypes if majority_celltypes else [ct for ct, _ in ranked]
            ct_count_str = ";".join(f"{ct}:{cnt}" for ct, cnt in ranked)
            decoded_n = int(decoded_n_by_read.get(read_name, 0))
            decoded_regions_n = int(decoded_regions_by_read.get(read_name, 0))
            rows.append(
                {
                    "read_name": read_name,
                    "is_supporting": True,
                    "assignment_status": "assigned",
                    "is_assigned": True,
                    "assigned_celltypes": "|".join(celltypes),
                    "assigned_celltype_counts": ct_count_str,
                    "n_candidate_rows": raw_n,
                    "n_assignment_rows": decoded_n,
                    "n_assignment_regions": decoded_regions_n,
                }
            )
        else:
            status = "unassigned_unresolved_code" if raw_n > 0 else "unassigned_no_overlap_rows"
            rows.append(
                {
                    "read_name": read_name,
                    "is_supporting": True,
                    "assignment_status": status,
                    "is_assigned": False,
                    "assigned_celltypes": "",
                    "assigned_celltype_counts": "",
                    "n_candidate_rows": raw_n,
                    "n_assignment_rows": 0,
                    "n_assignment_regions": 0,
                }
            )

    summary_df = pd.DataFrame(rows, columns=summary_cols)
    summary_df = summary_df.sort_values(["is_assigned", "read_name"], ascending=[False, True], kind="stable", ignore_index=True)
    return summary_df, decoded_detail


def _build_linked_ctdmr_callouts(
    *,
    bed_path: str | None,
    decoded_assignment_df: pd.DataFrame,
    sv_assignment_row: pd.Series | None,
    region_start: int,
    region_end: int,
    max_callouts: int = 4,
) -> pd.DataFrame:
    base_cols = [
        "chr",
        "start",
        "end",
        "best_group",
        "best_group_leaves",
        "best_dir",
        "name",
        "label",
        "chr_norm",
        "callout_side",
        "callout_support_count",
        "callout_support_reads",
        "callout_assigned_celltypes",
        "callout_distance_bp",
    ]
    if decoded_assignment_df.empty:
        return pd.DataFrame(columns=base_cols)

    winning_celltypes: list[str] = []
    if sv_assignment_row is not None:
        winning_celltypes = _parse_pipe_values(sv_assignment_row.get("linked_celltypes", ""))
        if not winning_celltypes:
            winning_celltypes = _parse_pipe_values(sv_assignment_row.get("majority_linked_celltypes", ""))
        if not winning_celltypes:
            primary = str(sv_assignment_row.get("primary_celltype", "")).strip()
            if primary:
                winning_celltypes = [primary]

    if winning_celltypes:
        winners = set(winning_celltypes)
        selected = decoded_assignment_df[
            decoded_assignment_df["assigned_celltypes"].astype(str).map(
                lambda token: any(ct in winners for ct in _parse_pipe_values(token))
            )
        ].copy()
    else:
        selected = decoded_assignment_df.copy()

    if selected.empty:
        return pd.DataFrame(columns=base_cols)

    region_counts = (
        selected.groupby(["chr", "chr_norm", "start", "end"], sort=False)
        .agg(
            callout_support_count=("read_name", lambda s: int(pd.Index(s.astype(str)).nunique())),
            callout_support_reads=("read_name", lambda s: "|".join(sorted(pd.Index(s.astype(str)).unique()))),
            callout_assigned_celltypes=("assigned_celltypes", lambda s: "|".join(sorted({x for token in s.astype(str) for x in _parse_pipe_values(token)}))),
        )
        .reset_index()
    )
    if region_counts.empty:
        return pd.DataFrame(columns=base_cols)

    region_counts["callout_side"] = np.where(
        region_counts["end"] <= int(region_start),
        "left",
        np.where(region_counts["start"] >= int(region_end), "right", "inside"),
    )
    region_counts = region_counts[region_counts["callout_side"].isin(["left", "right"])].copy()
    if region_counts.empty:
        return pd.DataFrame(columns=base_cols)

    region_counts["callout_distance_bp"] = np.where(
        region_counts["callout_side"].eq("left"),
        int(region_start) - region_counts["end"].astype(int),
        region_counts["start"].astype(int) - int(region_end),
    ).astype(int)

    full_dmrs = _load_ctdmr_table_cached(_path_cache_key(bed_path)) if (_path_cache_key(bed_path) is not None) else pd.DataFrame()
    if full_dmrs.empty:
        out = region_counts.copy()
        for col in ["best_group", "best_group_leaves", "best_dir", "name", "label"]:
            out[col] = ""
        return out[base_cols]

    merged = full_dmrs.merge(
        region_counts,
        on=["chr", "chr_norm", "start", "end"],
        how="inner",
    )
    if merged.empty:
        out = region_counts.copy()
        for col in ["best_group", "best_group_leaves", "best_dir", "name", "label"]:
            out[col] = ""
        return out[base_cols]

    merged = merged.sort_values(
        ["callout_support_count", "callout_distance_bp", "callout_side", "start", "end"],
        ascending=[False, True, True, True, True],
        kind="stable",
        ignore_index=True,
    )
    if max_callouts > 0 and len(merged) > max_callouts:
        keep_parts = []
        per_side_cap = max(1, max_callouts // 2)
        for side in ("left", "right"):
            side_df = merged[merged["callout_side"] == side].head(per_side_cap)
            if not side_df.empty:
                keep_parts.append(side_df)
        kept = pd.concat(keep_parts, ignore_index=True) if keep_parts else merged.head(max_callouts)
        if len(kept) < max_callouts:
            remainder = merged.loc[~merged.index.isin(kept.index)].head(max_callouts - len(kept))
            kept = pd.concat([kept, remainder], ignore_index=True)
        merged = kept.drop_duplicates(
            subset=["chr", "start", "end", "callout_side"],
            keep="first",
            ignore_index=True,
        )
    return merged[base_cols + [c for c in merged.columns if isinstance(c, str) and c.startswith("mean_") and c not in base_cols]]


def _extend_region_to_first_informative_ctdmr(
    region_start: int,
    region_end: int,
    linked_ctdmr_callouts: pd.DataFrame,
) -> tuple[int, int, pd.Series | None]:
    if linked_ctdmr_callouts.empty:
        return int(region_start), int(region_end), None

    ranked = linked_ctdmr_callouts.sort_values(
        ["callout_distance_bp", "callout_support_count", "start", "end"],
        ascending=[True, False, True, True],
        kind="stable",
        ignore_index=True,
    )
    row = ranked.iloc[0]
    side = str(row.get("callout_side", "")).strip().lower()

    new_start = int(region_start)
    new_end = int(region_end)
    region_width = max(1, int(region_end) - int(region_start))
    dmr_width = max(1, int(row["end"]) - int(row["start"]))
    edge_pad_bp = max(500, min(5000, max(region_width // 8, dmr_width // 2)))
    if side == "left":
        new_start = min(new_start, int(row["start"]) - edge_pad_bp)
    elif side == "right":
        new_end = max(new_end, int(row["end"]) + edge_pad_bp)
    if new_start < 0:
        new_start = 0
    return new_start, new_end, row


def _expand_interval_for_visibility(
    start: int,
    end: int,
    *,
    region_start: int,
    region_end: int,
    min_span_bp: int,
) -> tuple[int, int]:
    start_i = int(start)
    end_i = int(end)
    if end_i <= start_i:
        end_i = start_i + 1

    span_bp = end_i - start_i
    target_span = max(span_bp, int(min_span_bp))
    if target_span <= span_bp:
        return start_i, end_i

    center = 0.5 * (start_i + end_i)
    new_start = int(np.floor(center - (target_span / 2.0)))
    new_end = int(np.ceil(center + (target_span / 2.0)))

    if new_start < int(region_start):
        shift = int(region_start) - new_start
        new_start += shift
        new_end += shift
    if new_end > int(region_end):
        shift = new_end - int(region_end)
        new_start -= shift
        new_end -= shift

    new_start = max(int(region_start), new_start)
    new_end = min(int(region_end), max(new_start + 1, new_end))
    return new_start, new_end


def _estimate_ctdmr_label_half_width_bp(
    *,
    text: str,
    region_span_bp: int,
    font_size: float,
) -> float:
    span = max(1.0, float(region_span_bp))
    n_chars = max(2, len(str(text)))
    # Approximate label width in data coordinates for collision detection.
    return max(
        6.0,
        span * (0.0055 + (0.0017 * n_chars) + (0.00011 * max(0.0, font_size - 8.0))),
    )


def _assign_ctdmr_label_lanes(
    entries: list[dict[str, object]],
    *,
    region_span_bp: int,
    font_size: float,
    lane_gap_bp: float | None = None,
) -> tuple[list[int], int]:
    if not entries:
        return [], 1

    gap_bp = float(lane_gap_bp) if lane_gap_bp is not None else max(8.0, 0.002 * max(1, region_span_bp))
    indexed: list[tuple[int, float, float]] = []
    for idx, entry in enumerate(entries):
        x_center = float(entry["x_center"])
        half_width = _estimate_ctdmr_label_half_width_bp(
            text=str(entry["text"]),
            region_span_bp=region_span_bp,
            font_size=font_size,
        )
        indexed.append((idx, x_center - half_width, x_center + half_width))

    indexed.sort(key=lambda item: (item[1], item[2]))
    lane_end_positions: list[float] = []
    assigned = [0] * len(entries)
    for idx, est_start, est_end in indexed:
        lane_idx = 0
        while lane_idx < len(lane_end_positions) and est_start <= (lane_end_positions[lane_idx] + gap_bp):
            lane_idx += 1
        if lane_idx == len(lane_end_positions):
            lane_end_positions.append(est_end)
        else:
            lane_end_positions[lane_idx] = est_end
        assigned[idx] = lane_idx
    return assigned, max(1, len(lane_end_positions))


def _compute_supporting_read_ctdmr_methylation(
    *,
    sv_id: str,
    bam_path: str,
    reference_path: str | None,
    dmrs: pd.DataFrame,
    support_assignment_df: pd.DataFrame,
    decoded_assignment_df: pd.DataFrame,
    logger: logging.Logger,
    bam_handle: pysam.AlignmentFile | None = None,
    fasta_handle: pysam.FastaFile | None = None,
) -> pd.DataFrame:
    cols = [
        "sv_id",
        "read_name",
        "assignment_status",
        "is_assigned",
        "assigned_celltypes",
        "dmr_assigned_celltypes",
        "dmr_was_assigned",
        "chr",
        "start",
        "end",
        "label",
        "best_group",
        "best_group_leaves",
        "best_dir",
        "mean_methylation",
        "n_cpg_observed",
        "n_cpg_in_dmr",
    ]
    if support_assignment_df.empty or dmrs.empty:
        return pd.DataFrame(columns=cols)
    if reference_path is None:
        logger.warning(
            "Skipping per-read methylation on ctDMRs because --reference was not provided."
        )
        return pd.DataFrame(columns=cols)

    support_df = support_assignment_df.copy()
    support_df["read_name"] = support_df["read_name"].astype(str)
    support_lookup = support_df.set_index("read_name", drop=False)
    supporting_reads = list(support_lookup.index.unique())
    support_meta = (
        support_lookup[["assignment_status", "is_assigned", "assigned_celltypes"]]
        .to_dict(orient="index")
    )

    dmr_assign_map: dict[tuple[str, str, int, int], str] = {}
    if not decoded_assignment_df.empty:
        grouped = decoded_assignment_df.groupby(["read_name", "chr_norm", "start", "end"], sort=False)["assigned_celltypes"]
        for key, values in grouped:
            merged = []
            seen = set()
            for token in values:
                for ct in [x.strip() for x in str(token).split("|") if x.strip()]:
                    if ct not in seen:
                        seen.add(ct)
                        merged.append(ct)
            dmr_assign_map[(str(key[0]), str(key[1]), int(key[2]), int(key[3]))] = "|".join(merged)

    bam_for_mm = bam_handle if bam_handle is not None else _get_thread_bam_handle(bam_path)
    fa_for_mm = fasta_handle if fasta_handle is not None else _get_thread_fasta_handle(reference_path)

    per_chr_methyl: dict[str, tuple[pd.DataFrame, np.ndarray]] = {}
    for chr_name, chr_dmrs in dmrs.groupby("chr", sort=False):
        chr_text = str(chr_name)
        fetch_start = int(pd.to_numeric(chr_dmrs["start"], errors="coerce").min())
        fetch_end = int(pd.to_numeric(chr_dmrs["end"], errors="coerce").max())
        read_matrix = pd.DataFrame()
        cpg_positions = np.asarray([], dtype=np.int64)
        try:
            mm, cpgs = methyl_matrix_from_bam(
                bam_path,
                reference_path,
                chrom=chr_text,
                start=fetch_start,
                end=fetch_end,
                return_positions=True,
                read_name_whitelist=set(supporting_reads),
                bam_handle=bam_for_mm,
                fasta_handle=fa_for_mm,
            )
            cpg_positions = np.asarray(cpgs, dtype=np.int64)
            if not mm.empty:
                if isinstance(mm.index, pd.MultiIndex) and "read_name" in mm.index.names:
                    read_names = mm.index.get_level_values("read_name").astype(str)
                else:
                    read_names = mm.index.astype(str)
                mm2 = mm.copy()
                mm2["_read_name"] = read_names
                value_cols = [c for c in mm2.columns if c != "_read_name"]
                if value_cols:
                    read_matrix = mm2.groupby("_read_name", sort=False)[value_cols].mean()
        except Exception:
            logger.exception(
                "Failed methylation extraction for %s:%d-%d; writing NA values.",
                chr_text,
                fetch_start,
                fetch_end,
            )
        per_chr_methyl[chr_text] = (read_matrix, cpg_positions)

    out_rows = []
    for dmr in dmrs.itertuples(index=False):
        dmr_chr = str(dmr.chr)
        dmr_start = int(dmr.start)
        dmr_end = int(dmr.end)
        dmr_chr_norm = _norm_chr(dmr_chr)
        read_matrix, cpg_positions = per_chr_methyl.get(
            dmr_chr, (pd.DataFrame(), np.asarray([], dtype=np.int64))
        )
        in_dmr = (cpg_positions >= dmr_start) & (cpg_positions < dmr_end)
        dmr_cols = [int(x) for x in cpg_positions[in_dmr]]
        dmr_cpg_count = int(len(dmr_cols))

        read_mean = pd.Series(dtype="float64")
        read_n_obs = pd.Series(dtype="int64")
        if dmr_cols and (not read_matrix.empty):
            dmr_matrix = read_matrix.loc[:, dmr_cols]
            read_mean = dmr_matrix.mean(axis=1, skipna=True)
            read_n_obs = dmr_matrix.notna().sum(axis=1).astype(int)

        for read_name in supporting_reads:
            base = support_meta.get(read_name, {})
            mean_val = read_mean.get(read_name, np.nan)
            n_obs_val = int(read_n_obs.get(read_name, 0)) if read_name in read_n_obs.index else 0

            dmr_assigned = dmr_assign_map.get((read_name, dmr_chr_norm, dmr_start, dmr_end), "")
            out_rows.append(
                {
                    "sv_id": sv_id,
                    "read_name": read_name,
                    "assignment_status": str(base.get("assignment_status", "")),
                    "is_assigned": bool(base.get("is_assigned", False)),
                    "assigned_celltypes": str(base.get("assigned_celltypes", "")),
                    "dmr_assigned_celltypes": dmr_assigned,
                    "dmr_was_assigned": bool(dmr_assigned),
                    "chr": dmr_chr,
                    "start": dmr_start,
                    "end": dmr_end,
                    "label": str(dmr.label),
                    "best_group": str(dmr.best_group),
                    "best_group_leaves": str(dmr.best_group_leaves),
                    "best_dir": str(dmr.best_dir),
                    "mean_methylation": float(mean_val) if pd.notna(mean_val) else np.nan,
                    "n_cpg_observed": n_obs_val,
                    "n_cpg_in_dmr": dmr_cpg_count,
                }
            )

    return pd.DataFrame(out_rows, columns=cols)


def _build_methylation_heatmap_matrix(
    methyl_df: pd.DataFrame,
    max_reads: int = 35,
    max_dmrs: int = 20,
) -> pd.DataFrame:
    if methyl_df.empty:
        return pd.DataFrame()
    m = methyl_df.copy()
    m = m[pd.notna(m["mean_methylation"])].copy()
    if m.empty:
        return pd.DataFrame()

    m["dmr_key"] = m["label"].astype(str)
    duplicated = m["dmr_key"].duplicated(keep=False)
    m.loc[duplicated, "dmr_key"] = (
        m.loc[duplicated, "label"].astype(str)
        + "["
        + m.loc[duplicated, "start"].astype(str)
        + "-"
        + m.loc[duplicated, "end"].astype(str)
        + "]"
    )

    pivot = m.pivot_table(index="read_name", columns="dmr_key", values="mean_methylation", aggfunc="mean")
    if pivot.empty:
        return pd.DataFrame()

    read_order = (
        m[["read_name", "is_assigned"]]
        .drop_duplicates(ignore_index=True)
        .sort_values(["is_assigned", "read_name"], ascending=[False, True], kind="stable")
        ["read_name"]
        .tolist()
    )
    read_order = [r for r in read_order if r in pivot.index]
    if read_order:
        pivot = pivot.reindex(read_order)

    col_order = (
        m.groupby("dmr_key", sort=False)["mean_methylation"]
        .apply(lambda s: int(s.notna().sum()))
        .sort_values(ascending=False, kind="stable")
        .index.tolist()
    )
    col_order = [c for c in col_order if c in pivot.columns]
    if col_order:
        pivot = pivot.reindex(columns=col_order)

    if len(pivot.index) > max_reads:
        pivot = pivot.iloc[:max_reads, :]
    if len(pivot.columns) > max_dmrs:
        pivot = pivot.iloc[:, :max_dmrs]
    return pivot


def _build_dmr_methylation_stats_map(methyl_df: pd.DataFrame) -> dict[tuple[str, int, int], tuple[float, int]]:
    if methyl_df.empty:
        return {}
    stats = (
        methyl_df.groupby(["chr", "start", "end"], sort=False)
        .agg(
            mean_supporting_methylation=("mean_methylation", "mean"),
            supporting_reads_with_signal=("n_cpg_observed", lambda s: int(np.count_nonzero(np.asarray(s) > 0))),
        )
        .reset_index()
    )
    out: dict[tuple[str, int, int], tuple[float, int]] = {}
    for row in stats.itertuples(index=False):
        if pd.notna(row.mean_supporting_methylation):
            out[(str(row.chr), int(row.start), int(row.end))] = (
                float(row.mean_supporting_methylation),
                int(row.supporting_reads_with_signal),
            )
    return out


def _build_read_methylation_map(methyl_df: pd.DataFrame) -> dict[str, float]:
    if methyl_df.empty:
        return {}
    m = methyl_df[pd.notna(methyl_df["mean_methylation"])].copy()
    if m.empty:
        return {}
    stats = (
        m.groupby("read_name", sort=False)["mean_methylation"]
        .mean()
        .astype(float)
    )
    return {str(k): float(v) for k, v in stats.items() if pd.notna(v)}


def _reference_celltype_mean_columns(dmrs: pd.DataFrame) -> list[str]:
    if dmrs.empty:
        return []
    skip = {"mean_best_value", "mean_rest_value", "mean_margin"}
    all_cols = [
        c for c in dmrs.columns
        if isinstance(c, str) and c.startswith("mean_") and c not in skip
    ]
    if not all_cols:
        return []
    return all_cols


def _plot_sv_panel(
    *,
    sv: dict,
    shown_reads: pd.DataFrame,
    all_reads: pd.DataFrame,
    large_indels: pd.DataFrame,
    dmrs: pd.DataFrame,
    linked_ctdmr_callouts: pd.DataFrame,
    support_assignment_df: pd.DataFrame,
    methyl_df: pd.DataFrame,
    region_start: int,
    region_end: int,
    window: int,
    indel_min_bp: int,
    output_path: Path,
    dpi: int,
    applied_support_haplotype: int | None = None,
) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib import transforms
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch, Rectangle
    except ImportError as e:
        raise ImportError(
            "sniffcell viz requires matplotlib. Install it with `pip install matplotlib`."
        ) from e

    if dmrs.empty:
        dmr_panel_df = linked_ctdmr_callouts.copy()
    elif linked_ctdmr_callouts.empty:
        dmr_panel_df = dmrs.copy()
    else:
        dmr_panel_df = pd.concat([dmrs, linked_ctdmr_callouts], ignore_index=True, sort=False)

    ref_mean_cols = _reference_celltype_mean_columns(dmr_panel_df)
    ref_celltypes = [c[len("mean_"):] for c in ref_mean_cols]

    fig_height = max(6.4, 5.0 + 0.060 * max(1, len(shown_reads)) + 0.035 * max(1, len(ref_mean_cols)))
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(18, fig_height),
        sharex=True,
        constrained_layout=False,
        gridspec_kw={"height_ratios": [4.2, 1.8]},
    )
    fig.subplots_adjust(left=0.080, right=0.995, top=0.76, bottom=0.08, hspace=0.08)
    ax_reads, ax_dmrs = axes

    title_size = 17
    subtitle_size = 11
    axis_label_size = 16
    tick_label_size = 13

    sv_start = int(sv["start"])
    sv_end = int(sv["end"])
    chrom = str(sv["chrom"])

    assign_map: dict[str, bool] = {}
    assign_status_map: dict[str, str] = {}
    if not support_assignment_df.empty:
        assign_map = (
            support_assignment_df.set_index("read_name")["is_assigned"]
            .astype(bool)
            .to_dict()
        )
        assign_status_map = (
            support_assignment_df.set_index("read_name")["assignment_status"]
            .astype(str)
            .to_dict()
        )

    read_meth_map = _build_read_methylation_map(methyl_df)
    methyl_cmap = plt.cm.bwr
    low_methylation_color = methyl_cmap(0.0)
    mid_methylation_color = methyl_cmap(0.5)
    high_methylation_color = methyl_cmap(1.0)
    insertion_color = "#7b3294"
    deletion_edge = "#e66101"
    deletion_face = "#fff7bc"

    for ax in [ax_reads, ax_dmrs]:
        ax.axvspan(sv_start, sv_end, color="#3182bd", alpha=0.20, zorder=0)
        ax.axvline(sv_start, color="#225ea8", linewidth=1.2, alpha=0.85)
        ax.axvline(sv_end, color="#225ea8", linewidth=1.2, alpha=0.85)

    plotted_read_methylation = False
    read_track_map: dict[str, tuple[int, int, int, bool]] = {}
    distal_extension_reads: dict[str, set[str]] = {"left": set(), "right": set()}
    if not linked_ctdmr_callouts.empty:
        for row in linked_ctdmr_callouts.itertuples(index=False):
            side = str(getattr(row, "callout_side", "")).strip().lower()
            if side not in distal_extension_reads:
                continue
            distal_extension_reads[side].update(
                _parse_support_read_names(getattr(row, "callout_support_reads", ""))
            )
    if shown_reads.empty:
        ax_reads.text(0.01, 0.5, "No reads in selected window", transform=ax_reads.transAxes, va="center", fontsize=axis_label_size)
        ax_reads.set_ylim(0, 1)
    else:
        for i, row in enumerate(shown_reads.itertuples(index=False), start=1):
            read_name = str(row.read_name)
            r_start = int(row.start)
            r_end = int(row.end)
            is_supporting = bool(row.is_supporting)
            read_track_map[read_name] = (i, r_start, r_end, is_supporting)

            if bool(row.is_supporting):
                status = assign_status_map.get(read_name, "assigned")
                is_assigned = bool(assign_map.get(read_name, status == "assigned"))
                if is_assigned:
                    color = "#d62728"
                    marker = ">"
                    linestyle = "-"
                elif status.startswith("unassigned"):
                    color = "#ff7f0e"
                    marker = "X"
                    linestyle = "--"
                else:
                    color = "#ff7f0e"
                    marker = "D"
                    linestyle = ":"
                width = 2.6
                alpha = 0.95

                # Assignment pointer marker at read start.
                ax_reads.scatter(
                    [int(row.start)],
                    [i],
                    marker=marker,
                    s=130,
                    c=[color],
                    edgecolors="#111111",
                    linewidths=0.9,
                    zorder=8,
                )
            else:
                color = "#9e9e9e"
                width = 1.3
                alpha = 0.60
                linestyle = "-"
            ax_reads.hlines(i, r_start, r_end, color=color, linewidth=width, alpha=alpha, linestyles=linestyle, zorder=3)

            if is_supporting and read_name in distal_extension_reads["left"] and r_start > region_start:
                ax_reads.hlines(
                    i,
                    region_start,
                    r_start,
                    color=color,
                    linewidth=max(1.2, width - 0.8),
                    alpha=0.70,
                    linestyles=(0, (4, 3)),
                    zorder=4,
                )
            if is_supporting and read_name in distal_extension_reads["right"] and r_end < region_end:
                ax_reads.hlines(
                    i,
                    r_end,
                    region_end,
                    color=color,
                    linewidth=max(1.2, width - 0.8),
                    alpha=0.70,
                    linestyles=(0, (4, 3)),
                    zorder=4,
                )

        # Overlay ctDMR-specific methylation directly on supporting-read segments.
        if not methyl_df.empty:
            mdf = methyl_df[
                pd.notna(methyl_df["mean_methylation"])
                & (pd.to_numeric(methyl_df["n_cpg_observed"], errors="coerce").fillna(0) > 0)
            ].copy()
            for mrow in mdf.itertuples(index=False):
                read_name = str(mrow.read_name)
                if read_name not in read_track_map:
                    continue
                y, r_start, r_end, is_supporting = read_track_map[read_name]
                if not is_supporting:
                    continue
                seg_start = max(r_start, int(mrow.start))
                seg_end = min(r_end, int(mrow.end))
                if seg_end <= seg_start:
                    continue
                meth = float(np.clip(float(mrow.mean_methylation), 0.0, 1.0))
                seg_color = methyl_cmap(meth)
                # Black outline for visibility, then colored methylation segment on top.
                ax_reads.hlines(
                    y,
                    seg_start,
                    seg_end,
                    color="#111111",
                    linewidth=6.8,
                    alpha=0.95,
                    zorder=6,
                )
                ax_reads.hlines(
                    y,
                    seg_start,
                    seg_end,
                    color=seg_color,
                    linewidth=5.0,
                    alpha=0.92,
                    zorder=7,
                )
                mid_x = 0.5 * (seg_start + seg_end)
                ax_reads.scatter(
                    [mid_x],
                    [y],
                    marker="o",
                    s=26,
                    c=[seg_color],
                    edgecolors="#111111",
                    linewidths=0.35,
                    zorder=9,
                )
                plotted_read_methylation = True

        # Fallback single marker when ctDMR segment overlap is absent.
        if (not plotted_read_methylation) and read_meth_map:
            read_marker_offset = max(1, int(0.004 * max(1, region_end - region_start)))
            for read_name, meth in read_meth_map.items():
                if read_name not in read_track_map:
                    continue
                y, r_start, r_end, is_supporting = read_track_map[read_name]
                if not is_supporting:
                    continue
                marker_x = max(r_start, r_end - read_marker_offset)
                marker_color = methyl_cmap(float(np.clip(float(meth), 0.0, 1.0)))
                ax_reads.scatter(
                    [marker_x],
                    [y],
                    marker="o",
                    s=48,
                    c=[marker_color],
                    edgecolors="#000000",
                    linewidths=0.35,
                    zorder=9,
                )
                plotted_read_methylation = True

        if not large_indels.empty:
            annotate_indels = len(large_indels) <= 20
            for indel in large_indels.itertuples(index=False):
                read_name = str(indel.read_name)
                if read_name not in read_track_map:
                    continue
                y, _, _, _ = read_track_map[read_name]
                indel_len = int(indel.length)
                if str(indel.event_type) == "INS":
                    x = int(indel.pos)
                    bar_width_bp = _insertion_bar_width_bp(
                        indel_len=indel_len,
                        region_span=max(1, region_end - region_start),
                    )
                    half_width = 0.5 * bar_width_bp
                    x0 = max(float(region_start), float(x) - half_width)
                    x1 = min(float(region_end), float(x) + half_width)
                    if x1 <= x0:
                        x1 = x0 + 1.0
                    outer_height = 0.24
                    inner_height = 0.15
                    cap_half_height = 0.19
                    outer_rect = Rectangle(
                        (x0, y - (outer_height / 2.0)),
                        max(1.0, x1 - x0),
                        outer_height,
                        facecolor="#111111",
                        edgecolor="#111111",
                        linewidth=0.0,
                        alpha=0.96,
                        zorder=9,
                    )
                    inner_rect = Rectangle(
                        (x0, y - (inner_height / 2.0)),
                        max(1.0, x1 - x0),
                        inner_height,
                        facecolor=insertion_color,
                        edgecolor=insertion_color,
                        linewidth=0.0,
                        alpha=0.98,
                        zorder=10,
                    )
                    ax_reads.add_patch(outer_rect)
                    ax_reads.add_patch(inner_rect)
                    ax_reads.vlines(
                        [x0, x1],
                        [y - cap_half_height, y - cap_half_height],
                        [y + cap_half_height, y + cap_half_height],
                        color="#111111",
                        linewidth=1.9,
                        alpha=0.95,
                        zorder=11,
                    )
                    ax_reads.vlines(
                        [x0, x1],
                        [y - (cap_half_height - 0.03), y - (cap_half_height - 0.03)],
                        [y + (cap_half_height - 0.03), y + (cap_half_height - 0.03)],
                        color=insertion_color,
                        linewidth=1.0,
                        alpha=0.98,
                        zorder=12,
                    )
                    if annotate_indels:
                        ax_reads.text(
                            0.5 * (x0 + x1),
                            y + 0.28,
                            f"{indel_len}I",
                            ha="center",
                            va="bottom",
                            fontsize=8.8,
                            fontweight="bold",
                            color=insertion_color,
                            zorder=13,
                        )
                elif str(indel.event_type) == "DEL":
                    x0 = int(indel.start)
                    x1 = max(x0 + 1, int(indel.end))
                    rect = Rectangle(
                        (x0, y - 0.18),
                        max(1, x1 - x0),
                        0.36,
                        facecolor=deletion_face,
                        edgecolor=deletion_edge,
                        linewidth=1.1,
                        alpha=0.95,
                        zorder=10,
                    )
                    ax_reads.add_patch(rect)
                    if annotate_indels:
                        ax_reads.text(
                            0.5 * (x0 + x1),
                            y + 0.24,
                            f"{indel_len}D",
                            ha="center",
                            va="bottom",
                            fontsize=8,
                            color=deletion_edge,
                            rotation=90,
                            zorder=11,
                        )
        has_haplotype_groups = (
            ("haplotype_label" in shown_reads.columns)
            and shown_reads["haplotype_label"].astype(str).ne("unphased").any()
        )
        if has_haplotype_groups:
            group_rows: list[tuple[str, int, int]] = []
            prev_label: str | None = None
            group_start_y: int | None = None
            prev_y: int | None = None
            for y, row in enumerate(shown_reads.itertuples(index=False), start=1):
                label = str(getattr(row, "haplotype_label", "unphased"))
                if prev_label is None:
                    prev_label = label
                    group_start_y = y
                    prev_y = y
                    continue
                if label != prev_label:
                    assert group_start_y is not None and prev_y is not None
                    group_rows.append((prev_label, group_start_y, prev_y))
                    boundary_y = prev_y + 0.5
                    ax_reads.hlines(
                        boundary_y,
                        region_start,
                        region_end,
                        color="#bdbdbd",
                        linewidth=1.0,
                        linestyles=(0, (4, 3)),
                        alpha=0.85,
                        zorder=4,
                    )
                    prev_label = label
                    group_start_y = y
                prev_y = y
            if prev_label is not None and group_start_y is not None and prev_y is not None:
                group_rows.append((prev_label, group_start_y, prev_y))

            y_axis_transform = transforms.blended_transform_factory(ax_reads.transAxes, ax_reads.transData)
            for label, y0, y1 in group_rows:
                ax_reads.text(
                    -0.018,
                    0.5 * (y0 + y1),
                    label,
                    transform=y_axis_transform,
                    ha="right",
                    va="center",
                    fontsize=max(9, tick_label_size - 1),
                    fontweight="bold",
                    color="#4d4d4d",
                    clip_on=False,
                    zorder=12,
                )
        ax_reads.set_ylim(0, len(shown_reads) + 1)

    ax_reads.set_ylabel("Reads", fontsize=axis_label_size)
    ax_reads.set_yticks([])
    ax_reads.grid(axis="x", alpha=0.2)
    ax_reads.tick_params(axis="x", labelsize=tick_label_size)

    plotted_reference_methylation = False
    if dmr_panel_df.empty:
        ax_dmrs.text(0.01, 0.5, "No linked ctDMRs are available for this SV", transform=ax_dmrs.transAxes, va="center", fontsize=axis_label_size)
        ax_dmrs.set_ylim(0, 1)
        ax_dmrs.set_yticks([])
    elif not ref_mean_cols:
        ax_dmrs.text(
            0.01,
            0.5,
            "ctDMRs found, but no input cell-type methylation columns (mean_<celltype>) were found in BED.",
            transform=ax_dmrs.transAxes,
            va="center",
            fontsize=axis_label_size - 1,
        )
        ax_dmrs.set_ylim(0, 1)
        ax_dmrs.set_yticks([])
    else:
        dmrs_plot = dmrs.reset_index(drop=True)
        cmap = methyl_cmap
        n_celltypes = len(ref_mean_cols)
        min_dmr_display_bp = max(30, min(300, int(round(max(1, region_end - region_start) * 0.008))))
        if dmrs_plot.empty:
            ax_dmrs.text(
                0.01,
                0.5,
                "No ctDMRs overlap this window; linked distal ctDMR callouts are shown at the panel edge.",
                transform=ax_dmrs.transAxes,
                va="center",
                fontsize=axis_label_size - 1,
                bbox={"facecolor": "#ffffff", "edgecolor": "none", "alpha": 0.7},
            )
        for ct_idx, mean_col in enumerate(ref_mean_cols):
            y0 = ct_idx + 0.1
            y_center = ct_idx + 0.5
            label_entries: list[dict[str, object]] = []
            for _, row in dmrs_plot.iterrows():
                start = int(row["start"])
                end = int(row["end"])
                if end <= start:
                    continue
                draw_start, draw_end = _expand_interval_for_visibility(
                    start,
                    end,
                    region_start=region_start,
                    region_end=region_end,
                    min_span_bp=min_dmr_display_bp,
                )
                best_dir = str(row.get("best_dir", "")).lower()
                edge_color = "#2ca25f" if best_dir == "hyper" else "#756bb1"
                value = pd.to_numeric(row.get(mean_col, pd.NA), errors="coerce")
                if pd.notna(value):
                    v = float(np.clip(float(value), 0.0, 1.0))
                    face_color = cmap(v)
                    text_value = f"{v:.2f}"
                    text_color = "#111111"
                    plotted_reference_methylation = True
                else:
                    face_color = "#d9d9d9"
                    text_value = "NA"
                    text_color = "#111111"

                ax_dmrs.broken_barh(
                    [(draw_start, draw_end - draw_start)],
                    (y0, 0.8),
                    facecolors=face_color,
                    edgecolors=edge_color,
                    linewidth=0.8,
                    alpha=0.82,
                    zorder=2,
                )

                center_x = 0.5 * (draw_start + draw_end)
                label_entries.append(
                    {
                        "x_center": center_x,
                        "text": text_value,
                        "text_color": text_color,
                    }
                )

            label_font_size = max(10.0, tick_label_size - 1.0)
            lane_assignments, lane_count = _assign_ctdmr_label_lanes(
                label_entries,
                region_span_bp=max(1, region_end - region_start),
                font_size=label_font_size,
            )
            if lane_count <= 1:
                lane_positions = [y_center]
            else:
                lane_positions = np.linspace(y0 + 0.22, y0 + 0.58, lane_count).tolist()
            adjusted_font_size = max(8.2, label_font_size - (0.45 * max(0, lane_count - 1)))
            for entry, lane_idx in zip(label_entries, lane_assignments):
                ax_dmrs.text(
                    float(entry["x_center"]),
                    float(lane_positions[lane_idx]),
                    str(entry["text"]),
                    ha="center",
                    va="center",
                    fontsize=adjusted_font_size,
                    color=str(entry["text_color"]),
                    zorder=3,
                    bbox={"facecolor": (1, 1, 1, 0.22), "edgecolor": "none", "pad": 0.18},
                )

        if not linked_ctdmr_callouts.empty:
            edge_transform = transforms.blended_transform_factory(ax_dmrs.transAxes, ax_dmrs.transData)
            y_center_all = 0.5 * max(1, n_celltypes)
            side_slots = {
                "left": linked_ctdmr_callouts[linked_ctdmr_callouts["callout_side"] == "left"].reset_index(drop=True),
                "right": linked_ctdmr_callouts[linked_ctdmr_callouts["callout_side"] == "right"].reset_index(drop=True),
            }
            slot_width = 0.070
            slot_gap = 0.050
            for side, side_df in side_slots.items():
                for slot_idx, row in enumerate(side_df.itertuples(index=False)):
                    if side == "left":
                        x1 = -0.045 - (slot_idx * (slot_width + slot_gap))
                        x0 = x1 - slot_width
                        connector_x = x1
                        label_x = x0 - 0.018
                        label_rotation = 90
                    else:
                        x0 = 1.045 + (slot_idx * (slot_width + slot_gap))
                        x1 = x0 + slot_width
                        connector_x = x0
                        label_x = x1 + 0.018
                        label_rotation = 270

                    best_dir = str(getattr(row, "best_dir", "")).lower()
                    edge_color = "#2ca25f" if best_dir == "hyper" else "#756bb1"
                    for ct_idx, mean_col in enumerate(ref_mean_cols):
                        y0 = ct_idx + 0.1
                        value = pd.to_numeric(getattr(row, mean_col, pd.NA), errors="coerce")
                        if pd.notna(value):
                            v = float(np.clip(float(value), 0.0, 1.0))
                            face_color = cmap(v)
                            plotted_reference_methylation = True
                        else:
                            face_color = "#d9d9d9"
                        rect = Rectangle(
                            (x0, y0),
                            slot_width,
                            0.8,
                            transform=edge_transform,
                            facecolor=face_color,
                            edgecolor=edge_color,
                            linewidth=0.9,
                            alpha=0.96,
                            clip_on=False,
                            zorder=3,
                        )
                        ax_dmrs.add_patch(rect)

                    ax_dmrs.plot(
                        [0.0 if side == "left" else 1.0, connector_x],
                        [y_center_all, y_center_all],
                        transform=edge_transform,
                        color="#6b6b6b",
                        linewidth=1.3,
                        linestyle=(0, (4, 3)),
                        alpha=0.95,
                        clip_on=False,
                        zorder=2,
                    )
                    coord_label = f"{str(row.chr)}:{int(row.start) + 1}-{int(row.end)}"
                    support_n = int(getattr(row, "callout_support_count", 0))
                    support_label = "read" if support_n == 1 else "reads"
                    distance_bp = int(getattr(row, "callout_distance_bp", 0))
                    ax_dmrs.text(
                        label_x,
                        y_center_all,
                        f"{coord_label}\n{support_n} {support_label} | {distance_bp} bp",
                        transform=edge_transform,
                        ha="center",
                        va="center",
                        fontsize=max(7.8, tick_label_size - 4.5),
                        rotation=label_rotation,
                        color="#333333",
                        clip_on=False,
                        zorder=4,
                    )

        ax_dmrs.set_ylim(0, max(1, n_celltypes))
        ax_dmrs.set_yticks([i + 0.5 for i in range(n_celltypes)])
        ax_dmrs.set_yticklabels(ref_celltypes, fontsize=max(8, tick_label_size - 1))

    ax_dmrs.set_ylabel("Cell types", fontsize=axis_label_size)
    ax_dmrs.set_xlabel(f"{chrom} coordinate (bp)", fontsize=axis_label_size)
    ax_dmrs.grid(axis="x", alpha=0.2)
    ax_dmrs.tick_params(axis="x", labelsize=tick_label_size)
    ax_dmrs.set_xlim(region_start, region_end)

    n_support_listed = int(len(set(sv["supporting_reads"])))
    n_support_in_window = int(all_reads["is_supporting"].sum()) if not all_reads.empty else 0
    n_support_shown = int(shown_reads["is_supporting"].sum()) if not shown_reads.empty else 0
    n_assigned = int(support_assignment_df["is_assigned"].sum()) if not support_assignment_df.empty else 0
    n_total_support = int(len(support_assignment_df)) if not support_assignment_df.empty else n_support_listed
    n_unassigned = max(0, n_total_support - n_assigned)

    sv_len_signed: int | None
    try:
        sv_len_signed = int(float(sv.get("svlen", pd.NA)))
    except Exception:
        sv_len_signed = None
    if sv_len_signed is None:
        sv_len_signed = int(sv_end - sv_start)
    sv_len_abs = abs(int(sv_len_signed))
    sv_len_text = f"{sv_len_signed} bp (abs {sv_len_abs} bp)"

    sv_locus_igv = f"{chrom}:{sv_start + 1}-{sv_end}"
    window_locus_igv = f"{chrom}:{region_start + 1}-{region_end}"
    dmr_coords_igv = [
        f"{str(row.chr)}:{int(row.start) + 1}-{int(row.end)}"
        for row in dmrs.itertuples(index=False)
    ]
    callout_coords_igv = [
        f"{str(row.chr)}:{int(row.start) + 1}-{int(row.end)}"
        for row in linked_ctdmr_callouts.itertuples(index=False)
    ]
    dmr_preview = ", ".join(dmr_coords_igv[:8])
    if len(dmr_coords_igv) > 8:
        dmr_preview += f", ... (+{len(dmr_coords_igv) - 8} more)"
    if not dmr_preview:
        dmr_preview = "none"
    else:
        dmr_preview = shorten(dmr_preview, width=150, placeholder=" ...")
    callout_preview = ", ".join(callout_coords_igv[:4])
    if len(callout_coords_igv) > 4:
        callout_preview += f", ... (+{len(callout_coords_igv) - 4} more)"
    if not callout_preview:
        callout_preview = "none"
    else:
        callout_preview = shorten(callout_preview, width=120, placeholder=" ...")

    variant_class_label = str(sv.get("variant_class", "SV")).strip().upper() or "VARIANT"
    variant_type_label = str(sv.get("svtype", "")).strip() or variant_class_label

    title = (
        f"{variant_class_label} {sv['id']} ({variant_type_label}) at {chrom}:{sv_start + 1}-{sv_end}"
    )
    fig.suptitle(title, fontsize=title_size, y=0.992)
    subtitle = (
        f"size: {sv_len_text} | window +/-{window} bp | "
        f"supporting reads in VCF: {n_support_listed} | in BAM window: {n_support_in_window} | shown: {n_support_shown} | "
        f"assigned: {n_assigned} | unassigned: {n_unassigned} | "
        f"display haplotypes: {'HP' + str(applied_support_haplotype) if applied_support_haplotype is not None else 'all'}"
    )
    fig.text(0.01, 0.952, subtitle, fontsize=subtitle_size, ha="left", va="center")
    fig.text(0.01, 0.930, f"Locus (IGV): {sv_locus_igv} | Window (IGV): {window_locus_igv}", fontsize=10, ha="left", va="center")
    fig.text(0.01, 0.909, f"ctDMRs in-window (IGV, {len(dmr_coords_igv)}): {dmr_preview}", fontsize=9, ha="left", va="center")
    if callout_coords_igv:
        fig.text(
            0.01,
            0.888,
            f"winning linked ctDMR callouts outside window ({len(callout_coords_igv)}): {callout_preview}",
            fontsize=9,
            ha="left",
            va="center",
        )

    legend_handles = [
        Line2D([0], [0], color="#d62728", lw=2.6, linestyle="-", label="Supporting read (assigned)"),
        Line2D([0], [0], color="#ff7f0e", lw=2.6, linestyle="--", label="Supporting read (unassigned)"),
        Line2D([0], [0], color="#9e9e9e", lw=1.3, label="Other read"),
        Line2D([0], [0], linestyle="None", marker=">", markerfacecolor="#d62728", markeredgecolor="#111111", markersize=9, label="Assignment pointer"),
        Line2D([0], [0], color="#6b6b6b", lw=1.3, linestyle=(0, (4, 3)), label="Distal linked ctDMR extension"),
        Patch(facecolor=insertion_color, edgecolor="#111111", linewidth=1.1, label=f"Insertion >= {indel_min_bp} bp"),
        Patch(facecolor=deletion_face, edgecolor=deletion_edge, linewidth=1.1, label=f"Deletion >= {indel_min_bp} bp"),
        Line2D([0], [0], color="#6baed6", lw=5.0, label="Read methylation on ctDMR overlap"),
        Line2D([0], [0], linestyle="None", marker="o", markerfacecolor=low_methylation_color, markeredgecolor="#111111", markersize=7, label="Methylation 0.0"),
        Line2D([0], [0], linestyle="None", marker="o", markerfacecolor=mid_methylation_color, markeredgecolor="#111111", markersize=7, label="Methylation 0.5"),
        Line2D([0], [0], linestyle="None", marker="o", markerfacecolor=high_methylation_color, markeredgecolor="#111111", markersize=7, label="Methylation 1.0"),
        Patch(facecolor="#3182bd", alpha=0.20, label="Variant interval"),
        Patch(facecolor="#f7f7f7", edgecolor="#2ca25f", alpha=0.95, label="ctDMR hyper (edge)"),
        Patch(facecolor="#f7f7f7", edgecolor="#756bb1", alpha=0.95, label="ctDMR hypo/other (edge)"),
    ]
    ax_reads.legend(handles=legend_handles, loc="upper right", fontsize=tick_label_size, frameon=False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)


def viz_main(args) -> None:
    logger = logging.getLogger("sniffcell.viz")
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    max_reads = int(args.max_reads)
    if max_reads <= 0:
        raise ValueError("max_reads must be > 0")
    dpi = int(getattr(args, "dpi", 300))
    if dpi <= 0:
        raise ValueError("dpi must be > 0")
    indel_min_bp = int(getattr(args, "indel_min_bp", 40))
    if indel_min_bp < 0:
        raise ValueError("indel_min_bp must be >= 0")
    skip_methylation_overlay = bool(getattr(args, "skip_methylation_overlay", False))
    support_haplotype_only = bool(getattr(args, "support_haplotype_only", True))
    linked_ctdmr_mode = str(getattr(args, "linked_ctdmr_mode", "distal")).strip().lower()
    if linked_ctdmr_mode not in {"distal", "extend", "strict"}:
        raise ValueError("linked_ctdmr_mode must be one of: distal, extend, strict")

    resolved = _resolve_viz_runtime_inputs(args, logger)
    window = int(resolved["window"])
    if window < 0:
        raise ValueError("window must be >= 0")
    output_path = resolved["output_path"]
    logger.debug(
        "Resolved viz inputs: bam=%s vcf=%s ref=%s bed=%s read_assignment=%s sv_assignment=%s output=%s window=%d assignment_window=%d",
        resolved["bam_path"],
        resolved["vcf_path"],
        resolved["reference_path"],
        resolved["bed_path"],
        resolved["read_assignment_path"],
        resolved["sv_assignment_path"],
        output_path,
        window,
        int(resolved["assignment_window"]),
    )
    sv = _get_sv_payload(resolved["vcf_path"], args.sv_id)

    override_support = _load_kanpig_supporting_reads(resolved["kanpig_read_names"], args.sv_id)
    if override_support:
        sv["supporting_reads"] = override_support

    region_start = max(0, int(sv["start"]) - window)
    region_end = int(sv["end"]) + window

    assignment_df = _load_read_assignment_table(resolved["read_assignment_path"])
    sv_assignment_row = _get_sv_assignment_row(resolved["sv_assignment_path"], str(sv["id"]))
    support_assignment_df, decoded_assignment_df = _summarize_supporting_read_assignments(
        assignment_df,
        supporting_reads=set(sv["supporting_reads"]),
        sv_chrom=sv["chrom"],
        sv_start=int(sv["start"]),
        sv_end=int(sv["end"]),
        window=int(resolved["assignment_window"]),
        region_start=region_start,
        region_end=region_end,
        assignment_available=(resolved["read_assignment_path"] is not None),
        clip_to_region=False,
    )

    linked_ctdmr_candidates = _build_linked_ctdmr_callouts(
        bed_path=resolved["bed_path"],
        decoded_assignment_df=decoded_assignment_df,
        sv_assignment_row=sv_assignment_row,
        region_start=region_start,
        region_end=region_end,
    )
    if linked_ctdmr_mode == "extend":
        region_start, region_end, _ = _extend_region_to_first_informative_ctdmr(
            region_start=region_start,
            region_end=region_end,
            linked_ctdmr_callouts=linked_ctdmr_candidates,
        )
        linked_ctdmr_callouts = pd.DataFrame()
    elif linked_ctdmr_mode == "distal":
        linked_ctdmr_callouts = linked_ctdmr_candidates
    else:
        linked_ctdmr_callouts = pd.DataFrame()

    bam_handle = _get_thread_bam_handle(resolved["bam_path"])
    fasta_handle = _get_thread_fasta_handle(resolved["reference_path"])

    shown_reads, all_reads = _fetch_reads(
        resolved["bam_path"],
        sv["chrom"],
        region_start,
        region_end,
        supporting_reads=set(sv["supporting_reads"]),
        max_reads=max_reads,
        bam_handle=bam_handle,
        support_haplotype_only=support_haplotype_only,
    )
    applied_support_haplotype = shown_reads.attrs.get("applied_support_haplotype")
    large_indels = _fetch_large_indels(
        resolved["bam_path"],
        sv["chrom"],
        region_start,
        region_end,
        read_names=set(shown_reads["read_name"].astype(str)) if not shown_reads.empty else set(),
        min_indel_bp=indel_min_bp,
        bam_handle=bam_handle,
    )
    dmrs = _read_ctdmrs(resolved["bed_path"], sv["chrom"], region_start, region_end)
    export_tables = bool(getattr(args, "export_tables", False))
    if export_tables:
        overlap_summary = _summarize_ctdmr_overlap(dmrs, all_reads, int(sv["start"]), int(sv["end"]))
    else:
        overlap_summary = pd.DataFrame()
    if skip_methylation_overlay:
        methyl_df = pd.DataFrame(
            columns=[
                "sv_id",
                "read_name",
                "assignment_status",
                "is_assigned",
                "assigned_celltypes",
                "dmr_assigned_celltypes",
                "dmr_was_assigned",
                "chr",
                "start",
                "end",
                "label",
                "best_group",
                "best_group_leaves",
                "best_dir",
                "mean_methylation",
                "n_cpg_observed",
                "n_cpg_in_dmr",
            ]
        )
    else:
        methyl_df = _compute_supporting_read_ctdmr_methylation(
            sv_id=str(sv["id"]),
            bam_path=resolved["bam_path"],
            reference_path=resolved["reference_path"],
            dmrs=dmrs,
            support_assignment_df=support_assignment_df,
            decoded_assignment_df=decoded_assignment_df,
            logger=logger,
            bam_handle=bam_handle,
            fasta_handle=fasta_handle,
        )

    _plot_sv_panel(
        sv=sv,
        shown_reads=shown_reads,
        all_reads=all_reads,
        large_indels=large_indels,
        dmrs=dmrs,
        linked_ctdmr_callouts=linked_ctdmr_callouts,
        support_assignment_df=support_assignment_df,
        methyl_df=methyl_df,
        region_start=region_start,
        region_end=region_end,
        window=window,
        indel_min_bp=indel_min_bp,
        output_path=output_path,
        dpi=dpi,
        applied_support_haplotype=applied_support_haplotype,
    )

    logger.debug("Wrote SV visualization: %s", output_path)
    if export_tables:
        summary_path = output_path.with_name(f"{output_path.stem}.summary.tsv")
        assign_path = output_path.with_name(f"{output_path.stem}.supporting_reads_assignment.tsv")
        methyl_path = output_path.with_name(f"{output_path.stem}.supporting_reads_ctdmr_methylation.tsv")
        if not overlap_summary.empty:
            overlap_summary = overlap_summary.copy()
            overlap_summary["outside_window_callout"] = False
        if not linked_ctdmr_callouts.empty:
            callout_summary = linked_ctdmr_callouts[
                [
                    "chr",
                    "start",
                    "end",
                    "label",
                    "best_group",
                    "best_group_leaves",
                    "best_dir",
                    "callout_support_count",
                    "callout_support_reads",
                    "callout_assigned_celltypes",
                    "callout_side",
                    "callout_distance_bp",
                ]
            ].copy()
            callout_summary = callout_summary.rename(columns={"callout_support_count": "supporting_read_overlap_count"})
            callout_summary["non_supporting_read_overlap_count"] = 0
            callout_summary["read_overlap_count"] = callout_summary["supporting_read_overlap_count"]
            callout_summary["overlaps_sv_core"] = False
            callout_summary["sv_core_overlap_bp"] = 0
            callout_summary["outside_window_callout"] = True
            if overlap_summary.empty:
                overlap_summary = callout_summary
            else:
                overlap_summary = pd.concat([overlap_summary, callout_summary], ignore_index=True, sort=False)
        overlap_summary.to_csv(summary_path, sep="\t", index=False)
        support_assignment_df.to_csv(assign_path, sep="\t", index=False)
        methyl_df.to_csv(methyl_path, sep="\t", index=False)
        logger.info("Wrote ctDMR overlap summary: %s", summary_path)
        logger.info("Wrote supporting-read assignment table: %s", assign_path)
        logger.info("Wrote supporting-read ctDMR methylation table: %s", methyl_path)
