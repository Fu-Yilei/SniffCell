from __future__ import annotations

import json
import logging
import multiprocessing as mp
import os
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from sniffcell.anno.anno import _one_dmr
from sniffcell.anno.variant_assignment import (
    _build_group_leaf_sets,
    _compute_per_read_consensus_df,
    _decode_linked_celltypes_from_row,
    _normalize_binary_code,
    _resolve_hierarchy_labels,
    _split_code_token_schema_bits,
    _summarize_celltype_links,
)
from sniffcell.deconv.bam_split import _sanitize_group_label, _write_requested_split_group_outputs


def _resolve_output_paths(output_arg: str) -> dict[str, str]:
    output_arg = str(output_arg)
    path = Path(output_arg)

    if output_arg.endswith(os.sep) or output_arg.endswith("/"):
        output_dir = path
        summary_path = output_dir / "deconv_summary.tsv"
    elif path.exists() and path.is_dir():
        output_dir = path
        summary_path = output_dir / "deconv_summary.tsv"
    elif path.suffix:
        output_dir = path.parent if str(path.parent) else Path(".")
        summary_path = path
    else:
        output_dir = path
        summary_path = output_dir / "deconv_summary.tsv"

    return {
        "output_dir": str(output_dir),
        "summary": str(summary_path),
        "reads": str(output_dir / "deconv_reads_classification.tsv"),
        "blocks": str(output_dir / "deconv_blocks_classification.tsv"),
        "read_summary": str(output_dir / "deconv_read_summary.tsv"),
        "group_dir": str(output_dir / "deconv_reads_by_group"),
        "requested_split_dir": str(output_dir / "deconv_requested_group_splits"),
        "manifest": str(output_dir / "deconv_run_manifest.json"),
    }


def _write_deconv_run_manifest(
    *,
    output_path: str,
    bam: str,
    reference: str,
    bed: str,
    threads: int,
    read_assignment_mode: str,
    split_bam_groups: str | None,
    per_read_min_agreement: float,
    outputs: dict[str, str],
) -> str:
    manifest_path = outputs["manifest"]
    payload = {
        "command": "deconv",
        "version": "v1",
        "inputs": {
            "bam": os.path.abspath(bam),
            "reference": os.path.abspath(reference),
            "bed": os.path.abspath(bed),
        },
        "runtime": {
            "threads": int(threads),
            "read_assignment_mode": str(read_assignment_mode),
            "split_bam_groups": split_bam_groups,
            "per_read_min_agreement": float(per_read_min_agreement),
        },
        "outputs": {
            "requested_output": os.path.abspath(output_path),
            "output_dir": os.path.abspath(outputs["output_dir"]),
            "summary": os.path.abspath(outputs["summary"]),
            "reads_classification": os.path.abspath(outputs["reads"]),
            "blocks_classification": os.path.abspath(outputs["blocks"]),
            "read_summary": os.path.abspath(outputs["read_summary"]),
            "group_split_dir": os.path.abspath(outputs["group_dir"]),
            "requested_group_split_dir": os.path.abspath(outputs["requested_split_dir"]),
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return manifest_path


def _load_ctdmr_bed(bed_path: str) -> pd.DataFrame:
    bed = pd.read_csv(bed_path, sep="\t")
    if not bed.empty and isinstance(bed.columns[0], str) and bed.columns[0].startswith("#"):
        bed.rename(columns={bed.columns[0]: bed.columns[0].lstrip("#")}, inplace=True)
    bed = bed.drop_duplicates(ignore_index=True)
    bed = bed.sort_values(["chr", "start"], ignore_index=True)

    required = ["chr", "start", "end", "best_group", "best_dir"]
    missing = [col for col in required if col not in bed.columns]
    if missing:
        raise ValueError(f"BED missing required columns: {missing}")
    return bed


def _write_empty_reads_table(path: str) -> None:
    empty_reads = pd.DataFrame(
        columns=[
            "chr",
            "start",
            "end",
            "cpgstart",
            "cpgend",
            "best_group",
            "other_group",
            "is_best_group",
            "code_order",
            "best_group_leaves",
            "other_group_leaves",
            "hyper_group_leaves",
            "hypo_group_leaves",
            "code",
        ]
    )
    empty_reads.index.name = "readname"
    empty_reads.to_csv(path, sep="\t", index=True, header=True)


def _write_empty_blocks_table(path: str) -> None:
    pd.DataFrame(columns=["chr", "start", "end", "cpgstart", "cpgend"]).to_csv(
        path,
        sep="\t",
        index=False,
        header=True,
    )


_STREAM_WRITE_BUFFER = 200   # flush to disk every N ctDMR results
_STREAM_CONCAT_BATCH = 1000  # merge staging read DFs every N items to bound peak memory


def _stream_ctdmr_classification(
    *,
    bed_df: pd.DataFrame,
    input_bam: str,
    reference: str,
    threads: int,
    read_assignment_mode: str,
    reads_out: str,
    blocks_out: str,
) -> pd.DataFrame:
    """Stream ctDMR classification, write TSV outputs, and return the full
    reads DataFrame in memory to avoid a costly re-read of the 2+ GB file."""
    tasks = [(dict(row), input_bam, reference, read_assignment_mode) for _, row in bed_df.iterrows()]
    open(reads_out, "w").close()
    open(blocks_out, "w").close()

    reads_header_written = False
    blocks_header_written = False
    blocks_cols_locked: list[str] | None = None

    # Batched disk-write buffers (reduces fsync/append syscalls by ~100×)
    read_write_buf: list[pd.DataFrame] = []
    block_write_buf: list[pd.DataFrame] = []
    buf_count = 0

    # In-memory accumulation — avoids re-reading the 2+ GB TSV downstream.
    # Staged in two levels to keep concat overhead low.
    read_staging: list[pd.DataFrame] = []
    read_chunks: list[pd.DataFrame] = []

    def _flush_to_disk() -> None:
        nonlocal reads_header_written, blocks_header_written, blocks_cols_locked
        if read_write_buf:
            combined = pd.concat(read_write_buf, axis=0)
            combined.to_csv(
                reads_out,
                sep="\t",
                index=True,
                mode="a",
                header=not reads_header_written,
            )
            reads_header_written = True
            read_write_buf.clear()
        if block_write_buf:
            combined = pd.concat(block_write_buf, axis=0)
            if not blocks_header_written:
                blocks_cols_locked = list(combined.columns)
                combined.to_csv(blocks_out, sep="\t", index=False, mode="a", header=True)
                blocks_header_written = True
            else:
                assert blocks_cols_locked is not None
                combined.reindex(columns=blocks_cols_locked).to_csv(
                    blocks_out, sep="\t", index=False, mode="a", header=False,
                )
            block_write_buf.clear()

    with mp.Pool(threads) as pool:
        for res in tqdm(
            pool.imap(_one_dmr, tasks, chunksize=50),
            total=len(tasks),
            desc="Processing ctDMRs",
        ):
            if res is None:
                continue
            read_df, block_df = res

            if read_df is not None and not read_df.empty:
                read_write_buf.append(read_df)
                read_staging.append(read_df)
                if len(read_staging) >= _STREAM_CONCAT_BATCH:
                    read_chunks.append(pd.concat(read_staging, axis=0))
                    read_staging.clear()

            if block_df is not None and not block_df.empty:
                block_write_buf.append(block_df)

            buf_count += 1
            if buf_count >= _STREAM_WRITE_BUFFER:
                _flush_to_disk()
                buf_count = 0

        _flush_to_disk()  # final flush

    if not reads_header_written:
        _write_empty_reads_table(reads_out)
    if not blocks_header_written:
        _write_empty_blocks_table(blocks_out)

    # Final in-memory concat and return
    if read_staging:
        read_chunks.append(pd.concat(read_staging, axis=0))
    if read_chunks:
        result = pd.concat(read_chunks, axis=0)
        result.index.name = "readname"
        return result
    return pd.DataFrame(
        columns=[
            "chr", "start", "end", "cpgstart", "cpgend",
            "best_group", "other_group", "is_best_group",
            "code_order", "best_group_leaves", "other_group_leaves",
            "hyper_group_leaves", "hypo_group_leaves", "code",
        ],
        index=pd.Index([], name="readname"),
    )


def _prepare_read_assignment_df(read_assignment_df: pd.DataFrame) -> pd.DataFrame:
    assignment = read_assignment_df.copy()
    assignment.index = assignment.index.astype(str)

    if "code" not in assignment.columns:
        raise ValueError("read_assignment_df must contain a 'code' column")

    assignment["code"] = assignment["code"].astype("string")
    for col in ("best_group", "code_order", "best_group_leaves"):
        if col not in assignment.columns:
            assignment[col] = pd.Series(pd.NA, index=assignment.index, dtype="string")
        else:
            assignment[col] = assignment[col].astype("string")
    if "is_best_group" not in assignment.columns:
        assignment["is_best_group"] = False
    else:
        assignment["is_best_group"] = assignment["is_best_group"].fillna(False).astype(bool)

    # Vectorized replacement for the row-by-row apply of _normalize_binary_code.
    # The apply over millions of rows was a significant bottleneck.
    code_raw = assignment["code"]
    na_mask = code_raw.isna()
    code_str = code_raw.astype(str).str.strip()
    # Strip ".0" suffix from integer-valued floats (e.g. pandas reads "1" as "1.0")
    dot_zero = code_str.str.endswith(".0") & code_str.str[:-2].str.isdigit()
    code_str = code_str.where(~dot_zero, code_str.str[:-2])
    # Zero-pad binary codes to schema width when code_order is present
    has_schema = (
        assignment["code_order"].notna()
        & assignment["code_order"].astype(str).str.strip().str.len().gt(0)
    )
    if has_schema.any():
        n_labels = (
            assignment["code_order"].where(has_schema)
            .astype(str).str.split("|").str.len()
            .fillna(0).astype(int)
        )
        is_binary = code_str.str.fullmatch(r"[01]+").fillna(False)
        needs_pad = has_schema & is_binary & (code_str.str.len() < n_labels)
        if needs_pad.any():
            pad_lens = n_labels[needs_pad]
            code_str = code_str.copy()
            code_str[needs_pad] = [s.zfill(l) for s, l in zip(code_str[needs_pad], pad_lens)]
    code_str_result = code_str.astype("string")
    code_str_result[na_mask] = pd.NA
    assignment["code"] = code_str_result

    non_empty_orders = assignment["code_order"].dropna()
    non_empty_orders = non_empty_orders[non_empty_orders.str.len() > 0]
    use_qualified_code = non_empty_orders.nunique(dropna=True) > 1
    assignment["code_token"] = assignment["code"]
    if use_qualified_code:
        has_order = assignment["code_order"].str.len().fillna(0) > 0
        assignment.loc[has_order, "code_token"] = (
            assignment.loc[has_order, "code_order"] + "::" + assignment.loc[has_order, "code"]
        )
    return assignment


def _assignment_reset_index(assignment: pd.DataFrame, *, read_col: str) -> pd.DataFrame:
    out = assignment.reset_index()
    first_col = out.columns[0]
    if first_col != read_col:
        out.rename(columns={first_col: read_col}, inplace=True)
    out[read_col] = out[read_col].astype(str)
    return out


def _build_read_summary(
    read_assignment_df: pd.DataFrame,
    *,
    _prepared: bool = False,
    per_read_min_agreement: float = 0.66,
) -> pd.DataFrame:
    cols = [
        "readname",
        "n_ctdmrs",
        "consensus_code",
        "agreement_pct",
        "is_mixed",
        "used_intersection",
        "primary_celltype",
        "linked_celltypes",
        "linked_celltype_count",
        "linked_leaf_celltypes",
        "linked_leaf_celltype_count",
        "is_multi_celltype_link",
        "code_counts",
    ]
    if read_assignment_df.empty:
        return pd.DataFrame(columns=cols)

    assignment = read_assignment_df.copy() if _prepared else _prepare_read_assignment_df(read_assignment_df)
    assignment_reset = _assignment_reset_index(assignment, read_col="readname")
    group_leaf_sets = _build_group_leaf_sets(assignment_reset)

    ar = assignment_reset.assign(code_token_str=assignment_reset["code_token"].astype(str))
    # Capture original readname order so we can restore it at the end.
    readname_order = ar["readname"].drop_duplicates().tolist()

    # 1. code_counts string per readname (from raw ctDMR codes, for auditability)
    token_counts = (
        ar.groupby(["readname", "code_token_str"], sort=False, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["count", "code_token_str"], ascending=[False, True], kind="stable")
    )
    code_counts_s = (
        token_counts
        .assign(pair=token_counts["code_token_str"] + ":" + token_counts["count"].astype(str))
        .groupby("readname", sort=False)["pair"]
        .agg(";".join)
        .reset_index(name="code_counts")
    )

    # 2. Per-read consensus via bitwise intersection (replaces majority vote)
    per_read = _compute_per_read_consensus_df(assignment, per_read_min_agreement=per_read_min_agreement)

    # 3. One representative row per read (first occurrence) for code_order / best_group_leaves lookup
    decode_key_cols = [c for c in ("code_order", "code", "code_token_str", "best_group",
                                    "best_group_leaves", "is_best_group") if c in ar.columns]
    rep_rows = ar.drop_duplicates(subset=["readname"], keep="first")[["readname"] + decode_key_cols].copy()

    # 4. Override code / code_token_str in rep_rows with the consensus code for non-mixed reads
    consensus_map = per_read.set_index("readname")[
        ["read_consensus_code", "read_is_mixed", "read_agreement_pct",
         "read_n_ctdmrs", "read_used_intersection"]
    ]
    rep_rows = rep_rows.merge(consensus_map, left_on="readname", right_index=True, how="left")
    non_mixed_mask = ~rep_rows["read_is_mixed"].fillna(False)
    consensus_codes = rep_rows["read_consensus_code"].where(rep_rows["read_consensus_code"].notna(), "")
    rep_rows.loc[non_mixed_mask, "code_token_str"] = consensus_codes[non_mixed_mask]
    rep_rows.loc[non_mixed_mask, "code"] = consensus_codes[non_mixed_mask].map(
        lambda x: _split_code_token_schema_bits(x)[1] if x else pd.NA
    )
    # Mixed reads: clear code fields so decoding returns empty strings
    rep_rows.loc[~non_mixed_mask, "code_token_str"] = ""
    rep_rows.loc[~non_mixed_mask, "code"] = pd.NA

    # 5. Decode linked celltypes — deduplicate unique key combinations (typically ~few rows)
    unique_keys = rep_rows.loc[non_mixed_mask, decode_key_cols].drop_duplicates()
    key_decode_rows: list[dict[str, object]] = []
    for _, urow in unique_keys.iterrows():
        row_as_series = urow.rename({"code_token_str": "code_token"}) if "code_token_str" in urow.index else urow
        _schema, linked = _decode_linked_celltypes_from_row(row_as_series)
        resolved = _resolve_hierarchy_labels(linked, group_leaf_sets) if linked else []
        key_decode_rows.append(
            {
                **{c: urow[c] for c in decode_key_cols},
                "primary_celltype": resolved[0] if resolved else "",
                "linked_celltypes": "|".join(resolved),
                "linked_celltype_count": len(resolved),
                "linked_leaf_celltypes": "|".join(linked),
                "linked_leaf_celltype_count": len(linked),
                "is_multi_celltype_link": len(resolved) > 1 if resolved else pd.NA,
            }
        )

    decode_map_df = pd.DataFrame(key_decode_rows) if key_decode_rows else pd.DataFrame(
        columns=decode_key_cols + ["primary_celltype", "linked_celltypes",
                                    "linked_celltype_count", "linked_leaf_celltypes",
                                    "linked_leaf_celltype_count", "is_multi_celltype_link"]
    )
    linked_df = (
        rep_rows[["readname"] + decode_key_cols]
        .merge(decode_map_df, on=decode_key_cols, how="left")
        [["readname", "primary_celltype", "linked_celltypes", "linked_celltype_count",
          "linked_leaf_celltypes", "linked_leaf_celltype_count", "is_multi_celltype_link"]]
    )

    # 6. Assemble result
    result = (
        per_read
        .rename(columns={
            "readname": "readname",
            "read_consensus_code": "consensus_code",
            "read_agreement_pct": "agreement_pct",
            "read_is_mixed": "is_mixed",
            "read_n_ctdmrs": "n_ctdmrs",
            "read_used_intersection": "used_intersection",
        })
        .merge(code_counts_s, on="readname", how="left")
        .merge(linked_df, on="readname", how="left")
    )
    for c in ("primary_celltype", "linked_celltypes", "linked_leaf_celltypes", "code_counts"):
        if c in result.columns:
            result[c] = result[c].fillna("")
    # Mixed reads: ensure cell type fields are empty strings
    mixed_mask = result["is_mixed"].fillna(False)
    for c in ("primary_celltype", "linked_celltypes", "linked_leaf_celltypes"):
        result.loc[mixed_mask, c] = ""
    result.loc[mixed_mask, "linked_celltype_count"] = 0
    result.loc[mixed_mask, "linked_leaf_celltype_count"] = 0

    result = result.set_index("readname").loc[readname_order].reset_index()
    return result[cols]


def _build_deconv_summary(
    read_assignment_df: pd.DataFrame,
    *,
    _prepared: bool = False,
    per_read_min_agreement: float = 0.66,
) -> pd.DataFrame:
    cols = [
        "summary_mode",
        "n_assignment_rows",
        "n_evidence_units",
        "n_unique_reads",
        "primary_celltype",
        "linked_celltypes",
        "linked_celltype_counts",
        "linked_celltype_fractions",
        "is_multi_celltype_link",
    ]

    n_assignment_rows = int(read_assignment_df.shape[0])
    n_unique_reads = int(read_assignment_df.index.astype(str).nunique()) if not read_assignment_df.empty else 0
    if read_assignment_df.empty:
        return pd.DataFrame(
            [
                {
                    "summary_mode": "all_rows",
                    "n_assignment_rows": 0,
                    "n_evidence_units": 0,
                    "n_unique_reads": 0,
                    "primary_celltype": "",
                    "linked_celltypes": "",
                    "linked_celltype_counts": "",
                    "linked_celltype_fractions": "",
                    "is_multi_celltype_link": pd.NA,
                },
                {
                    "summary_mode": "per_read",
                    "n_assignment_rows": 0,
                    "n_evidence_units": 0,
                    "n_unique_reads": 0,
                    "primary_celltype": "",
                    "linked_celltypes": "",
                    "linked_celltype_counts": "",
                    "linked_celltype_fractions": "",
                    "is_multi_celltype_link": pd.NA,
                },
            ],
            columns=cols,
        )

    assignment = read_assignment_df.copy() if _prepared else _prepare_read_assignment_df(read_assignment_df)
    evidence = _assignment_reset_index(assignment, read_col="read")
    evidence["sample_id"] = "sample"

    # Per-read consensus: used for the per_read summary mode
    per_read_consensus = _compute_per_read_consensus_df(
        assignment, per_read_min_agreement=per_read_min_agreement
    )
    n_mixed_reads = int(per_read_consensus["read_is_mixed"].sum())
    non_mixed = per_read_consensus[~per_read_consensus["read_is_mixed"]]

    # Build synthetic evidence df for the per_read mode: one row per non-mixed read,
    # code_token overridden with the consensus code so _summarize_celltype_links decodes correctly.
    if not non_mixed.empty:
        non_mixed_reps = (
            evidence.drop_duplicates(subset=["read"], keep="first")
            .merge(
                non_mixed[["readname", "read_consensus_code"]],
                left_on="read", right_on="readname", how="inner",
            )
            .copy()
        )
        non_mixed_reps["code_token"] = non_mixed_reps["read_consensus_code"]
        non_mixed_reps["code"] = non_mixed_reps["read_consensus_code"].map(
            lambda x: _split_code_token_schema_bits(x)[1] if pd.notna(x) and str(x) else pd.NA
        )
        per_read_evidence = non_mixed_reps
    else:
        per_read_evidence = pd.DataFrame(columns=evidence.columns)

    rows = []
    for summary_mode in ("all_rows", "per_read"):
        if summary_mode == "all_rows":
            ev = evidence
            n_ev = n_assignment_rows
        else:
            ev = per_read_evidence
            n_ev = n_unique_reads - n_mixed_reads

        summary = _summarize_celltype_links(ev, "sample_id", unique_reads_for_overlap=False)
        if summary.empty:
            rows.append(
                {
                    "summary_mode": summary_mode,
                    "n_assignment_rows": n_assignment_rows,
                    "n_evidence_units": 0,
                    "n_unique_reads": n_unique_reads,
                    "primary_celltype": "",
                    "linked_celltypes": "",
                    "linked_celltype_counts": "",
                    "linked_celltype_fractions": "",
                    "is_multi_celltype_link": pd.NA,
                }
            )
            continue

        row = summary.iloc[0]
        rows.append(
            {
                "summary_mode": summary_mode,
                "n_assignment_rows": n_assignment_rows,
                "n_evidence_units": n_ev,
                "n_unique_reads": n_unique_reads,
                "primary_celltype": row.get("primary_celltype", ""),
                "linked_celltypes": row.get("linked_celltypes", ""),
                "linked_celltype_counts": row.get("linked_celltype_counts", ""),
                "linked_celltype_fractions": row.get("linked_celltype_fractions", ""),
                "is_multi_celltype_link": row.get("is_multi_celltype_link", pd.NA),
            }
        )

    return pd.DataFrame(rows, columns=cols)


def _write_group_split_reads(read_assignment_df: pd.DataFrame, output_dir: str) -> list[str]:
    split_dir = Path(output_dir)
    split_dir.mkdir(parents=True, exist_ok=True)
    if read_assignment_df.empty:
        return []

    assignment = read_assignment_df.copy()
    assignment_reset = _assignment_reset_index(assignment, read_col="readname")

    written_paths: list[str] = []
    used_names: dict[str, str] = {}

    groups = list(assignment_reset.groupby("best_group", sort=True, dropna=False))
    for group_value, group_df in tqdm(groups, desc="Writing group read tables"):
        group_label = group_value
        if pd.isna(group_label) or not str(group_label).strip():
            fallback = group_df["best_group_leaves"].dropna().astype(str)
            group_label = fallback.iloc[0] if not fallback.empty else "unlabeled"

        base_name = _sanitize_group_label(group_label)
        filename = base_name
        suffix = 2
        while filename in used_names and used_names[filename] != str(group_label):
            filename = f"{base_name}.{suffix}"
            suffix += 1
        used_names[filename] = str(group_label)

        subset = group_df.set_index("readname")
        subset.index.name = "readname"
        out_path = split_dir / f"{filename}.tsv"
        subset.to_csv(out_path, sep="\t", index=True)
        written_paths.append(str(out_path))

    return written_paths


def deconv_main(args):
    logger = logging.getLogger("sniffcell.deconv")

    bed_input = str(args.bed)
    input_bam = str(args.input)
    reference = str(args.reference)
    output_arg = str(args.output)
    threads = int(getattr(args, "threads", 1))
    read_assignment_mode = str(getattr(args, "read_assignment_mode", "closest_reference_mean")).strip().lower()
    split_bam_groups = getattr(args, "split_bam_groups", None)
    split_bam_groups = None if split_bam_groups is None else str(split_bam_groups).strip() or None
    per_read_min_agreement = float(getattr(args, "per_read_min_agreement", 0.66))
    skip_overall_summary = bool(getattr(args, "skip_overall_summary", False))

    if threads < 1:
        raise ValueError("threads must be >= 1")
    if read_assignment_mode not in {"closest_reference_mean", "kmeans"}:
        raise ValueError("read_assignment_mode must be one of: closest_reference_mean, kmeans")
    if not (0.0 <= per_read_min_agreement <= 1.0):
        raise ValueError("per_read_min_agreement must be in [0, 1]")

    outputs = _resolve_output_paths(output_arg)
    os.makedirs(outputs["output_dir"], exist_ok=True)

    logger.info(
        "Starting deconvolution: bed=%s bam=%s ref=%s threads=%d read_assignment_mode=%s "
        "per_read_min_agreement=%.3f out=%s",
        bed_input,
        input_bam,
        reference,
        threads,
        read_assignment_mode,
        per_read_min_agreement,
        output_arg,
    )

    manifest_path = _write_deconv_run_manifest(
        output_path=output_arg,
        bam=input_bam,
        reference=reference,
        bed=bed_input,
        threads=threads,
        read_assignment_mode=read_assignment_mode,
        split_bam_groups=split_bam_groups,
        per_read_min_agreement=per_read_min_agreement,
        outputs=outputs,
    )
    logger.info("Wrote deconv run manifest: %s", manifest_path)

    resume = bool(getattr(args, "resume", False))
    reads_tsv = outputs["reads"]
    read_summary_tsv = outputs["read_summary"]

    if resume and os.path.exists(reads_tsv):
        logger.info("--resume: loading existing reads classification from %s", reads_tsv)
        read_assign_df = pd.read_csv(reads_tsv, sep="\t", low_memory=False, index_col=0)
        logger.info(
            "Preparing read assignment dataframe (n_rows=%d, n_unique_reads=%d)...",
            len(read_assign_df),
            int(read_assign_df.index.astype(str).nunique()),
        )
        prepared_df = _prepare_read_assignment_df(read_assign_df)
        # Always recompute the read summary so that parameter changes (e.g. per_read_min_agreement)
        # take effect without re-running the slow ctDMR scan phase.
        logger.info("Building per-read summary (recomputed from cached classification)...")
        read_summary_df = _build_read_summary(
            prepared_df, _prepared=True, per_read_min_agreement=per_read_min_agreement
        )
        logger.info("Writing per-read summary: %s", read_summary_tsv)
        read_summary_df.to_csv(read_summary_tsv, sep="\t", index=False)
    else:
        if resume:
            logger.warning("--resume requested but %s not found; running full ctDMR phase", reads_tsv)

        bed_df = _load_ctdmr_bed(bed_input)
        logger.info("Loaded BED with %d unique ctDMRs", len(bed_df))

        # _stream_ctdmr_classification writes the TSV outputs AND returns the full
        # reads DataFrame in memory, eliminating the costly 2+ GB re-read.
        read_assign_df = _stream_ctdmr_classification(
            bed_df=bed_df,
            input_bam=input_bam,
            reference=reference,
            threads=threads,
            read_assignment_mode=read_assignment_mode,
            reads_out=outputs["reads"],
            blocks_out=outputs["blocks"],
        )

        logger.info(
            "Preparing read assignment dataframe (n_rows=%d, n_unique_reads=%d)...",
            len(read_assign_df),
            int(read_assign_df.index.astype(str).nunique()),
        )
        prepared_df = _prepare_read_assignment_df(read_assign_df)

        logger.info("Building per-read summary...")
        read_summary_df = _build_read_summary(
            prepared_df, _prepared=True, per_read_min_agreement=per_read_min_agreement
        )
        logger.info("Writing per-read summary: %s", outputs["read_summary"])
        read_summary_df.to_csv(outputs["read_summary"], sep="\t", index=False)

    logger.info("Writing per-group read tables: %s", outputs["group_dir"])
    group_paths = _write_group_split_reads(prepared_df, outputs["group_dir"])

    requested_split_manifest = None
    if split_bam_groups:
        logger.info("Splitting BAM by requested groups: %s", outputs["requested_split_dir"])
        requested_split_manifest = _write_requested_split_group_outputs(
            bam_path=input_bam,
            read_summary_df=read_summary_df,
            read_assignment_df=prepared_df,
            split_group_spec=split_bam_groups,
            output_dir=outputs["requested_split_dir"],
            threads=threads,
            _prepared=True,
        )

    if skip_overall_summary:
        logger.info("Skipping overall deconv summary by request: %s", outputs["summary"])
    else:
        logger.info("Building overall deconv summary...")
        summary_df = _build_deconv_summary(
            prepared_df, _prepared=True, per_read_min_agreement=per_read_min_agreement
        )
        summary_df.to_csv(outputs["summary"], sep="\t", index=False)

    logger.info("Wrote row-level read assignments: %s", outputs["reads"])
    logger.info("Wrote ctDMR block summary: %s", outputs["blocks"])
    logger.info("Wrote per-read deconvolution summary: %s", outputs["read_summary"])
    if skip_overall_summary:
        logger.info("Skipped overall deconvolution summary: %s", outputs["summary"])
    else:
        logger.info("Wrote overall deconvolution summary: %s", outputs["summary"])
    logger.info("Wrote %d per-group split read tables under: %s", len(group_paths), outputs["group_dir"])
    if requested_split_manifest is not None:
        unmatched = requested_split_manifest["unmatched_members"].fillna("").astype(str)
        for _, row in requested_split_manifest.loc[unmatched.ne("")].iterrows():
            logger.warning(
                "Requested split group %s had unmatched labels: %s",
                row["requested_group"],
                row["unmatched_members"],
            )
        logger.info(
            "Wrote %d requested BAM/TSV group splits under: %s",
            len(requested_split_manifest),
            outputs["requested_split_dir"],
        )
