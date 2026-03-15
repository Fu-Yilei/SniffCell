from __future__ import annotations

import json
import logging
import multiprocessing as mp
import os
import re
from pathlib import Path

import pandas as pd
import pysam
from tqdm import tqdm

from sniffcell.anno.anno import _one_dmr
from sniffcell.anno.variant_assignment import (
    _build_group_leaf_sets,
    _decode_linked_celltypes_from_row,
    _normalize_binary_code,
    _resolve_hierarchy_labels,
    _split_pipe_values,
    _summarize_celltype_links,
)


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


def _build_read_summary(read_assignment_df: pd.DataFrame, *, _prepared: bool = False) -> pd.DataFrame:
    cols = [
        "readname",
        "n_ctdmrs",
        "majority_code",
        "majority_code_count",
        "majority_pct",
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

    # --- vectorized aggregation (replaces Python loop + nested groupby) ---
    ar = assignment_reset
    # Normalise code_token to plain str for groupby keys
    ar = ar.assign(code_token_str=ar["code_token"].astype(str))
    # Capture original readname order (first appearance) so we can restore it at the end.
    readname_order = ar["readname"].drop_duplicates().tolist()

    # 1. Count per (readname, code_token) — single C-level groupby over all rows
    token_counts = (
        ar.groupby(["readname", "code_token_str"], sort=False, dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["count", "code_token_str"], ascending=[False, True], kind="stable")
    )

    # 2. Majority code = first row per readname after descending-count sort
    majority_df = (
        token_counts.drop_duplicates(subset=["readname"], keep="first")
        [["readname", "code_token_str", "count"]]
        .rename(columns={"code_token_str": "majority_code", "count": "majority_code_count"})
    )

    # 3. n_ctdmrs per readname
    n_ctdmrs_df = ar.groupby("readname", sort=False).size().reset_index(name="n_ctdmrs")

    # 4. code_counts string per readname
    code_counts_s = (
        token_counts
        .assign(pair=token_counts["code_token_str"] + ":" + token_counts["count"].astype(str))
        .groupby("readname", sort=False)["pair"]
        .agg(";".join)
        .reset_index(name="code_counts")
    )

    # 5. Merge base stats
    result = (
        majority_df
        .merge(n_ctdmrs_df, on="readname", how="left")
        .merge(code_counts_s, on="readname", how="left")
    )
    result["majority_pct"] = result["majority_code_count"] / result["n_ctdmrs"].replace(0, 1)

    # 6. Decode linked celltypes.
    #    _decode_linked_celltypes_from_row depends only on a small set of columns.
    #    For typical runs there are very few unique combinations (often just 2-10),
    #    so we deduplicate first, decode once per unique key, then join — turning
    #    7M iterrows() calls into O(n_unique_keys) calls + one vectorized merge.
    decode_key_cols = [c for c in ("code_order", "code", "code_token_str", "best_group",
                                    "best_group_leaves", "is_best_group") if c in ar.columns]

    rep_rows = (
        ar
        .merge(
            majority_df[["readname", "majority_code"]],
            left_on=["readname", "code_token_str"],
            right_on=["readname", "majority_code"],
            how="inner",
        )
        .drop_duplicates(subset=["readname"])
    )

    # Compute decode only for unique key combinations (typically ~2-10 rows).
    unique_keys = rep_rows[decode_key_cols].drop_duplicates()
    key_decode_rows: list[dict[str, object]] = []
    for _, urow in unique_keys.iterrows():
        # _decode_linked_celltypes_from_row reads "code_token" but we stored it as "code_token_str"
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

    # Vectorized join: each of the 7M reads maps to one of the ~few decoded keys.
    linked_df = (
        rep_rows[["readname"] + decode_key_cols]
        .merge(decode_map_df, on=decode_key_cols, how="left")
        [["readname", "primary_celltype", "linked_celltypes", "linked_celltype_count",
          "linked_leaf_celltypes", "linked_leaf_celltype_count", "is_multi_celltype_link"]]
    )

    result = result.merge(linked_df, on="readname", how="left")
    for c in ("primary_celltype", "linked_celltypes", "linked_leaf_celltypes", "code_counts"):
        if c in result.columns:
            result[c] = result[c].fillna("")
    # Restore original readname insertion order (matches original groupby sort=False behaviour).
    result = result.set_index("readname").loc[readname_order].reset_index()
    return result[cols]


def _build_deconv_summary(read_assignment_df: pd.DataFrame, *, _prepared: bool = False) -> pd.DataFrame:
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

    rows = []
    for summary_mode, unique_reads in (("all_rows", False), ("per_read", True)):
        summary = _summarize_celltype_links(
            evidence,
            "sample_id",
            unique_reads_for_overlap=unique_reads,
        )
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
                "n_evidence_units": n_unique_reads if unique_reads else n_assignment_rows,
                "n_unique_reads": n_unique_reads,
                "primary_celltype": row.get("primary_celltype", ""),
                "linked_celltypes": row.get("linked_celltypes", ""),
                "linked_celltype_counts": row.get("linked_celltype_counts", ""),
                "linked_celltype_fractions": row.get("linked_celltype_fractions", ""),
                "is_multi_celltype_link": row.get("is_multi_celltype_link", pd.NA),
            }
        )

    return pd.DataFrame(rows, columns=cols)


def _sanitize_group_label(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if not text:
        text = "unlabeled"
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("._")
    return text or "unlabeled"


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


def _normalize_split_label(value: object) -> str:
    text = "" if value is None or pd.isna(value) else str(value).strip().lower()
    return re.sub(r"[^a-z0-9]+", "", text)


def _parse_requested_split_groups(spec_text: str | None) -> list[dict[str, object]]:
    if spec_text is None:
        return []
    text = str(spec_text).strip()
    if not text:
        return []

    specs: list[dict[str, object]] = []
    used_names: set[str] = set()
    for idx, raw_group in enumerate(text.split(";"), start=1):
        raw_group = raw_group.strip()
        if not raw_group:
            continue

        explicit_name = None
        members_text = raw_group
        if "=" in raw_group:
            lhs, rhs = raw_group.split("=", 1)
            lhs = lhs.strip()
            rhs = rhs.strip()
            if not lhs or not rhs:
                raise ValueError(f"Invalid split group definition: {raw_group!r}")
            explicit_name = lhs
            members_text = rhs

        members = [token.strip() for token in members_text.split(",") if token.strip()]
        if not members:
            raise ValueError(f"Split group {raw_group!r} does not contain any labels")

        base_name = explicit_name or "_".join(members)
        file_stub = _sanitize_group_label(base_name)
        if not file_stub:
            file_stub = f"group_{idx}"
        deduped = file_stub
        suffix = 2
        while deduped in used_names:
            deduped = f"{file_stub}.{suffix}"
            suffix += 1
        used_names.add(deduped)

        specs.append(
            {
                "order": idx,
                "raw_group": raw_group,
                "name": explicit_name or ",".join(members),
                "members": members,
                "file_stub": deduped,
            }
        )

    if not specs:
        raise ValueError("split_bam_groups did not contain any usable group definitions")
    return specs


def _build_label_leaf_map(read_assignment_df: pd.DataFrame, *, _prepared: bool = False) -> dict[str, set[str]]:
    if read_assignment_df.empty:
        return {}

    assignment = read_assignment_df.copy() if _prepared else _prepare_read_assignment_df(read_assignment_df)
    assignment_reset = _assignment_reset_index(assignment, read_col="readname")
    label_leaf_map = _build_group_leaf_sets(assignment_reset)

    for col in ("code_order", "best_group_leaves", "other_group_leaves"):
        if col not in assignment_reset.columns:
            continue
        for value in assignment_reset[col].dropna():
            for label in _split_pipe_values(value):
                label_leaf_map.setdefault(label, {label})

    normalized: dict[str, set[str]] = {}
    for label, leaves in label_leaf_map.items():
        norm_label = _normalize_split_label(label)
        if norm_label:
            normalized.setdefault(norm_label, set()).update(str(leaf) for leaf in leaves if str(leaf).strip())
        for leaf in leaves:
            leaf_text = str(leaf).strip()
            norm_leaf = _normalize_split_label(leaf_text)
            if norm_leaf:
                normalized.setdefault(norm_leaf, set()).add(leaf_text)
    return normalized


def _resolve_requested_split_group_targets(
    split_specs: list[dict[str, object]],
    read_assignment_df: pd.DataFrame,
    *,
    _prepared: bool = False,
) -> list[dict[str, object]]:
    label_leaf_map = _build_label_leaf_map(read_assignment_df, _prepared=_prepared)
    resolved_specs: list[dict[str, object]] = []

    for spec in split_specs:
        target_leaves: set[str] = set()
        unmatched_members: list[str] = []
        matched_members: list[str] = []

        for member in spec["members"]:
            norm_member = _normalize_split_label(member)
            leaves = label_leaf_map.get(norm_member)
            if not leaves:
                unmatched_members.append(str(member))
                continue
            target_leaves.update(leaves)
            matched_members.append(str(member))

        if not target_leaves:
            raise ValueError(
                f"Requested split group {spec['raw_group']!r} did not match any ctDMR labels or leaf cell types"
            )

        updated = dict(spec)
        updated["target_leaves"] = sorted(target_leaves)
        updated["matched_members"] = matched_members
        updated["unmatched_members"] = unmatched_members
        resolved_specs.append(updated)

    return resolved_specs


def _plan_requested_split_group_outputs(
    read_summary_df: pd.DataFrame,
    split_specs: list[dict[str, object]],
) -> list[dict[str, object]]:
    if "readname" not in read_summary_df.columns:
        raise ValueError("read_summary_df must contain a 'readname' column")

    summary = read_summary_df.copy()
    if "linked_leaf_celltypes" not in summary.columns:
        summary["linked_leaf_celltypes"] = summary.get("linked_celltypes", "")
    summary["readname"] = summary["readname"].astype(str)
    summary["_linked_leaf_set"] = summary["linked_leaf_celltypes"].map(lambda x: set(_split_pipe_values(x)))

    planned: list[dict[str, object]] = []
    for spec in split_specs:
        target_leaves = set(str(x) for x in spec["target_leaves"])
        subset = summary.loc[
            summary["_linked_leaf_set"].map(lambda leaf_set: bool(leaf_set.intersection(target_leaves)))
        ].copy()
        subset["requested_group"] = spec["name"]
        subset["requested_group_members"] = ",".join(str(x) for x in spec["members"])
        subset["requested_group_target_leaves"] = "|".join(sorted(target_leaves))
        subset["matched_requested_leaves"] = subset["_linked_leaf_set"].map(
            lambda leaf_set: "|".join(sorted(leaf_set.intersection(target_leaves)))
        )
        subset.drop(columns=["_linked_leaf_set"], inplace=True)

        updated = dict(spec)
        updated["summary_df"] = subset
        updated["readnames"] = set(subset["readname"].astype(str))
        planned.append(updated)

    return planned


def _write_requested_split_group_outputs(
    *,
    bam_path: str,
    read_summary_df: pd.DataFrame,
    read_assignment_df: pd.DataFrame,
    split_group_spec: str,
    output_dir: str,
    threads: int = 1,
    _prepared: bool = False,
) -> pd.DataFrame:
    split_specs = _parse_requested_split_groups(split_group_spec)
    resolved_specs = _resolve_requested_split_group_targets(split_specs, read_assignment_df, _prepared=_prepared)
    planned_specs = _plan_requested_split_group_outputs(read_summary_df, resolved_specs)

    split_dir = Path(output_dir)
    split_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    tmp_dir = Path(os.path.expanduser("~/tmp"))
    tmp_dir.mkdir(parents=True, exist_ok=True)

    for spec in tqdm(planned_specs, desc="Splitting BAM by group"):
        bam_out = split_dir / f"{spec['file_stub']}.bam"
        tsv_out = split_dir / f"{spec['file_stub']}.read_summary.tsv"
        spec["summary_df"].to_csv(tsv_out, sep="\t", index=False)

        names_file = tmp_dir / f"sniffcell_readnames_{spec['file_stub']}.txt"
        names_file.write_text("\n".join(str(r) for r in spec["readnames"]))
        try:
            view_args = ["-N", str(names_file), "-b", "-o", str(bam_out)]
            if threads > 1:
                view_args += ["--threads", str(threads - 1)]
            view_args.append(str(bam_path))
            pysam.view(*view_args, catch_stdout=False)
            index_args = ["-@", str(threads - 1), str(bam_out)] if threads > 1 else [str(bam_out)]
            pysam.index(*index_args)
        finally:
            try:
                names_file.unlink()
            except Exception:
                pass

        manifest_rows.append(
            {
                "requested_group": spec["name"],
                "requested_group_members": ",".join(str(x) for x in spec["members"]),
                "matched_members": ",".join(str(x) for x in spec["matched_members"]),
                "unmatched_members": ",".join(str(x) for x in spec["unmatched_members"]),
                "target_leaves": "|".join(str(x) for x in spec["target_leaves"]),
                "n_reads": int(len(spec["readnames"])),
                "bam_path": str(bam_out),
                "read_summary_path": str(tsv_out),
            }
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(split_dir / "requested_group_splits.tsv", sep="\t", index=False)
    return manifest_df


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

    if threads < 1:
        raise ValueError("threads must be >= 1")
    if read_assignment_mode not in {"closest_reference_mean", "kmeans"}:
        raise ValueError("read_assignment_mode must be one of: closest_reference_mean, kmeans")

    outputs = _resolve_output_paths(output_arg)
    os.makedirs(outputs["output_dir"], exist_ok=True)

    logger.info(
        "Starting deconvolution: bed=%s bam=%s ref=%s threads=%d read_assignment_mode=%s out=%s",
        bed_input,
        input_bam,
        reference,
        threads,
        read_assignment_mode,
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
        outputs=outputs,
    )
    logger.info("Wrote deconv run manifest: %s", manifest_path)

    resume = bool(getattr(args, "resume", False))
    reads_tsv = outputs["reads"]
    read_summary_tsv = outputs["read_summary"]

    if resume and os.path.exists(reads_tsv) and os.path.exists(read_summary_tsv):
        logger.info("--resume: loading existing reads classification from %s", reads_tsv)
        read_assign_df = pd.read_csv(reads_tsv, sep="\t", low_memory=False)
        logger.info(
            "Preparing read assignment dataframe (n_rows=%d, n_unique_reads=%d)...",
            len(read_assign_df),
            int(read_assign_df["readname"].nunique()) if "readname" in read_assign_df.columns else 0,
        )
        prepared_df = _prepare_read_assignment_df(read_assign_df)
        logger.info("--resume: loading existing read summary from %s", read_summary_tsv)
        read_summary_df = pd.read_csv(read_summary_tsv, sep="\t", low_memory=False)
    else:
        if resume:
            logger.warning("--resume requested but existing TSVs not found; running full ctDMR phase")

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
        read_summary_df = _build_read_summary(prepared_df, _prepared=True)
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

    logger.info("Building overall deconv summary...")
    summary_df = _build_deconv_summary(prepared_df, _prepared=True)
    summary_df.to_csv(outputs["summary"], sep="\t", index=False)

    logger.info("Wrote row-level read assignments: %s", outputs["reads"])
    logger.info("Wrote ctDMR block summary: %s", outputs["blocks"])
    logger.info("Wrote per-read deconvolution summary: %s", outputs["read_summary"])
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
