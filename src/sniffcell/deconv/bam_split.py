from __future__ import annotations

import concurrent.futures
from contextlib import ExitStack
import os
import re
from pathlib import Path

import pandas as pd
import pysam
from tqdm import tqdm

from sniffcell.anno.variant_assignment import _build_group_leaf_sets, _split_pipe_values

_MAX_THREADS_PER_SPLIT_GROUP = 4


def _sanitize_group_label(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    if not text:
        text = "unlabeled"
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("._")
    return text or "unlabeled"


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
    from sniffcell.deconv.deconv import _assignment_reset_index, _prepare_read_assignment_df

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


def _run_group_bam_split(
    bam_path: str,
    names_file_path: str,
    bam_out_path: str,
    threads_per_group: int,
) -> None:
    """Worker: extract reads listed in names_file from bam_path → bam_out_path, then index.

    Must be a top-level function so ProcessPoolExecutor can pickle it.
    Uses samtools view -N (C-level, multi-threaded) for maximum throughput.
    """
    view_args = ["-N", names_file_path, "-b", "-o", bam_out_path]
    if threads_per_group > 1:
        view_args += ["--threads", str(threads_per_group - 1)]
    view_args.append(bam_path)
    pysam.view(*view_args, catch_stdout=False)

    index_args = [bam_out_path]
    if threads_per_group > 1:
        index_args = ["-@", str(threads_per_group - 1), bam_out_path]
    pysam.index(*index_args)


def _index_group_bam(bam_out_path: str, threads: int) -> None:
    index_args = [bam_out_path]
    if threads > 1:
        index_args = ["-@", str(threads - 1), bam_out_path]
    pysam.index(*index_args)


def _compute_split_parallelism(threads: int, n_groups: int) -> tuple[int, int]:
    """Choose a conservative worker layout for requested BAM splits."""
    if n_groups <= 0:
        raise ValueError("n_groups must be positive")

    total_threads = max(1, int(threads))
    max_workers = min(n_groups, total_threads)
    threads_per_group = max(1, total_threads // max_workers)
    threads_per_group = min(threads_per_group, _MAX_THREADS_PER_SPLIT_GROUP)
    return max_workers, threads_per_group


def _split_bam_single_pass(
    bam_path: str,
    planned_specs: list[dict[str, object]],
    bam_out_paths: list[Path],
    threads: int = 1,
) -> None:
    """Split a BAM in one sequential read pass and write all matching outputs.

    This avoids rescanning the full input BAM once per requested group, which is
    often the dominant cost on large files and networked storage.
    """
    read_to_group_indexes: dict[str, list[int]] = {}
    for idx, spec in enumerate(planned_specs):
        for readname in spec["readnames"]:
            read_to_group_indexes.setdefault(str(readname), []).append(idx)

    with pysam.AlignmentFile(bam_path, "rb") as bam_in:
        with ExitStack() as stack:
            bam_out_handles = [
                stack.enter_context(pysam.AlignmentFile(str(path), "wb", template=bam_in))
                for path in bam_out_paths
            ]
            for record in tqdm(bam_in.fetch(until_eof=True), desc="Splitting BAM (single pass)"):
                target_indexes = read_to_group_indexes.get(record.query_name)
                if not target_indexes:
                    continue
                for idx in target_indexes:
                    bam_out_handles[idx].write(record)

    max_index_workers = min(len(bam_out_paths), max(1, threads))
    index_threads = max(1, threads // max_index_workers)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_index_workers) as pool:
        futures = [
            pool.submit(_index_group_bam, str(path), index_threads)
            for path in bam_out_paths
        ]
        for fut in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="Indexing split BAMs",
        ):
            fut.result()


def _split_bam_parallel(
    bam_path: str,
    planned_specs: list[dict[str, object]],
    bam_out_paths: list[Path],
    names_dir: Path,
    threads: int = 1,
) -> None:
    """Split a BAM into N groups in parallel using one samtools process per group.

    Each group runs ``samtools view -N <names_file>`` concurrently. For large
    BAMs, giving each worker a smaller thread allotment is usually faster and
    more stable than saturating the full thread budget per group.
    """
    n_groups = len(planned_specs)
    max_workers, threads_per_group = _compute_split_parallelism(threads, n_groups)

    # Write readnames files to disk (one per group) before spawning workers.
    names_files: list[str] = []
    for spec, bam_out in zip(planned_specs, bam_out_paths):
        stub = spec["file_stub"]
        names_path = names_dir / f"{stub}.readnames.txt"
        with open(str(names_path), "w") as fh:
            for rname in spec["readnames"]:
                fh.write(str(rname) + "\n")
        names_files.append(str(names_path))

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(
                _run_group_bam_split,
                bam_path,
                names_files[i],
                str(bam_out_paths[i]),
                threads_per_group,
            ): i
            for i in range(n_groups)
        }
        for fut in tqdm(
            concurrent.futures.as_completed(futures),
            total=n_groups,
            desc="Splitting BAM (parallel)",
        ):
            fut.result()  # re-raise any worker exception

    # Clean up names files.
    for nf in names_files:
        try:
            os.remove(nf)
        except OSError:
            pass


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

    # Write per-group TSV summaries and collect output paths.
    bam_out_paths: list[Path] = []
    manifest_rows: list[dict[str, object]] = []
    for spec in planned_specs:
        bam_out = split_dir / f"{spec['file_stub']}.bam"
        tsv_out = split_dir / f"{spec['file_stub']}.read_summary.tsv"
        spec["summary_df"].to_csv(tsv_out, sep="\t", index=False)
        bam_out_paths.append(bam_out)
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

    if planned_specs:
        _split_bam_parallel(bam_path, planned_specs, bam_out_paths, split_dir, threads=threads)

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(split_dir / "requested_group_splits.tsv", sep="\t", index=False)
    return manifest_df
