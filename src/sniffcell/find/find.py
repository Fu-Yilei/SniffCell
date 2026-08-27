import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

import numpy as np
import pandas as pd

from sniffcell.find import ctdmr
from sniffcell.find.ctdmr import means_from_mapping


def _dedupe_keep_order(values: Iterable[str]) -> list[str]:
    out: list[str] = []
    seen = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _normalize_group_mapping(mapping: Mapping[str, Any]) -> Dict[str, list[str]]:
    out: Dict[str, list[str]] = {}
    for group, samples in mapping.items():
        if not isinstance(samples, (list, tuple)):
            continue
        cleaned = [str(x).strip() for x in samples if str(x).strip()]
        if cleaned:
            out[str(group)] = _dedupe_keep_order(cleaned)
    return out


def _extract_mapping_keys(atlas_payload: Mapping[str, Any]) -> list[str]:
    keys: list[str] = []
    for key, value in atlas_payload.items():
        if str(key).startswith("__"):
            continue
        if isinstance(value, Mapping):
            keys.append(str(key))
    return sorted(keys)


def _resolve_celltype_mapping(celltypes_key: str, atlas_payload: Mapping[str, Any]) -> Dict[str, list[str]]:
    if celltypes_key not in atlas_payload:
        raise KeyError(
            f"Cell type key '{celltypes_key}' not found. "
            f"Available keys: {_extract_mapping_keys(atlas_payload)}"
        )
    raw = atlas_payload[celltypes_key]
    if not isinstance(raw, Mapping):
        raise TypeError(f"Cell type key '{celltypes_key}' must map to an object of group->sample_list.")
    mapping = _normalize_group_mapping(raw)
    if len(mapping) < 2:
        raise ValueError(f"Cell type key '{celltypes_key}' must contain at least 2 non-empty groups.")
    return mapping


def _build_empty_output(group_order: list[str]) -> pd.DataFrame:
    cols = [
        "chr",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "n_rows",
        "n_cpgs",
        "bp_len",
        "best_group",
        "other_group",
        "best_dir",
        "mean_margin",
        "second_best_margin",
        "rest_std_mean",
        "mean_best_value",
        "mean_rest_value",
        "best_group_leaves",
        "other_group_leaves",
        "hyper_group_leaves",
        "hypo_group_leaves",
        "code_order",
    ] + [f"mean_{ct}" for ct in group_order]
    return pd.DataFrame(columns=cols)


def _resolve_mdb_assays(value: str) -> list[str]:
    """Resolve a single assay, comma-separated assays, or the dual shorthand."""
    raw = str(value).strip()
    if not raw:
        raise ValueError("--assay cannot be empty")
    if raw.lower() == "dual":
        return ["5mC", "5hmC"]
    assays = _dedupe_keep_order(part.strip() for part in raw.split(",") if part.strip())
    if not assays:
        raise ValueError("--assay must name at least one MDB assay")
    return assays


def _load_mdb_index(mdb_path: str) -> pd.DataFrame:
    groups_path = os.path.join(mdb_path, "groups.npz")
    if not os.path.isfile(groups_path):
        raise ValueError(
            f"MDB atlas must contain groups.npz with the SniffCell Loyfer row definitions: {mdb_path}"
        )
    groups = np.load(groups_path, allow_pickle=True)
    required = {
        "chroms",
        "chrom_offsets",
        "reference_start",
        "reference_end",
        "source_row_start",
        "source_row_end",
    }
    missing = required - set(groups.files)
    if missing:
        raise ValueError(f"MDB groups.npz is missing: {', '.join(sorted(missing))}")

    chroms = [str(x) for x in groups["chroms"]]
    offsets = np.asarray(groups["chrom_offsets"], dtype=np.int64)
    starts = np.asarray(groups["reference_start"], dtype=np.int64)
    ends = np.asarray(groups["reference_end"], dtype=np.int64)
    source_starts = np.asarray(groups["source_row_start"], dtype=np.int64)
    source_ends = np.asarray(groups["source_row_end"], dtype=np.int64)
    n_rows = len(starts)
    if not (len(ends) == len(source_starts) == len(source_ends) == n_rows):
        raise ValueError("MDB groups.npz row arrays have inconsistent lengths")
    chrom_ends = np.asarray(
        [int(offsets[i + 1]) if i + 1 < len(offsets) else n_rows for i in range(len(chroms))],
        dtype=np.int64,
    )
    codes = np.empty(n_rows, dtype=np.int16)
    for code, (lo, hi) in enumerate(zip(offsets, chrom_ends, strict=False)):
        codes[int(lo) : int(hi)] = code
    return pd.DataFrame(
        {
            "chr": pd.Categorical.from_codes(codes, categories=chroms, ordered=True),
            "start": starts,
            "end": ends,
            "startCpG": source_starts,
            "endCpG": source_ends,
        }
    )


def _means_from_mdb_reader(
    reader,
    sample_ids: list[str],
    mapping: Mapping[str, list[str]],
    batch_rows: int,
) -> Dict[str, pd.Series]:
    sample_to_index = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    group_indices: Dict[str, list[int]] = {}
    for group, declared_samples in mapping.items():
        indices = [sample_to_index[sample_id] for sample_id in declared_samples if sample_id in sample_to_index]
        if indices:
            group_indices[group] = indices
    if len(group_indices) < 2:
        raise ValueError(
            "Need at least two cell-type groups with samples present in the selected MDB view. "
            f"Available samples: {len(sample_ids)}"
        )

    means = {
        group: np.full(reader.n_rows, np.nan, dtype=np.float32)
        for group in group_indices
    }
    batch_rows = max(int(batch_rows), 1)
    for start in range(0, reader.n_rows, batch_rows):
        stop = min(start + batch_rows, reader.n_rows)
        block = reader.get_block(slice(start, stop))
        for group, indices in group_indices.items():
            selected = block[:, indices]
            if len(indices) == 1:
                means[group][start:stop] = selected[:, 0]
                continue
            finite = np.isfinite(selected)
            count = finite.sum(axis=1)
            total = np.nansum(selected, axis=1, dtype=np.float64)
            out = np.full(stop - start, np.nan, dtype=np.float32)
            np.divide(total, count, out=out, where=count > 0)
            means[group][start:stop] = out
    return {group: pd.Series(values, copy=False) for group, values in means.items()}


def _paired_stability_mask(
    reader,
    sample_ids: list[str],
    mapping: Mapping[str, list[str]],
    *,
    metadata_path: str,
    batch_rows: int,
    median_effect: float,
    support_effect: float,
    min_support: int,
    min_effect: float,
    min_donors: int,
) -> tuple[np.ndarray, int]:
    if len(mapping) != 2:
        raise ValueError("Paired stability filtering currently requires exactly two declared cell-type groups")
    metadata = pd.read_csv(metadata_path, sep=None, engine="python", dtype=str)
    id_col = "sample_id" if "sample_id" in metadata.columns else "id" if "id" in metadata.columns else None
    required = {"donor", "cell_type"}
    if id_col is None or not required.issubset(metadata.columns):
        raise ValueError("Paired metadata must contain sample_id or id, donor, and cell_type columns")
    metadata = metadata.set_index(id_col, drop=False)
    sample_to_index = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    groups = list(mapping)
    present = {
        group: [sample_id for sample_id in mapping[group] if sample_id in sample_to_index]
        for group in groups
    }
    donors = sorted(
        {
            str(metadata.loc[sample_id, "donor"])
            for group in groups
            for sample_id in present[group]
            if sample_id in metadata.index
        }
    )
    pairs: list[tuple[int, int]] = []
    for donor in donors:
        hits = []
        for group in groups:
            group_hits = [
                sample_id
                for sample_id in present[group]
                if sample_id in metadata.index and str(metadata.loc[sample_id, "donor"]) == donor
            ]
            if len(group_hits) > 1:
                raise ValueError(f"Donor {donor} has multiple samples in paired group {group}")
            hits.append(group_hits)
        if all(len(group_hits) == 1 for group_hits in hits):
            pairs.append((sample_to_index[hits[0][0]], sample_to_index[hits[1][0]]))
    if len(pairs) < int(min_donors):
        raise ValueError(f"Paired stability filter requires {min_donors} complete donors; found {len(pairs)}")
    if min_support > len(pairs):
        raise ValueError(f"paired_min_support={min_support} exceeds the {len(pairs)} complete donor pairs")

    mask = np.zeros(reader.n_rows, dtype=bool)
    batch_rows = max(int(batch_rows), 1)
    for start in range(0, reader.n_rows, batch_rows):
        stop = min(start + batch_rows, reader.n_rows)
        block = reader.get_block(slice(start, stop))
        deltas = np.column_stack([block[:, left] - block[:, right] for left, right in pairs])
        finite = np.isfinite(deltas).all(axis=1)
        positive = np.all(deltas > 0, axis=1)
        negative = np.all(deltas < 0, axis=1)
        safe_deltas = np.where(finite[:, None], deltas, 0.0)
        abs_deltas = np.abs(safe_deltas)
        median_abs = np.abs(np.median(safe_deltas, axis=1))
        support_count = np.sum(abs_deltas >= float(support_effect), axis=1)
        minimum_abs = np.min(abs_deltas, axis=1)
        mask[start:stop] = (
            finite
            & (positive | negative)
            & (median_abs >= float(median_effect))
            & (support_count >= int(min_support))
            & (minimum_abs >= float(min_effect))
        )
    return mask, len(pairs)


def _load_mdb_atlas(
    mdb_path: str,
    *,
    assay: str,
    haplotype: str,
    strand: str,
    mapping: Mapping[str, list[str]],
    batch_rows: int,
    paired_options: dict | None = None,
):
    try:
        from mdb.schema import TrackKey
        from mdb.storage import load_view_reader
    except ImportError as exc:
        raise RuntimeError(
            "MDB-backed sniffcell find requires the mdb package. Install mdb or use the legacy --npy input."
        ) from exc

    key = TrackKey(assay=assay, haplotype=haplotype, strand=strand)
    reader, columns, _ = load_view_reader(mdb_path, key)
    try:
        cpg_index = _load_mdb_index(mdb_path)
        if reader.n_rows != len(cpg_index):
            raise ValueError(
                f"MDB view rows {reader.n_rows:,} do not match groups.npz rows {len(cpg_index):,}"
            )
        sample_ids = list(columns["sample_id"])
        mean_by_group = _means_from_mdb_reader(reader, sample_ids, mapping, batch_rows)
        if paired_options is None:
            row_eligible_mask = None
            paired_donors = 0
        else:
            row_eligible_mask, paired_donors = _paired_stability_mask(
                reader,
                sample_ids,
                mapping,
                batch_rows=batch_rows,
                **paired_options,
            )
    finally:
        reader.close()
    return cpg_index, mean_by_group, sample_ids, key, row_eligible_mask, paired_donors


def _write_igv_bed(annotated_dmrs: pd.DataFrame, output_path: Path) -> None:
    igv_out = output_path.with_suffix(output_path.suffix + ".igv.bed")
    if annotated_dmrs.empty:
        igv_out.write_text("", encoding="utf-8")
        return

    bed9 = annotated_dmrs[["chr", "start", "end", "name", "score", "strand"]].copy()
    bed9["thickStart"] = bed9["start"]
    bed9["thickEnd"] = bed9["end"]
    bed9["itemRgb"] = 0
    bed9.to_csv(igv_out, sep="\t", index=False, header=False)


def find_main(args):
    logger = logging.getLogger(__name__)
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    celltypes_file = args.celltypes_file
    celltypes_key = args.celltypes_keys
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Loading cell type definitions from %s", celltypes_file)
    with open(celltypes_file, "r", encoding="utf-8") as f:
        atlas_mapping = json.load(f)

    mapping = _resolve_celltype_mapping(celltypes_key, atlas_mapping)
    logger.info("Using key '%s' with %d declared groups", celltypes_key, len(mapping))
    using_mdb = bool(getattr(args, "mdb", None))
    paired_filter = bool(getattr(args, "paired_stability_filter", False))
    if paired_filter and not using_mdb:
        raise ValueError("--paired-stability-filter requires --mdb")
    if paired_filter and not getattr(args, "paired_metadata", None):
        raise ValueError("--paired-stability-filter requires --paired-metadata")
    if using_mdb:
        assays = _resolve_mdb_assays(args.assay)
        assay_results: list[pd.DataFrame] = []
        for assay in assays:
            logger.info(
                "Loading MDB atlas %s view=%s__%s__%s",
                args.mdb,
                assay,
                args.haplotype,
                args.strand,
            )
            paired_options = None
            if paired_filter:
                paired_options = {
                    "metadata_path": args.paired_metadata,
                    "median_effect": args.diff_threshold,
                    "support_effect": args.paired_support_effect,
                    "min_support": args.paired_min_support,
                    "min_effect": args.paired_min_effect,
                    "min_donors": args.paired_min_donors,
                }
            cpg_index, mean_by_group, all_celltypes, mdb_key, row_eligible_mask, paired_donors = _load_mdb_atlas(
                args.mdb,
                assay=assay,
                haplotype=args.haplotype,
                strand=args.strand,
                mapping=mapping,
                batch_rows=args.mdb_batch_rows,
                paired_options=paired_options,
            )
            logger.info("MDB rows=%d samples=%d", len(cpg_index), len(all_celltypes))
            if row_eligible_mask is not None:
                logger.info(
                    "Paired stability filter retained %d/%d rows across %d donors",
                    int(row_eligible_mask.sum()),
                    len(row_eligible_mask),
                    paired_donors,
                )

            group_order = list(mean_by_group.keys())
            logger.info("Usable groups with available samples: %s", "|".join(group_order))
            if len(group_order) < 2:
                raise ValueError(
                    f"Need at least 2 usable groups after matching metadata columns for key '{celltypes_key}'."
                )
            dmrs = ctdmr.call_ct_combination_dmrs(
                idx_df=cpg_index,
                mean_by_group=mean_by_group,
                diff_threshold=args.diff_threshold,
                min_rows=args.min_rows,
                min_cpgs=args.min_cpgs,
                min_bp=0,
                direction="both",
                max_gap_bp=args.max_gap_bp,
                row_eligible_mask=row_eligible_mask,
                bed_out=None,
            )
            if dmrs.empty:
                assay_dmrs = _build_empty_output(group_order)
            else:
                assay_dmrs = dmrs.copy()
            assay_dmrs["modification"] = mdb_key.assay
            assay_dmrs["atlas_source"] = str(Path(args.mdb).resolve())
            assay_dmrs["paired_stability_filter"] = paired_filter
            assay_dmrs["paired_donors"] = paired_donors
            if paired_filter:
                assay_dmrs["paired_median_effect_min"] = float(args.diff_threshold)
                assay_dmrs["paired_support_effect"] = float(args.paired_support_effect)
                assay_dmrs["paired_min_support"] = int(args.paired_min_support)
                assay_dmrs["paired_min_effect"] = float(args.paired_min_effect)
            assay_results.append(assay_dmrs)
            logger.info("Assay %s produced %d ctDMRs", assay, len(assay_dmrs))

        annotated_dmrs = pd.concat(assay_results, ignore_index=True, sort=False)
    else:
        npy_file = args.npy
        index_file = args.index
        meta_file = args.meta
        logger.warning("Using legacy NPY atlas input; use --mdb for new atlas construction and dual assays")
        logger.info("Loading legacy atlas matrix from %s", npy_file)
        all_celltype_blocks = np.load(npy_file)
        logger.info("Atlas matrix shape: %s", getattr(all_celltype_blocks, "shape", None))

        logger.info("Loading sample names from %s", meta_file)
        with open(meta_file, "r", encoding="utf-8") as meta_f:
            all_celltypes = [line.strip() for line in meta_f if line.strip()]
        logger.info("Sample name count: %d", len(all_celltypes))

        logger.info("Loading CpG index from %s", index_file)
        cpg_index = pd.read_csv(
            index_file,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "startCpG", "endCpG"],
        )
        logger.info("CpG rows: %d", len(cpg_index))

        M_df = pd.DataFrame(all_celltype_blocks, columns=all_celltypes, index=cpg_index.index)
        logger.info("Methylation matrix dataframe shape: %s", M_df.shape)
        mean_by_group = means_from_mapping(M_df, mapping)
        group_order = list(mean_by_group.keys())
        logger.info("Usable groups with available samples: %s", "|".join(group_order))
        if len(group_order) < 2:
            raise ValueError(
                f"Need at least 2 usable groups after matching metadata columns for key '{celltypes_key}'."
            )
        annotated_dmrs = ctdmr.call_ct_combination_dmrs(
            idx_df=cpg_index,
            mean_by_group=mean_by_group,
            diff_threshold=args.diff_threshold,
            min_rows=args.min_rows,
            min_cpgs=args.min_cpgs,
            min_bp=0,
            direction="both",
            max_gap_bp=args.max_gap_bp,
            row_eligible_mask=None,
            bed_out=None,
        )
        if annotated_dmrs.empty:
            annotated_dmrs = _build_empty_output(group_order)

    sort_cols = [
        c for c in ["chr", "start", "end", "modification", "best_group"]
        if c in annotated_dmrs.columns
    ]
    if sort_cols:
        annotated_dmrs = annotated_dmrs.sort_values(sort_cols, kind="stable", ignore_index=True)

    annotated_dmrs.to_csv(output_path, sep="\t", index=False)
    _write_igv_bed(annotated_dmrs, output_path)

    logger.info("Wrote annotation-ready ctDMR BED/TSV: %s", output_path)
    logger.info(
        "Wrote IGV BED9 companion file: %s",
        output_path.with_suffix(output_path.suffix + ".igv.bed"),
    )
    logger.info("find_main completed (key=%s, total_dmrs=%d)", celltypes_key, len(annotated_dmrs))
    return annotated_dmrs
