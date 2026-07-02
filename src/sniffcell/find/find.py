import json
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping

import numpy as np
import pandas as pd

from sniffcell.find import ctdmr
from sniffcell.find.ctdmr import means_from_mapping
from sniffcell.tissue_atlas import (
    resolve_tissue_key,
    write_tissue_catalog_manifest,
)


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

    npy_file = args.npy
    index_file = args.index
    meta_file = args.meta
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    celltypes_file = args.celltypes_file
    requested_celltypes_key = args.celltypes_keys

    logger.info("Loading atlas matrix from %s", npy_file)
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

    logger.info("Loading cell type definitions from %s", celltypes_file)
    with open(celltypes_file, "r", encoding="utf-8") as f:
        atlas_mapping = json.load(f)

    celltypes_key, tissue_row = resolve_tissue_key(requested_celltypes_key, atlas_mapping)
    if celltypes_key != requested_celltypes_key:
        logger.info("Resolved tissue key '%s' to atlas key '%s'", requested_celltypes_key, celltypes_key)

    mapping = _resolve_celltype_mapping(celltypes_key, atlas_mapping)
    logger.info("Using key '%s' with %d declared groups", celltypes_key, len(mapping))

    mean_by_group = means_from_mapping(M_df, mapping)
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
        bed_out=None,
    )

    if dmrs.empty:
        annotated_dmrs = _build_empty_output(group_order)
    else:
        annotated_dmrs = dmrs.copy()
        sort_cols = [c for c in ["chr", "start", "end", "best_group"] if c in annotated_dmrs.columns]
        if sort_cols:
            annotated_dmrs = annotated_dmrs.sort_values(sort_cols, kind="stable", ignore_index=True)

    annotated_dmrs.to_csv(output_path, sep="\t", index=False)
    _write_igv_bed(annotated_dmrs, output_path)
    if tissue_row is not None:
        write_tissue_catalog_manifest(
            output_path=output_path,
            celltypes_file=celltypes_file,
            requested_key=requested_celltypes_key,
            resolved_key=celltypes_key,
            tissue_row=tissue_row,
        )

    logger.info("Wrote annotation-ready ctDMR BED/TSV: %s", output_path)
    logger.info(
        "Wrote IGV BED9 companion file: %s",
        output_path.with_suffix(output_path.suffix + ".igv.bed"),
    )
    logger.info("find_main completed (key=%s, total_dmrs=%d)", celltypes_key, len(annotated_dmrs))
    return annotated_dmrs
