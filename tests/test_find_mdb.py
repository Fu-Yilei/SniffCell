from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from mdb.schema import TrackKey
from mdb.storage import create_cohort_store, create_view_store
from sniffcell.find.find import find_main


def _toy_mdb(path: Path) -> None:
    chroms = ["chr1"]
    offsets = np.asarray([0], dtype=np.int64)
    create_cohort_store(
        str(path),
        chroms=chroms,
        chrom_offsets=offsets,
        pos0=np.asarray([100, 200, 400], dtype=np.uint32),
        backend="zarr",
        block_size=2,
        zarr_row_chunk=2,
    )
    np.savez(
        path / "groups.npz",
        method=np.asarray("sniffcell_loyfer", dtype=object),
        chroms=np.asarray(chroms, dtype=object),
        chrom_offsets=offsets,
        reference_start=np.asarray([100, 200, 400], dtype=np.uint32),
        reference_end=np.asarray([150, 250, 450], dtype=np.uint32),
        source_row_start=np.asarray([0, 2, 4], dtype=np.int64),
        source_row_end=np.asarray([2, 4, 6], dtype=np.int64),
    )
    writer = create_view_store(
        str(path),
        key=TrackKey("5hmC", "combined", "combined"),
        chroms=chroms,
        chrom_offsets=offsets,
        n_rows=3,
        sample_ids=["n1", "o1"],
        bundle_paths=["n.smdb", "o.smdb"],
        platforms=["ont", "ont"],
        source_paths=["n.bed.gz", "o.bed.gz"],
        input_tags=["combined", "combined"],
        block_size=2,
    )
    try:
        writer.write_dense_chrom_block(
            "chr1",
            0,
            3,
            0,
            2,
            np.asarray([[8001, 2001], [7501, 2501], [5501, 4501]], dtype=np.uint16),
        )
    finally:
        writer.close()


def _copy_toy_view(path: Path, assay: str) -> None:
    writer = create_view_store(
        str(path),
        key=TrackKey(assay, "combined", "combined"),
        chroms=["chr1"],
        chrom_offsets=np.asarray([0], dtype=np.int64),
        n_rows=3,
        sample_ids=["n1", "o1"],
        bundle_paths=["n.smdb", "o.smdb"],
        platforms=["ont", "ont"],
        source_paths=["n.bed.gz", "o.bed.gz"],
        input_tags=["combined", "combined"],
        block_size=2,
    )
    try:
        writer.write_dense_chrom_block(
            "chr1",
            0,
            3,
            0,
            2,
            np.asarray([[8001, 2001], [7501, 2501], [5501, 4501]], dtype=np.uint16),
        )
    finally:
        writer.close()


def test_find_reads_mdb_view_and_writes_modification_columns(tmp_path: Path):
    atlas = tmp_path / "atlas.mmdb"
    _toy_mdb(atlas)
    mapping = tmp_path / "celltypes.json"
    mapping.write_text(json.dumps({"brain": {"Neuron": ["n1"], "Oligodendrocyte": ["o1"]}}))
    output = tmp_path / "ctdmr.tsv"

    result = find_main(
        SimpleNamespace(
            mdb=str(atlas),
            assay="5hmC",
            haplotype="combined",
            strand="combined",
            mdb_batch_rows=2,
            paired_stability_filter=False,
            paired_metadata=None,
            npy=None,
            index=None,
            meta=None,
            celltypes_file=str(mapping),
            celltypes_keys="brain",
            output=str(output),
            diff_threshold=0.4,
            min_rows=2,
            min_cpgs=3,
            max_gap_bp=100,
        )
    )

    assert len(result) == 1
    row = result.iloc[0]
    assert row["start"] == 100
    assert row["end"] == 250
    assert row["modification"] == "5hmC"
    assert row["atlas_source"] == str(atlas.resolve())
    written = pd.read_csv(output, sep="\t")
    assert written["modification"].tolist() == ["5hmC"]


def test_find_dual_reads_both_channels_without_collapsing(tmp_path: Path):
    atlas = tmp_path / "atlas.mmdb"
    _toy_mdb(atlas)
    _copy_toy_view(atlas, "5mC")
    mapping = tmp_path / "celltypes.json"
    mapping.write_text(json.dumps({"brain": {"Neuron": ["n1"], "Oligodendrocyte": ["o1"]}}))
    output = tmp_path / "dual.tsv"

    result = find_main(
        SimpleNamespace(
            mdb=str(atlas),
            assay="dual",
            haplotype="combined",
            strand="combined",
            mdb_batch_rows=2,
            paired_stability_filter=False,
            paired_metadata=None,
            npy=None,
            index=None,
            meta=None,
            celltypes_file=str(mapping),
            celltypes_keys="brain",
            output=str(output),
            diff_threshold=0.4,
            min_rows=2,
            min_cpgs=3,
            max_gap_bp=100,
        )
    )

    assert len(result) == 2
    assert result["modification"].tolist() == ["5hmC", "5mC"]
    assert set(zip(result["chr"], result["start"], result["end"])) == {("chr1", 100, 250)}


def test_find_mdb_paired_stability_filter(tmp_path: Path):
    atlas = tmp_path / "paired.mmdb"
    chroms = ["chr1"]
    offsets = np.asarray([0], dtype=np.int64)
    create_cohort_store(
        str(atlas),
        chroms=chroms,
        chrom_offsets=offsets,
        pos0=np.asarray([100, 200], dtype=np.uint32),
        backend="zarr",
        block_size=8,
        zarr_row_chunk=2,
    )
    np.savez(
        atlas / "groups.npz",
        method=np.asarray("sniffcell_loyfer", dtype=object),
        chroms=np.asarray(chroms, dtype=object),
        chrom_offsets=offsets,
        reference_start=np.asarray([100, 200], dtype=np.uint32),
        reference_end=np.asarray([150, 250], dtype=np.uint32),
        source_row_start=np.asarray([0, 2], dtype=np.int64),
        source_row_end=np.asarray([2, 4], dtype=np.int64),
    )
    sample_ids = [f"d{i}_{cell}" for i in range(1, 5) for cell in ("n", "o")]
    writer = create_view_store(
        str(atlas),
        key=TrackKey("5hmC", "combined", "combined"),
        chroms=chroms,
        chrom_offsets=offsets,
        n_rows=2,
        sample_ids=sample_ids,
        bundle_paths=["x"] * 8,
        platforms=["ont"] * 8,
        source_paths=["x"] * 8,
        input_tags=["combined"] * 8,
        block_size=8,
    )
    # Row 1 passes all paired criteria. Row 2 has one donor with only a 0.10 effect.
    values = np.asarray(
        [
            [0.80, 0.20, 0.82, 0.22, 0.78, 0.18, 0.81, 0.21],
            [0.80, 0.20, 0.82, 0.22, 0.78, 0.18, 0.55, 0.45],
        ],
        dtype=np.float32,
    )
    try:
        writer.write_dense_chrom_block("chr1", 0, 2, 0, 8, (values * 10000 + 1).astype(np.uint16))
    finally:
        writer.close()

    mapping = tmp_path / "paired.json"
    mapping.write_text(
        json.dumps(
            {
                "brain": {
                    "Neuron": [f"d{i}_n" for i in range(1, 5)],
                    "Oligodendrocyte": [f"d{i}_o" for i in range(1, 5)],
                }
            }
        )
    )
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(
        [
            {"id": f"d{i}_{cell}", "donor": f"d{i}", "cell_type": "Neuron" if cell == "n" else "Oligodendrocyte"}
            for i in range(1, 5)
            for cell in ("n", "o")
        ]
    ).to_csv(metadata, sep="\t", index=False)

    result = find_main(
        SimpleNamespace(
            mdb=str(atlas),
            assay="5hmC",
            haplotype="combined",
            strand="combined",
            mdb_batch_rows=2,
            paired_stability_filter=True,
            paired_metadata=str(metadata),
            paired_min_donors=4,
            paired_support_effect=0.30,
            paired_min_support=3,
            paired_min_effect=0.15,
            npy=None,
            index=None,
            meta=None,
            celltypes_file=str(mapping),
            celltypes_keys="brain",
            output=str(tmp_path / "paired.tsv"),
            diff_threshold=0.4,
            min_rows=1,
            min_cpgs=2,
            max_gap_bp=100,
        )
    )
    assert len(result) == 1
    assert result.iloc[0]["start"] == 100
    assert result.iloc[0]["paired_donors"] == 4
    assert bool(result.iloc[0]["paired_stability_filter"])
