"""Behavioral coverage of paired display evidence and legacy dispatch."""

import logging
from unittest.mock import patch

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from sniffcell.viz import assay, viz


def catalog():
    return pd.DataFrame(
        [
            dict(
                chr="chr1",
                chr_norm="1",
                start=10,
                end=20,
                label="Neuron",
                best_group="Neuron",
                best_group_leaves="Neuron",
                best_dir="hyper",
                name=mod,
                modification=mod,
                mean_Neuron=m,
                mean_Oligodendrocyte=1 - m,
            )
            for mod, m in [("5mC", 0.3), ("5hmC", 0.6)]
        ]
    )


def test_loader_and_decode_preserve_same_coordinate_assays(tmp_path):
    path = tmp_path / "atlas.tsv"
    catalog().to_csv(path, sep="\t", index=False)
    rows = viz._load_ctdmr_table_cached(str(path))
    assert list(rows.modification) == ["5mC", "5hmC"]
    rows["read_name"] = "r"
    rows["code_order"] = "Neuron|Oligodendrocyte"
    rows["code"] = "10"
    detail = viz._decode_read_assignment_rows(rows)
    assert list(detail.modification) == ["5mC", "5hmC"]
    assert len(rows.drop_duplicates(assay.marker_columns(rows))) == 2


def test_paired_means_use_same_observed_cpgs_and_missing_is_not_zero():
    # Extra m-only CpG must not bias the paired m:h summary.
    matrices = {
        "5mC": pd.DataFrame(
            [[0.2, 0.9], [0.8, 0.7]], index=["r", "missing_h"], columns=[11, 15]
        ),
        "5hmC": pd.DataFrame([[0.6, np.nan]], index=["r"], columns=[11, 15]),
    }
    meta = pd.DataFrame(
        [
            dict(
                read_name=r,
                assignment_status="assigned",
                is_assigned=True,
                assigned_celltypes="Neuron",
            )
            for r in ["r", "missing_h"]
        ]
    )
    with patch.object(
        assay,
        "methyl_matrix_from_bam",
        side_effect=lambda *a, **k: (matrices[k["modification"]], [11, 15]),
    ):
        result = assay.extract_assay_means(
            sv_id="v",
            bam_path="b",
            reference_path="f",
            dmrs=catalog(),
            support_assignment_df=meta,
            decoded_assignment_df=pd.DataFrame(),
            logger=logging.getLogger(),
        )
    observed = result[result.read_name.eq("r")].set_index("modification")
    assert observed.loc["5mC", "mean_methylation"] == 0.2
    assert observed.loc["5hmC", "mean_methylation"] == 0.6
    assert observed.n_cpg_observed.eq(1).all()
    missing = result[result.read_name.eq("missing_h")]
    assert missing.n_cpg_observed.eq(0).all() and missing.mean_methylation.isna().all()


def test_fraction_geometry_height_and_missing_channels():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 100)
    assay._draw_boxes(
        ax, 10, 90, 1, assay.READ_BOX_HEIGHT, {"5mC": 0.2, "5hmC": 0.6}, paired=True
    )
    assert ax.patches[0].get_width() == pytest.approx(20)
    assert ax.patches[1].get_width() == pytest.approx(60)
    assay._draw_boxes(ax, 10, 90, 2, assay.READ_BOX_HEIGHT, {"5mC": 0.2}, paired=False)
    assert all(p.get_height() == assay.READ_BOX_HEIGHT for p in ax.patches)
    n = len(ax.patches)
    assay._draw_boxes(
        ax, 10, 90, 3, assay.READ_BOX_HEIGHT, {"5mC": 0.2, "5hmC": np.nan}, paired=True
    )
    assert len(ax.patches) == n
    assert assay.paired_fractions(0, 0) is None
    plt.close(fig)


def test_partial_overlaps_have_distinct_lanes_but_exact_pair_shares_lane():
    frame = pd.DataFrame({"start": [10, 10, 15, 30], "end": [20, 20, 25, 40]})
    lanes, count = assay.interval_lanes(frame)
    assert count == 2 and len(lanes) == 3
    assert lanes[(10, 20)] != lanes[(15, 25)]
    assert lanes[(30, 40)] == lanes[(10, 20)]


def test_legacy_renderer_and_extractor_dispatch_are_unchanged():
    for frame in [
        catalog().drop(columns="modification"),
        catalog().assign(modification="modifiedC"),
    ]:
        with (
            patch.object(viz, "_plot_legacy_sv_panel", return_value="legacy") as old,
            patch.object(assay, "plot_assay_panel") as new,
        ):
            assert (
                viz._plot_sv_panel(dmrs=frame, linked_ctdmr_callouts=pd.DataFrame())
                == "legacy"
            )
            old.assert_called_once()
            new.assert_not_called()
        with patch.object(
            viz, "_compute_legacy_read_ctdmr_methylation", return_value="legacy"
        ) as old:
            assert (
                viz._compute_supporting_read_ctdmr_methylation(dmrs=frame) == "legacy"
            )
            old.assert_called_once()


def test_fallback_uses_consensus_threshold_and_saved_calls_take_priority(tmp_path):
    rows = pd.DataFrame(
        [
            dict(
                code=code,
                code_order="Neuron|Oligodendrocyte",
                best_group="Neuron",
                best_group_leaves="Neuron",
            )
            for code in ["10", "10", "10", "01", "01"]
        ],
        index=pd.Index(["r"] * 5, name="readname"),
    )
    calls = assay.read_consensus(rows)
    assert calls.iloc[0].is_mixed
    path = tmp_path / "summary.tsv"
    pd.DataFrame(
        [dict(readname="r", is_mixed=False, linked_leaf_celltypes="Oligodendrocyte")]
    ).to_csv(path, sep="\t", index=False)
    calls = assay.read_consensus(rows, str(path))
    meta = pd.DataFrame(
        [
            dict(
                read_name="r",
                is_assigned=True,
                assigned_celltypes="Neuron",
                assignment_status="assigned",
            )
        ]
    )
    result = assay.apply_consensus(meta, calls)
    assert result.iloc[0].assigned_celltypes == "Oligodendrocyte"


def test_genomic_window_counts_distinct_intervals_not_assay_rows():
    frame = catalog()
    assert assay.genomic_window(frame, "chr1", 30, 40, 1) == (10, 40)
    assert assay.genomic_window(frame, "chr2", 30, 40, 1) == (30, 40)
    with pytest.raises(ValueError):
        assay.genomic_window(frame, "chr1", 30, 40, -1)
