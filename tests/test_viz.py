import unittest
import tempfile
import json
import logging
from pathlib import Path

import pandas as pd
import pysam

from sniffcell.viz.viz import (
    _assign_ctdmr_label_lanes,
    _build_variant_payload_from_table_row,
    _collect_large_indels_from_cigar,
    _build_methylation_heatmap_matrix,
    _expand_interval_for_visibility,
    _fetch_reads,
    _extend_region_to_first_informative_ctdmr,
    _reference_celltype_mean_columns,
    _resolve_viz_runtime_inputs,
    _summarize_supporting_read_assignments,
)


class TestVizVariantPayload(unittest.TestCase):
    def test_harmonized_tr_payload_normalizes_medaka_read_headers(self):
        payload = _build_variant_payload_from_table_row(
            pd.Series(
                {
                    "chrom": "chr13",
                    "start": 102161539,
                    "end": 102161869,
                    "variant_class": "TR",
                    "variant_id": "chr13_102161539_102161869",
                    "variant_subtype": "expansion_all",
                    "group_a_read_names": json.dumps(
                        [
                            "uuidA_chr13_102161539_102161869_pad_102161289_102162119_fwd_hap2_phased-set1_ploidy2",
                            "uuidB_chr13_102161539_102161869_pad_102161289_102162119_rev_hap2_phased-set1_ploidy2",
                        ]
                    ),
                    "group_b_read_names": "[]",
                }
            )
        )

        self.assertEqual(payload["variant_class"], "TR")
        self.assertEqual(payload["supporting_reads"], {"uuidA", "uuidB"})

    def test_non_tr_payload_preserves_read_names_containing_chr(self):
        payload = _build_variant_payload_from_table_row(
            pd.Series(
                {
                    "chrom": "chr1",
                    "start": 100,
                    "end": 200,
                    "variant_class": "SV",
                    "variant_id": "sv1",
                    "group_a_read_names": json.dumps(["sample_chr_read"]),
                    "group_b_read_names": "[]",
                }
            )
        )

        self.assertEqual(payload["supporting_reads"], {"sample_chr_read"})


class TestVizSupportingReadAssignment(unittest.TestCase):
    def test_supporting_reads_assigned_and_unassigned(self):
        assignment_df = pd.DataFrame(
            {
                "chr": ["1", "1", "1"],
                "start": [100, 120, 5000],
                "end": [140, 160, 5060],
                "code_order": ["A|B", "A|B", "A|B"],
                "code": ["10", "xyz", "01"],
                "best_group": ["A", "A", "B"],
                "best_group_leaves": ["A", "A", "B"],
                "is_best_group": [True, False, True],
            },
            index=pd.Index(["r1", "r2", "rX"], name="readname"),
        )
        assignment_df["read_name"] = assignment_df.index.astype(str)

        summary_df, detail_df = _summarize_supporting_read_assignments(
            assignment_df,
            supporting_reads={"r1", "r2", "r3"},
            sv_chrom="1",
            sv_start=1000,
            sv_end=1010,
            window=1000,
            region_start=0,
            region_end=1000,
            assignment_available=True,
        )

        summary = summary_df.set_index("read_name")
        self.assertTrue(bool(summary.loc["r1", "is_assigned"]))
        self.assertEqual(summary.loc["r1", "assignment_status"], "assigned")
        self.assertEqual(summary.loc["r1", "assigned_celltypes"], "A")
        self.assertFalse(bool(summary.loc["r2", "is_assigned"]))
        self.assertEqual(summary.loc["r2", "assignment_status"], "unassigned_unresolved_code")
        self.assertFalse(bool(summary.loc["r3", "is_assigned"]))
        self.assertEqual(summary.loc["r3", "assignment_status"], "unassigned_no_overlap_rows")

        self.assertEqual(len(detail_df), 1)
        self.assertEqual(detail_df.iloc[0]["read_name"], "r1")
        self.assertEqual(detail_df.iloc[0]["assigned_celltypes"], "A")

    def test_supporting_reads_can_remain_assigned_outside_local_plot_window(self):
        assignment_df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1200],
                "end": [1245],
                "code_order": ["Neuron|Oligodendrocyte"],
                "code": ["10"],
                "best_group": ["Neuron"],
                "best_group_leaves": ["Neuron"],
                "is_best_group": [True],
            },
            index=pd.Index(["r1"], name="readname"),
        )
        assignment_df["read_name"] = assignment_df.index.astype(str)

        summary_df, detail_df = _summarize_supporting_read_assignments(
            assignment_df,
            supporting_reads={"r1"},
            sv_chrom="1",
            sv_start=2000,
            sv_end=2010,
            window=1000,
            region_start=1900,
            region_end=2100,
            assignment_available=True,
            clip_to_region=False,
        )

        self.assertEqual(len(summary_df), 1)
        self.assertTrue(bool(summary_df.iloc[0]["is_assigned"]))
        self.assertEqual(summary_df.iloc[0]["assignment_status"], "assigned")
        self.assertEqual(summary_df.iloc[0]["assigned_celltypes"], "Neuron")
        self.assertEqual(len(detail_df), 1)

    def test_prescreened_assignment_rows_are_not_distance_filtered_again(self):
        assignment_df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [12000],
                "end": [12045],
                "code_order": ["Neuron|Oligodendrocyte"],
                "code": ["10"],
                "best_group": ["Neuron"],
                "best_group_leaves": ["Neuron"],
                "is_best_group": [True],
            },
            index=pd.Index(["r1"], name="readname"),
        )
        assignment_df["read_name"] = assignment_df.index.astype(str)

        summary_df, detail_df = _summarize_supporting_read_assignments(
            assignment_df,
            supporting_reads={"r1"},
            sv_chrom="1",
            sv_start=2000,
            sv_end=2010,
            window=1000,
            region_start=1900,
            region_end=2100,
            assignment_available=True,
            clip_to_region=False,
            prescreened_assignment_rows=True,
        )

        self.assertTrue(bool(summary_df.iloc[0]["is_assigned"]))
        self.assertEqual(summary_df.iloc[0]["assigned_celltypes"], "Neuron")
        self.assertEqual(len(detail_df), 1)


class TestVizHeatmapMatrix(unittest.TestCase):
    def test_heatmap_matrix_limits_and_orders(self):
        methyl_df = pd.DataFrame(
            {
                "read_name": ["r1", "r2", "r3", "r1"],
                "is_assigned": [True, False, True, True],
                "label": ["DMR_A", "DMR_A", "DMR_A", "DMR_B"],
                "chr": ["1", "1", "1", "1"],
                "start": [100, 100, 100, 300],
                "end": [200, 200, 200, 360],
                "mean_methylation": [0.8, 0.2, 0.7, 0.9],
            }
        )

        heat = _build_methylation_heatmap_matrix(methyl_df, max_reads=2, max_dmrs=1)
        self.assertEqual(heat.shape, (2, 1))
        self.assertEqual(list(heat.index), ["r1", "r3"])

    def test_reference_celltype_mean_columns_uses_all_mean_columns(self):
        dmrs = pd.DataFrame(
            {
                "label": ["A", "B"],
                "start": [100, 200],
                "end": [150, 260],
                "mean_T-cell": [0.8, 0.3],
                "mean_B-cell": [0.2, 0.7],
                "mean_best_value": [0.8, 0.7],  # should be excluded
            }
        )
        cols = _reference_celltype_mean_columns(dmrs)
        self.assertEqual(cols, ["mean_T-cell", "mean_B-cell"])


class TestVizManifestResolution(unittest.TestCase):
    def test_resolve_inputs_from_anno_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            manifest = {
                "inputs": {
                    "bam": "/tmp/in.bam",
                    "vcf": "/tmp/in.vcf.gz",
                    "reference": "/tmp/ref.fa",
                    "bed": "/tmp/dmrs.tsv",
                },
                "runtime": {"window": 12000},
                "outputs": {"reads_classification": str(anno_dir / "reads_classification.tsv")},
            }
            (anno_dir / "anno_run_manifest.json").write_text(
                json.dumps(manifest),
                encoding="utf-8",
            )

            class Args:
                anno_output = str(anno_dir)
                input = None
                vcf = None
                reference = None
                bed = None
                read_assignment = None
                output = None
                format = "png"
                window = 5000
                sv_id = "sv1"
                kanpig_read_names = None

            resolved = _resolve_viz_runtime_inputs(Args(), logging.getLogger("test"))
            self.assertEqual(resolved["bam_path"], "/tmp/in.bam")
            self.assertEqual(resolved["vcf_path"], "/tmp/in.vcf.gz")
            self.assertEqual(resolved["reference_path"], "/tmp/ref.fa")
            self.assertEqual(resolved["bed_path"], "/tmp/dmrs.tsv")
            self.assertEqual(resolved["window"], 12000)
            self.assertTrue(str(resolved["output_path"]).endswith("sv1.viz.png"))

    def test_resolve_inputs_uses_variant_specific_lite_runtime(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            runtime_path = anno_dir / "lite_variant_runtime.tsv"
            pd.DataFrame(
                [
                    {
                        "id": "sv1",
                        "bam": "/tmp/sv1.bam",
                        "bed": "/tmp/sv1.ctdmr.tsv",
                        "reference": "/tmp/sv1.fa",
                    }
                ]
            ).to_csv(runtime_path, sep="\t", index=False)
            manifest = {
                "inputs": {
                    "bam": "/tmp/default.bam",
                    "vcf": "/tmp/variants.tsv",
                    "reference": "/tmp/default.fa",
                    "bed": "/tmp/default.ctdmr.tsv",
                },
                "outputs": {"lite_variant_runtime": str(runtime_path)},
            }
            (anno_dir / "anno_run_manifest.json").write_text(
                json.dumps(manifest),
                encoding="utf-8",
            )

            class Args:
                anno_output = str(anno_dir)
                input = None
                vcf = None
                reference = None
                bed = None
                read_assignment = None
                output = None
                format = "png"
                window = 5000
                sv_id = "sv1"
                kanpig_read_names = None

            resolved = _resolve_viz_runtime_inputs(Args(), logging.getLogger("test"))

        self.assertEqual(resolved["bam_path"], "/tmp/sv1.bam")
        self.assertEqual(resolved["bed_path"], "/tmp/sv1.ctdmr.tsv")
        self.assertEqual(resolved["reference_path"], "/tmp/sv1.fa")
        self.assertTrue(resolved["prescreened_assignment_rows"])


class TestVizSupportHaplotypeFiltering(unittest.TestCase):
    def _write_test_bam(self, bam_path: Path, rows: list[dict[str, object]]) -> None:
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
            for row in rows:
                seg = pysam.AlignedSegment()
                seg.query_name = str(row["read_name"])
                seg.query_sequence = "A" * 100
                seg.flag = 0
                seg.reference_id = 0
                seg.reference_start = int(row.get("start", 1000))
                seg.mapping_quality = 60
                seg.cigar = ((0, 100),)
                seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
                haplotype = row.get("haplotype")
                if haplotype is not None:
                    seg.set_tag("HP", int(haplotype))
                bam_out.write(seg)
        pysam.index(str(bam_path))

    def test_fetch_reads_filters_background_to_shared_support_haplotype(self):
        with tempfile.TemporaryDirectory() as td:
            bam_path = Path(td) / "reads.bam"
            self._write_test_bam(
                bam_path,
                [
                    {"read_name": "support1", "haplotype": 1},
                    {"read_name": "support2", "haplotype": 1, "start": 1010},
                    {"read_name": "background_hp1", "haplotype": 1, "start": 1020},
                    {"read_name": "background_hp2", "haplotype": 2, "start": 1030},
                ],
            )

            shown, all_reads = _fetch_reads(
                str(bam_path),
                "chr1",
                900,
                1300,
                supporting_reads={"support1", "support2"},
                max_reads=20,
                support_haplotype_only=True,
            )

            applied_haplotype = shown.attrs.get("applied_support_haplotype")
            self.assertEqual(applied_haplotype, 1)
            self.assertEqual(
                set(shown["read_name"].astype(str)),
                {"support1", "support2", "background_hp1"},
            )
            self.assertEqual(
                set(all_reads["read_name"].astype(str)),
                {"support1", "support2", "background_hp1", "background_hp2"},
            )

    def test_fetch_reads_keeps_unphased_supporting_reads_when_filtering(self):
        with tempfile.TemporaryDirectory() as td:
            bam_path = Path(td) / "reads.bam"
            self._write_test_bam(
                bam_path,
                [
                    {"read_name": "support1", "haplotype": 1},
                    {"read_name": "support_unphased", "haplotype": None, "start": 1010},
                    {"read_name": "background_hp1", "haplotype": 1, "start": 1020},
                    {"read_name": "background_hp2", "haplotype": 2, "start": 1030},
                ],
            )

            shown, _ = _fetch_reads(
                str(bam_path),
                "chr1",
                900,
                1300,
                supporting_reads={"support1", "support_unphased"},
                max_reads=20,
                support_haplotype_only=True,
            )

            applied_haplotype = shown.attrs.get("applied_support_haplotype")
            self.assertEqual(applied_haplotype, 1)
            self.assertEqual(
                set(shown["read_name"].astype(str)),
                {"support1", "support_unphased", "background_hp1"},
            )

    def test_fetch_reads_keeps_all_haplotypes_when_support_is_mixed(self):
        with tempfile.TemporaryDirectory() as td:
            bam_path = Path(td) / "reads.bam"
            self._write_test_bam(
                bam_path,
                [
                    {"read_name": "support1", "haplotype": 1},
                    {"read_name": "support2", "haplotype": 2, "start": 1010},
                    {"read_name": "background_hp1", "haplotype": 1, "start": 1020},
                    {"read_name": "background_hp2", "haplotype": 2, "start": 1030},
                ],
            )

            shown, _ = _fetch_reads(
                str(bam_path),
                "chr1",
                900,
                1300,
                supporting_reads={"support1", "support2"},
                max_reads=20,
                support_haplotype_only=True,
            )

            applied_haplotype = shown.attrs.get("applied_support_haplotype")
            self.assertIsNone(applied_haplotype)
            self.assertEqual(
                set(shown["read_name"].astype(str)),
                {"support1", "support2", "background_hp1", "background_hp2"},
            )


class TestVizLargeIndels(unittest.TestCase):
    def test_collect_large_indels_from_cigar_filters_and_clips(self):
        events = _collect_large_indels_from_cigar(
            cigartuples=[
                (0, 100),  # M
                (1, 55),   # INS
                (0, 20),
                (2, 80),   # DEL
                (0, 30),
                (1, 15),   # too small
            ],
            reference_start=1000,
            region_start=900,
            region_end=1200,
            min_indel_bp=40,
        )

        self.assertEqual(len(events), 2)
        self.assertEqual(events[0]["event_type"], "INS")
        self.assertEqual(events[0]["pos"], 1100)
        self.assertEqual(events[0]["length"], 55)
        self.assertEqual(events[1]["event_type"], "DEL")
        self.assertEqual(events[1]["start"], 1120)
        self.assertEqual(events[1]["end"], 1200)
        self.assertEqual(events[1]["length"], 80)


class TestVizLinkedCtdmrMode(unittest.TestCase):
    def test_expand_interval_for_visibility_enforces_minimum_span(self):
        new_start, new_end = _expand_interval_for_visibility(
            1000,
            1010,
            region_start=0,
            region_end=5000,
            min_span_bp=120,
        )
        self.assertEqual(new_end - new_start, 120)
        self.assertLessEqual(new_start, 1000)
        self.assertGreaterEqual(new_end, 1010)

    def test_extend_region_uses_nearest_informative_ctdmr(self):
        callouts = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [1200, 3100],
                "end": [1250, 3150],
                "callout_side": ["left", "right"],
                "callout_support_count": [1, 2],
                "callout_distance_bp": [650, 1100],
            }
        )

        new_start, new_end, row = _extend_region_to_first_informative_ctdmr(
            region_start=1900,
            region_end=2900,
            linked_ctdmr_callouts=callouts,
        )

        self.assertEqual(new_start, 700)
        self.assertEqual(new_end, 2900)
        self.assertIsNotNone(row)
        self.assertEqual(int(row["start"]), 1200)

    def test_assign_ctdmr_label_lanes_separates_overlapping_labels(self):
        entries = [
            {"x_center": 100.0, "text": "0.95"},
            {"x_center": 110.0, "text": "0.82"},
            {"x_center": 300.0, "text": "0.14"},
        ]

        lanes, lane_count = _assign_ctdmr_label_lanes(
            entries,
            region_span_bp=1000,
            font_size=12.0,
        )

        self.assertEqual(len(lanes), 3)
        self.assertGreaterEqual(lane_count, 2)
        self.assertNotEqual(lanes[0], lanes[1])
        self.assertEqual(lanes[2], 0)


if __name__ == "__main__":
    unittest.main()
