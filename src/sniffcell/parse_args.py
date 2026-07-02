#!/usr/bin/env python
import argparse
import os

from sniffcell.__init__ import __version__ as version
from sniffcell.tissue_atlas import default_tissue_atlas_path


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="sniffcell-lite",
        description="Lite SniffCell interface for tissue ctDMR discovery and variant annotation.",
        epilog=f"Version {version}",
    )
    parser.add_argument("-v", "--version", action="version", version=f"sniffcell-lite {version}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    atlas_dir = os.path.abspath("atlas")
    default_npy = os.path.join(atlas_dir, "all_celltypes_blocks.npy")
    default_index = os.path.join(atlas_dir, "all_celltypes_blocks.index.gz")
    default_meta = os.path.join(atlas_dir, "all_celltypes.txt")

    find_parser = subparsers.add_parser("find", help="Find cell type-specific DMRs.")
    find_parser.add_argument("-n", "--npy", default=default_npy, help=f"Input .npy matrix for finding cell type DMRs, default={default_npy}")
    find_parser.add_argument("-i", "--index", default=default_index, help=f"Index for CpGs in the npy matrix, default={default_index}")
    find_parser.add_argument("-cf", "--celltypes_file", default=default_tissue_atlas_path(), help=f"Cell type atlas JSON, default={default_tissue_atlas_path()}")
    find_parser.add_argument("-m", "--meta", default=default_meta, help=f"Metadata file for cell types in the npy matrix, default={default_meta}")
    find_parser.add_argument("-ck", "--celltypes_keys", required=True, help="Atlas key, tissue code, or tissue name to use from the cell type atlas JSON.")
    find_parser.add_argument("-o", "--output", required=True, help="Output BED/TSV for cell type DMRs.")
    find_parser.add_argument("--diff_threshold", type=float, default=0.35, help="Minimum difference threshold for calling DMRs, default=0.35")
    find_parser.add_argument("--min_rows", type=int, default=2, help="Minimum number of rows for calling DMRs, default=2")
    find_parser.add_argument("--min_cpgs", type=int, default=3, help="Minimum number of CpGs for calling DMRs, default=3")
    find_parser.add_argument("--max_gap_bp", type=int, default=2000, help="Maximum gap among groups for calling DMRs, default=2000")

    anno_parser = subparsers.add_parser("anno", help="Annotate variants by computing supporting-read methylation at ctDMRs carried by the supporting reads.")
    anno_parser.add_argument("-i", "--input", required=False, default=None, help="Input BAM for single-variant annotation.")
    anno_parser.add_argument("-r", "--reference", required=False, default=None, help="Reference FASTA for CpG lookup and BAM methylation extraction.")
    anno_parser.add_argument("-o", "--output", required=True, help="Output folder.")
    anno_parser.add_argument("--variant-name", default=None, help="Single-variant mode: variant identifier/name.")
    anno_parser.add_argument("--variant-location", default=None, help="Single-variant mode: variant location as chr:pos or chr:start-end, using 1-based coordinates.")
    anno_parser.add_argument("--supporting-reads", default=None, help="Single-variant mode: supporting read names as comma/pipe/semicolon-delimited text, JSON list, or @file.")
    anno_parser.add_argument("--catalog", default=None, help="Single-variant mode: ctDMR catalog TSV from sniffcell-lite find.")
    anno_parser.add_argument("--batch", default=None, help="Batch TSV/CSV with columns: variant_name, variant_location, supporting_reads, catalog, bam, and optionally reference.")
    anno_parser.add_argument("--evidence_mode", type=str, choices=["all_rows", "per_read"], default="per_read", help="How to aggregate ctDMR evidence. Default=per_read.")
    anno_parser.add_argument("--min_overlap_pct", type=float, default=0.0, help="Minimum overlapped evidence fraction required to keep assigned_code, default=0.0")
    anno_parser.add_argument("--min_agreement_pct", type=float, default=0.0, help="Minimum majority agreement fraction required to keep assigned_code, default=0.0")
    anno_parser.add_argument("--per_read_min_agreement", type=float, default=0.66, help="Minimum plurality fraction for conflicted per-read consensus, default=0.66")
    anno_parser.add_argument("--window", type=int, default=10000, help="BAM fetch padding around each variant for finding supporting-read alignments. ctDMR evidence is selected from the mapped support-read spans, not limited to this distance. Default=10000")
    anno_parser.add_argument("--breakpoint_exclusion_frac", type=float, default=0.0, help="Breakpoint exclusion fraction retained for assignment compatibility. Default=0.0")

    return parser.parse_args(argv)
