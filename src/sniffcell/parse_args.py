#!/usr/bin/env python
import argparse
import sys, os
from sniffcell.__init__ import __version__ as version


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="sniffcell",
        description="Annotating mosaic structural variants (SVs) with cell type-specific methylation information.",
        epilog=f"Version {version}",
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"sniffcell {version}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    valid_commands = ["find", "deconv", "anno", "svanno", "dmsv", "viz", "igvviz", "report"]
    # Subcommand: find
    find_parser = subparsers.add_parser("find", help="Find cell type-specific DMRs.")
    atlas_dir = os.path.abspath("atlas")
    default_npy   = os.path.join(atlas_dir, "all_celltypes_blocks.npy")
    default_ct    = os.path.join(atlas_dir, "index_to_major_celltypes.json")
    default_index = os.path.join(atlas_dir, "all_celltypes_blocks.index.gz")
    default_meta  = os.path.join(atlas_dir, "all_celltypes.txt")


    find_parser.add_argument("-n", "--npy", default=default_npy, help=f"Input .npy matrix for finding cell type DMRs, default={default_npy}")
    find_parser.add_argument("-i", "--index", default=default_index, help=f"Index for CpGs in the npy matrix, default={default_index}")
    find_parser.add_argument("-cf", "--celltypes_file", default=default_ct, help=f"Cell type json files mapped to the major cell types, default={default_ct}")
    find_parser.add_argument("-m", "--meta", default=default_meta, help=f"Metadata file for cell types in the npy matrix, default={default_meta}")
    find_parser.add_argument("-ck", "--celltypes_keys", required=True, help="keys for major cell types in the cell type json file")
    find_parser.add_argument("-o", "--output", required=True, help="Output BED files for cell type DMRs")

    find_parser.add_argument( "--diff_threshold", type=float, default=0.40, help="Minimum difference threshold for calling DMRs, default=0.40" )
    find_parser.add_argument( "--min_rows", type=int, default=2, help="Minimum number of rows (CpG groups in index) for calling DMRs, default=2")
    find_parser.add_argument( "--min_cpgs", type=int, default=3, help="Minimum number of CpGs for calling DMRs, default=3" )
    find_parser.add_argument( "--max_gap_bp", type=int, default=500, help="Maximum gap among groups for calling DMRs, default=500" )


    # Subcommand: deconv
    deconv_parser = subparsers.add_parser("deconv", help="Deconvolve cell-type composition from methylation data.")
    # add deconv-specific args here
    deconv_parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    deconv_parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file")
    deconv_parser.add_argument("-b", "--bed", required=True, help="Input BED file with DMR indications")
    deconv_parser.add_argument("-o", "--output", required=True, help="Output file")
    
    # Subcommand: anno
    anno_parser = subparsers.add_parser("anno", help="Annotate variants with cell-type-specific methylation.")
    # add anno-specific args here
    anno_parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    anno_parser.add_argument("-v", "--vcf", required=True, help="Input VCF file for variant annotation")
    anno_parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file")
    anno_parser.add_argument(
        "-b", "--bed",
        required=True,
        help="Input ctDMR BED/TSV file from sniffcell find.",
    )
    anno_parser.add_argument("-o", "--output", required=True, help="Output folder")
    anno_parser.add_argument( "-krn", "--kanpig_read_names", type=str, default=None, help="Read names TSV from kanpig output, will use Sniffles read names if not sepecified." )
    anno_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use, default=1")
    anno_parser.add_argument("-w", "--window", type=int, default=5000, help="Window size for filtering BED based on variants, default=5000")
    anno_parser.add_argument(
        "--read_assignment_mode",
        type=str,
        choices=["closest_reference_mean", "kmeans"],
        default="closest_reference_mean",
        help=(
            "How to assign each read to best_group vs other_group within each ctDMR. "
            "'closest_reference_mean' compares per-read methylation mean to mean_best_value/mean_rest_value, "
            "'kmeans' uses unsupervised clustering. Default=closest_reference_mean."
        ),
    )
    anno_parser.add_argument(
        "--evidence_mode",
        type=str,
        choices=["all_rows", "per_read"],
        default="all_rows",
        help=(
            "How to aggregate ctDMR evidence for SV assignment: "
            "'all_rows' uses every supporting-read x ctDMR row, "
            "'per_read' keeps one winning code per read. Default=all_rows."
        ),
    )
    anno_parser.add_argument(
        "--min_overlap_pct",
        type=float,
        default=0.0,
        help="Minimum overlapped evidence fraction required to keep assigned_code, default=0.0",
    )
    anno_parser.add_argument(
        "--min_agreement_pct",
        type=float,
        default=1.0,
        help="Minimum majority agreement fraction required to keep assigned_code, default=1.0",
    )
    
    svanno_parser = subparsers.add_parser("svanno", help="Use pre-annotated reads csv to annotate variants' cell types")
    svanno_parser.add_argument("-v", "--vcf", required=True, help="Input VCF file for variant annotation")
    svanno_parser.add_argument("-i", "--input", required=True, help="Input reads_classification.tsv file from anno step")
    svanno_parser.add_argument( "-krn", "--kanpig_read_names", type=str, default=None, help="Read names TSV from kanpig output, will use Sniffles read names if not sepecified." )
    svanno_parser.add_argument("-w", "--window", type=int, default=5000, help="Window size for SV-aware region matching, default=5000")
    svanno_parser.add_argument("-o", "--output", required=True, help="Output TSV file path for sv_assignment")
    svanno_parser.add_argument(
        "--evidence_mode",
        type=str,
        choices=["all_rows", "per_read"],
        default="all_rows",
        help=(
            "How to aggregate ctDMR evidence for SV assignment: "
            "'all_rows' uses every supporting-read x ctDMR row, "
            "'per_read' keeps one winning code per read. Default=all_rows."
        ),
    )
    svanno_parser.add_argument(
        "--min_overlap_pct",
        type=float,
        default=0.0,
        help="Minimum overlapped evidence fraction required to keep assigned_code, default=0.0",
    )
    svanno_parser.add_argument(
        "--min_agreement_pct",
        type=float,
        default=1.0,
        help="Minimum majority agreement fraction required to keep assigned_code, default=1.0",
    )

    dmsv_parser = subparsers.add_parser("dmsv", help="Find out which SV's supporting reads have differential methylation compared to non-supporting reads.")
    dmsv_parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    dmsv_parser.add_argument("-v", "--vcf", required=True, help="Input VCF file for variant annotation")
    dmsv_parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file")
    dmsv_parser.add_argument("-c", "--min_cpgs", type=int, default=5, help="Minimum number of CpGs in the flanking region to consider an SV is causing methylation changes, default=5")
    dmsv_parser.add_argument("-f", "--flank_size", type=int, default=1000, help="Number of base pairs to flank on both sides of the SV, default=1000")
    dmsv_parser.add_argument( "--test_type", type=str, default="t-test", choices=["t-test", "mannwhitneyu", "fisher"], help="Type of statistical test to perform. Options are 't-test', 'mannwhitneyu', 'fisher'. Default is 't-test'" ) 
    dmsv_parser.add_argument( "--haplotype_majority_threshold", type=float, default=0.7, help="Threshold for majority haplotype in supporting reads to perform statistical test, default=0.7" )
    dmsv_parser.add_argument("-m", "--min_supporting", type=int, default=3, help="Minimum supporting reads for SV. default 3.")
    dmsv_parser.add_argument("-o", "--output", required=True, help="Output folder")
    dmsv_parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use, default=4")

    viz_parser = subparsers.add_parser("viz", help="Visualize one SV with nearby reads and ctDMRs.")
    viz_parser.add_argument(
        "--anno_output",
        default=None,
        help="sniffcell anno output folder; if set, viz can auto-load BAM/VCF/REF/BED/read assignments from its run manifest.",
    )
    viz_parser.add_argument("-i", "--input", required=False, default=None, help="Input BAM file")
    viz_parser.add_argument("-v", "--vcf", required=False, default=None, help="Input VCF file")
    viz_parser.add_argument("-s", "--sv_id", required=True, help="SV ID to visualize")
    viz_parser.add_argument(
        "-r", "--reference",
        default=None,
        help="Reference FASTA (needed for per-read methylation on ctDMRs).",
    )
    viz_parser.add_argument(
        "-b", "--bed",
        default=None,
        help="Optional ctDMR BED/TSV file from sniffcell find.",
    )
    viz_parser.add_argument(
        "-a", "--read_assignment",
        default=None,
        help="Optional reads_classification.tsv from anno for assigned/unassigned supporting-read summaries.",
    )
    viz_parser.add_argument(
        "-krn", "--kanpig_read_names",
        type=str,
        default=None,
        help="Optional TSV mapping SV IDs to supporting read names.",
    )
    viz_parser.add_argument(
        "-w", "--window",
        type=int,
        default=5000,
        help="Window size around SV to plot, default=5000",
    )
    viz_parser.add_argument(
        "--exact_window",
        action="store_true",
        help=(
            "Use --window exactly as provided. "
            "By default, when --window=5000 and --anno_output is used, viz adopts anno manifest window if available."
        ),
    )
    viz_parser.add_argument(
        "-m", "--max_reads",
        type=int,
        default=250,
        help="Maximum reads to draw, prioritizing supporting reads, default=250",
    )
    viz_parser.add_argument(
        "-f", "--format",
        default="png",
        choices=["png", "pdf"],
        help="Output figure format, default=png",
    )
    viz_parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output figure DPI, default=300",
    )
    viz_parser.add_argument(
        "--skip_methylation_overlay",
        action="store_true",
        help="Skip per-read ctDMR methylation extraction/overlay for faster rendering.",
    )
    viz_parser.add_argument(
        "--export_tables",
        action="store_true",
        help="Write supplementary TSV tables in addition to the figure.",
    )
    viz_parser.add_argument(
        "-o", "--output",
        required=False,
        default=None,
        help="Output figure path or prefix. Defaults to <anno_output>/<sv_id>.viz.<format> when --anno_output is set.",
    )

    igvviz_parser = subparsers.add_parser(
        "igvviz",
        help="Use IGV batch mode to render screenshots of one SV across one or more BAM files.",
    )
    igvviz_parser.add_argument(
        "--anno_output",
        default=None,
        help="sniffcell anno output folder; if set, igvviz can auto-load BAM/VCF/REF/BED from anno_run_manifest.json.",
    )
    igvviz_parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=False,
        default=None,
        help="One or more input BAM files. If omitted, falls back to anno manifest BAM when --anno_output is set.",
    )
    igvviz_parser.add_argument("-v", "--vcf", required=False, default=None, help="Input VCF file")
    igvviz_parser.add_argument("-s", "--sv_id", required=True, help="SV ID to visualize")
    igvviz_parser.add_argument(
        "-r",
        "--reference",
        default=None,
        help="Reference FASTA file to load in IGV.",
    )
    igvviz_parser.add_argument(
        "-b",
        "--bed",
        default=None,
        help="Optional ctDMR BED/TSV file. Overlapping ctDMRs are loaded as an IGV marker track.",
    )
    igvviz_parser.add_argument(
        "-krn",
        "--kanpig_read_names",
        type=str,
        default=None,
        help="Optional TSV mapping SV IDs to supporting read names; overrides VCF RNAMES when provided.",
    )
    igvviz_parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=5000,
        help="Window size around SV to screenshot, default=5000",
    )
    igvviz_parser.add_argument(
        "--visibility_window",
        type=int,
        default=None,
        help="Optional IGV alignment visibility window (bp). Defaults to --window when omitted.",
    )
    igvviz_parser.add_argument(
        "--phase_tag",
        default="HP",
        help="Tag used for IGV grouping by phase (default=HP).",
    )
    igvviz_parser.add_argument(
        "--support_tag",
        default="SC",
        help="Temporary BAM tag used to mark supporting reads (default=SC).",
    )
    igvviz_parser.add_argument(
        "--igv_cmd",
        default="igv.sh",
        help="IGV executable command used to run batch mode (default=igv.sh).",
    )
    igvviz_parser.add_argument(
        "--batch_only",
        action="store_true",
        help="Only write batch/intermediate files; do not execute IGV.",
    )
    igvviz_parser.add_argument(
        "--keep_intermediates",
        action="store_true",
        help="Keep intermediate tagged BAM/BED files under output directory.",
    )
    igvviz_parser.add_argument(
        "--snapshot_format",
        choices=["png", "jpg", "svg"],
        default="png",
        help="Snapshot image format written by IGV, default=png",
    )
    igvviz_parser.add_argument(
        "--snapshot_width",
        type=int,
        default=3600,
        help="IGV window width in pixels for snapshots, default=3600",
    )
    igvviz_parser.add_argument(
        "--snapshot_height",
        type=int,
        default=1600,
        help="IGV window height in pixels for snapshots, default=1600",
    )
    igvviz_parser.add_argument(
        "--hide_methylation",
        action="store_true",
        help="Do not use IGV base-modification coloring in screenshots.",
    )
    igvviz_parser.add_argument(
        "-o",
        "--output",
        required=False,
        default=None,
        help="Output directory for snapshots and IGV batch files. Defaults to <anno_output>/igvviz or ./igvviz.",
    )

    report_parser = subparsers.add_parser(
        "report",
        help=(
            "Generate an HTML report for high-confidence SV assignments. "
            "By default this is figure-less (fast); add --with_figures to render viz panels."
        ),
    )
    report_parser.add_argument(
        "--anno_output",
        required=True,
        help="sniffcell anno output folder containing sv_assignment.tsv and anno_run_manifest.json.",
    )
    report_parser.add_argument(
        "--min_overlap_pct",
        type=float,
        default=0.8,
        help="Minimum overlap_pct threshold for including an SV in report, default=0.8",
    )
    report_parser.add_argument(
        "--min_majority_pct",
        type=float,
        default=1.0,
        help="Minimum majority_pct threshold for including an SV in report, default=1.0",
    )
    report_parser.add_argument(
        "--include_unassigned",
        action="store_true",
        help="Include SVs with empty assigned_code (default filters to assigned SVs only).",
    )
    report_parser.add_argument(
        "--allow_hard_conflict",
        action="store_true",
        help="Include SVs with has_hard_conflict=True (default excludes them).",
    )
    report_parser.add_argument(
        "--max_sv",
        type=int,
        default=0,
        help="Maximum number of SVs to include after filtering. 0 means no limit.",
    )
    report_parser.add_argument(
        "--with_figures",
        action="store_true",
        help=(
            "Render viz panel figures for selected SVs. "
            "Default is figure-less report (no viz rendering)."
        ),
    )
    report_parser.add_argument(
        "-w", "--window",
        type=int,
        default=5000,
        help="Window size passed through to viz; when default=5000 and anno manifest has a window, viz uses the manifest window.",
    )
    report_parser.add_argument(
        "-m", "--max_reads",
        type=int,
        default=250,
        help="Maximum reads per viz panel, default=250",
    )
    report_parser.add_argument(
        "-f", "--format",
        default="png",
        choices=["png", "pdf"],
        help="Figure format for SV panels, default=png",
    )
    report_parser.add_argument(
        "--figure_profile",
        choices=["fast", "full"],
        default="full",
        help="Figure rendering profile for report panels. 'full' keeps read-level methylation overlay; 'fast' skips it for speed.",
    )
    report_parser.add_argument(
        "--figure_dpi",
        type=int,
        default=160,
        help="Figure DPI for report panel rendering, default=160",
    )
    report_parser.add_argument(
        "--reuse_existing_viz",
        action="store_true",
        help="Reuse existing per-SV viz figure files when present instead of regenerating.",
    )
    report_parser.add_argument(
        "--figure_threads",
        type=int,
        default=1,
        help="Shared thread count for figure and igvviz rendering, default=1",
    )
    report_parser.add_argument(
        "--with_igvviz",
        action="store_true",
        help="Render IGV screenshots for selected SVs using sniffcell igvviz.",
    )
    report_parser.add_argument(
        "--igv_bams",
        nargs="+",
        default=None,
        help="One or more BAM files for igvviz. If omitted, igvviz uses anno manifest BAM.",
    )
    report_parser.add_argument(
        "--igv_cmd",
        default="igv.sh",
        help="IGV executable command for igvviz batch rendering, default=igv.sh",
    )
    report_parser.add_argument(
        "--igv_snapshot_format",
        choices=["png", "jpg", "svg"],
        default="png",
        help="Snapshot format for igvviz outputs, default=png",
    )
    report_parser.add_argument(
        "--igv_snapshot_width",
        type=int,
        default=3600,
        help="IGV window width in pixels for igvviz snapshots, default=3600",
    )
    report_parser.add_argument(
        "--igv_snapshot_height",
        type=int,
        default=1600,
        help="IGV window height in pixels for igvviz snapshots, default=1600",
    )
    report_parser.add_argument(
        "--reuse_existing_igvviz",
        action="store_true",
        help="Reuse existing igvviz manifests/snapshots when available instead of regenerating.",
    )
    report_parser.add_argument(
        "--with_igvreport",
        action="store_true",
        help="Generate an alternate igv-reports HTML page for the selected SVs.",
    )
    report_parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Report output directory or HTML file path. "
            "If output ends with .gz/.tgz/.tar.gz, report files are written to a sibling folder and "
            "a gzipped tar archive is also produced at the provided output path. "
            "Defaults to <anno_output>/report/."
        ),
    )


    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) >= 1 and argv[0] not in valid_commands:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "find":
        find_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "deconv":
        deconv_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "anno":
        anno_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "svanno":
        svanno_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "dmsv":
        dmsv_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "viz":
        viz_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "igvviz":
        igvviz_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 1 and argv[0] == "report":
        report_parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)
    return args
