#!/usr/bin/env python
import argparse
import sys, os
from sniffcell.__init__ import __version__ as version


def _build_discover_run_parent_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--deconv-dir", required=True, help="Path to the sample deconv directory.")
    parser.add_argument("--reference", required=True, help="Reference FASTA used by downstream tools.")
    parser.add_argument("--tr-bed", required=True, help="Tandem repeat BED for medaka tandem.")
    parser.add_argument("--sex", required=True, choices=["female", "male"], help="Sample sex for medaka tandem.")
    parser.add_argument(
        "--scheduler",
        default="local",
        choices=["local", "slurm"],
        help="Execution mode. local runs sequentially; slurm renders or submits HPC scripts.",
    )
    parser.add_argument("--slurm-account", default=None, help="Slurm account written into the generated submit_pipeline.sh. Optional.")
    parser.add_argument("--split-dir", default=None, help="Optional override for the deconv_requested_group_splits directory.")
    parser.add_argument("--sample-id", default=None, help="Optional sample ID override.")
    parser.add_argument("--groups", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--stages", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--run-id", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--dry-run", action="store_true", default=False, help="Write manifests and commands without executing.")
    parser.add_argument("--force", action="store_true", default=False, help="Rerun stages even if done markers already exist.")
    parser.add_argument("--rerun-failed", action="store_true", default=False, help="Allow rerunning tasks that previously failed.")
    parser.add_argument("--sniffles-bin", default=None, help="Optional Sniffles executable path.")
    parser.add_argument("--bcftools-bin", default=None, help="Optional bcftools executable path.")
    parser.add_argument("--bgzip-bin", default=None, help="Optional bgzip executable path.")
    parser.add_argument("--kanpig-bin", default=None, help="Optional Kanpig executable path.")
    parser.add_argument("--truvari-bin", default=None, help="Optional Truvari executable path.")
    parser.add_argument("--medaka-bin", default=None, help="Optional Medaka executable path.")
    parser.add_argument("--tdb-bin", default=None, help="Optional tdb executable path.")
    parser.add_argument("--modkit-bin", default=None, help="Optional modkit executable path.")
    parser.add_argument("--tabix-bin", default=None, help="Optional tabix executable path.")
    parser.add_argument("--threads", type=int, default=16, help="Threads used by all tools (sniffles, kanpig, medaka, modkit, clair3, clairS, tdb merge). Default=16.")
    parser.add_argument(
        "--sniffles-mosaic-filter-expression",
        default="INFO/MOSAIC=1",
        help="bcftools expression used with -f PASS to derive the Kanpig handoff VCF. Default=INFO/MOSAIC=1.",
    )
    parser.add_argument("--sniffles-cluster-merge-len", type=float, default=0.2, help="Sniffles --cluster-merge-len. Default=0.2.")
    parser.add_argument("--kanpig-seqsim", type=float, default=0.8, help="Kanpig seqsim. Default=0.8.")
    parser.add_argument("--kanpig-sizesim", type=float, default=0.85, help="Kanpig sizesim. Default=0.85.")
    parser.add_argument("--kanpig-passonly", action="store_true", default=True, help="Use --passonly for Kanpig.")
    parser.add_argument(
        "--kanpig-sample-name-template",
        default="{sample_id}_{group}",
        help="Sample naming template for Kanpig. Default={sample_id}_{group}.",
    )
    parser.add_argument("--collapse-use", choices=["kanpig", "sniffles"], default="kanpig", help="Which callset to collapse. Default=kanpig.")
    parser.add_argument("--truvari-refdist", type=int, default=500, help="Truvari refdist. Default=500.")
    parser.add_argument("--truvari-pctseq", type=float, default=0.95, help="Truvari pctseq. Default=0.95.")
    parser.add_argument("--truvari-pctsize", type=float, default=0.95, help="Truvari pctsize. Default=0.95.")
    parser.add_argument("--truvari-passonly", action="store_true", default=True, help="Use --passonly for truvari collapse.")
    parser.add_argument(
        "--medaka-model",
        default="dna_r10.4.1_e8.2_400bps_sup@v4.3.0:consensus",
        help="Model for medaka tandem.",
    )
    parser.add_argument("--medaka-padding", type=int, default=250, help="Padding for medaka tandem. Default=250.")
    parser.add_argument(
        "--medaka-phasing",
        choices=["hybrid", "abpoa", "unphased", "prephased"],
        default=None,
        help="Optional medaka tandem phasing mode. When unset, medaka uses its own default.",
    )
    parser.add_argument(
        "--medaka-sample-name-template",
        default="{sample_id}.{group}",
        help="Sample naming template for medaka tandem. Default={sample_id}.{group}.",
    )
    parser.add_argument(
        "--tail-expansion-rescue",
        action="store_true",
        default=False,
        help=(
            "Experimental: let tr_post_processing rescue loci with little/no TDB hap delta when one group "
            "shows a hap-specific expansion tail in the trimmed reads."
        ),
    )
    parser.add_argument(
        "--tail-expansion-require-sample-range-support",
        dest="tail_require_sample_range_support",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "When tail-expansion rescue is enabled, require the rescued haplotype to also show same-haplotype "
            "TDB sample-range support. Use --no-tail-expansion-require-sample-range-support to allow "
            "trimmed-read-only tail rescue."
        ),
    )
    parser.add_argument("--tdb-create-mem", type=int, default=4, help="Memory in GB for tdb create. Default=4.")
    parser.add_argument("--tdb-create-force", action="store_true", default=False, help="Pass --force to tdb create.")
    parser.add_argument(
        "--mods-mode",
        choices=["combined", "separate"],
        default="separate",
        help="Output mode label for modkit-derived methylation summaries. Default=separate.",
    )
    parser.add_argument("--clair3-bin", default=None, help="Optional run_clair3.sh executable path.")
    parser.add_argument("--clair3-platform", default="ont", help="Sequencing platform for Clair3 (e.g. ont, hifi). Default=ont.")
    parser.add_argument("--clair3-model-path", default=None, help="Path to Clair3 model directory. Required when running the clair3 stage.")
    return parser


def _normalize_discover_argv(argv: list[str]) -> list[str]:
    if not argv or argv[0] != "discover":
        return argv
    if len(argv) == 1:
        return argv
    if argv[1] in {"tools", "ctprocessing", "-h", "--help"}:
        return argv
    return ["discover", "tools", "run", *argv[1:]]


def parse_args(argv):
    from sniffcell.discover.envcheck import _build_parser as _build_discover_envcheck_parser
    from sniffcell.discover.harmonize_variants import _build_arg_parser as _build_discover_harmonize_parser
    from sniffcell.discover.snv_post_processing import _build_arg_parser as _build_discover_snv_post_parser
    from sniffcell.discover.sv_post_processing import _build_arg_parser as _build_discover_sv_post_parser
    from sniffcell.discover.tr_post_processing import _build_arg_parser as _build_discover_tr_post_parser
    from sniffcell.sv_discovery import build_parser as _build_discover_sv_discovery_parser

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
    valid_commands = ["find", "deconv", "anno", "svanno", "dmsv", "viz", "igvviz", "report", "discover"]
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
    deconv_parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    deconv_parser.add_argument("-r", "--reference", required=True, help="Reference FASTA file")
    deconv_parser.add_argument("-b", "--bed", required=True, help="Input ctDMR BED/TSV file from sniffcell find.")
    deconv_parser.add_argument(
        "-o",
        "--output",
        required=True,
        help=(
            "Output folder, or an explicit TSV path for the overall summary. "
            "deconv also writes row-level/read-level companion tables beside it."
        ),
    )
    deconv_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use, default=1")
    deconv_parser.add_argument(
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
    deconv_parser.add_argument(
        "--split_bam_groups",
        type=str,
        default=None,
        help=(
            "Optional user-defined BAM/TSV split groups after deconvolution. "
            "Use ';' to separate output groups and ',' to separate cell types/labels within a group. "
            "Example: 't_cell,b_cell,nk_cell;monocyte' or 'lymph=t_cell,b_cell,nk_cell;myeloid=monocyte'. "
            "Labels are matched case-insensitively after punctuation normalization, and parent labels expand to their leaf subtypes."
        ),
    )
    deconv_parser.add_argument(
        "--resume",
        action="store_true",
        default=False,
        help=(
            "Resume post-processing from existing output TSVs. "
            "Reads deconv_reads_classification.tsv and deconv_read_summary.tsv from the output directory "
            "and skips the ctDMR phase, jumping directly to group splitting and BAM output."
        ),
    )
    deconv_parser.add_argument(
        "--per_read_min_agreement",
        type=float,
        default=0.66,
        help=(
            "Per-read consensus: minimum plurality fraction to accept a read's majority code "
            "when its ctDMRs have a genuine conflict (bitwise intersection = all zeros). "
            "Set to 1.0 to mark all conflicted reads as mixed. Default=0.66."
        ),
    )
    deconv_parser.add_argument(
        "--skip_overall_summary",
        action="store_true",
        default=False,
        help=(
            "Skip writing deconv_summary.tsv. Useful when deconv is only being used to produce "
            "per-read outputs and requested BAM splits, and whole-sample aggregation would add runtime."
        ),
    )

    # Subcommand: anno
    anno_parser = subparsers.add_parser("anno", help="Annotate variants with cell-type-specific methylation.")
    # add anno-specific args here
    anno_parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    anno_parser.add_argument(
        "-v", "--vcf", "--variants",
        required=True,
        help="Variant input for annotation. Accepts a harmonized TSV from sniffcell discover or a VCF file.",
    )
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
        "--breakpoint_exclusion_frac",
        type=float,
        default=0.0,
        help=(
            "Expand the no-ctDMR zone around the SV core by this fraction of absolute SV length on each side. "
            "Example: 0.1 excludes ctDMRs within +/-10% of SV length around breakpoints. Default=0.0"
        ),
    )
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
    anno_parser.add_argument(
        "--per_read_min_agreement",
        type=float,
        default=0.66,
        help=(
            "Per-read consensus: minimum plurality fraction to accept a read's majority code "
            "when its ctDMRs have a genuine conflict (bitwise intersection = all zeros). "
            "Set to 1.0 to mark all conflicted reads as mixed. Default=0.66."
        ),
    )

    svanno_parser = subparsers.add_parser(
        "svanno",
        help=(
            "SV annotation. Full mode matches the historical sniffcell anno workflow "
            "(BAM + VCF + ref + ctDMR BED). Reassignment mode recomputes from an existing reads_classification.tsv."
        ),
    )
    svanno_parser.add_argument(
        "-v", "--vcf", "--variants",
        required=True,
        help="Variant input. Full mode expects a VCF; reassignment mode also accepts a harmonized TSV.",
    )
    svanno_parser.add_argument("-i", "--input", required=True, help="Input BAM in full mode, or reads_classification.tsv in reassignment mode")
    svanno_parser.add_argument("-r", "--reference", required=False, default=None, help="Reference FASTA for full annotation mode")
    svanno_parser.add_argument(
        "-b", "--bed",
        required=False,
        default=None,
        help="Input ctDMR BED/TSV from sniffcell find for full annotation mode.",
    )
    svanno_parser.add_argument( "-krn", "--kanpig_read_names", type=str, default=None, help="Read names TSV from kanpig output, will use Sniffles read names if not sepecified." )
    svanno_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use in full annotation mode, default=1")
    svanno_parser.add_argument("-w", "--window", type=int, default=5000, help="Window size for variant-aware region matching, default=5000")
    svanno_parser.add_argument(
        "--breakpoint_exclusion_frac",
        type=float,
        default=0.0,
        help=(
            "Expand the no-ctDMR zone around the SV core by this fraction of absolute SV length on each side. "
            "Example: 0.1 excludes ctDMRs within +/-10% of SV length around breakpoints. Default=0.0"
        ),
    )
    svanno_parser.add_argument(
        "--read_assignment_mode",
        type=str,
        choices=["closest_reference_mean", "kmeans"],
        default="closest_reference_mean",
        help=(
            "How to assign each read to best_group vs other_group within each ctDMR in full annotation mode. "
            "'closest_reference_mean' compares per-read methylation mean to mean_best_value/mean_rest_value, "
            "'kmeans' uses unsupervised clustering. Default=closest_reference_mean."
        ),
    )
    svanno_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output folder in full mode, or output folder / TSV path in reassignment mode.",
    )
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
    svanno_parser.add_argument(
        "--per_read_min_agreement",
        type=float,
        default=0.66,
        help=(
            "Per-read consensus: minimum plurality fraction to accept a read's majority code "
            "when its ctDMRs have a genuine conflict (bitwise intersection = all zeros). "
            "Set to 1.0 to mark all conflicted reads as mixed. Default=0.66."
        ),
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
        "--indel_min_bp",
        type=int,
        default=40,
        help="Show insertion/deletion events on reads when CIGAR length >= this threshold; set 0 to disable, default=40",
    )
    viz_parser.add_argument(
        "--skip_methylation_overlay",
        action="store_true",
        help="Skip per-read ctDMR methylation extraction/overlay for faster rendering.",
    )
    viz_parser.add_argument(
        "--support_haplotype_only",
        dest="support_haplotype_only",
        action="store_true",
        default=True,
        help=(
            "When phased SV-supporting reads agree on one haplotype, keep only reads from that "
            "haplotype in the panel while retaining supporting reads. default=on"
        ),
    )
    viz_parser.add_argument(
        "--show_all_haplotypes",
        dest="support_haplotype_only",
        action="store_false",
        help="Disable support-haplotype filtering and show reads from all haplotypes.",
    )
    viz_parser.add_argument(
        "--linked_ctdmr_mode",
        choices=["distal", "extend", "strict"],
        default="distal",
        help=(
            "How to handle winning linked ctDMRs outside the local display window: "
            "distal=show side callouts with dashed extensions; "
            "extend=expand to the nearest informative linked ctDMR; "
            "strict=keep the requested window and ignore off-window linked ctDMRs. "
            "default=distal"
        ),
    )
    viz_parser.add_argument(
        "--no_distal_ctdmr_callouts",
        dest="linked_ctdmr_mode",
        action="store_const",
        const="strict",
        help="Alias for `--linked_ctdmr_mode strict`.",
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
            "Generate an HTML report for high-confidence variant assignments. "
            "By default this is figure-less (fast); add --with_figures to render viz panels."
        ),
    )
    report_parser.add_argument(
        "--anno_output",
        required=True,
        help="sniffcell anno output folder containing variant_assignment.tsv and anno_run_manifest.json.",
    )
    report_parser.add_argument(
        "--min_overlap_pct",
        type=float,
        default=0.8,
        help="Base overlap threshold for including a variant in report, default=0.8",
    )
    report_parser.add_argument(
        "--overlap_filter_mode",
        choices=["gradient", "hard_clip"],
        default="gradient",
        help=(
            "How report filters overlap support. "
            "'gradient' scales the required n_overlapped sublinearly with n_supporting, "
            "so larger SV support can tolerate a lower overlap fraction. "
            "'hard_clip' uses the fixed overlap_pct threshold directly. Default=gradient"
        ),
    )
    report_parser.add_argument(
        "--overlap_gradient_exponent",
        type=float,
        default=0.5,
        help=(
            "Exponent used by report overlap_filter_mode=gradient. "
            "Required overlapped reads are ceil(min_overlap_pct * n_supporting^exponent). "
            "Lower values are more permissive for larger n_supporting. Default=0.5"
        ),
    )
    report_parser.add_argument(
        "--min_majority_pct",
        type=float,
        default=0.8,
        help="Minimum majority_pct threshold for including a variant in report, default=0.8",
    )
    report_parser.add_argument(
        "--include_unassigned",
        action="store_true",
        help="Include variants with empty assigned_code (default filters to assigned rows only, except linked TR rows).",
    )
    report_parser.add_argument(
        "--allow_hard_conflict",
        action="store_true",
        help="Include non-TR variants with has_hard_conflict=True (default excludes them).",
    )
    report_parser.add_argument(
        "--max_sv",
        type=int,
        default=0,
        help="Maximum number of variants to include after filtering. 0 means no limit.",
    )
    report_parser.add_argument(
        "--with_figures",
        action="store_true",
        help=(
            "Render viz panel figures for selected variants. "
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
        help="Figure format for report panels, default=png",
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
        help="Reuse existing per-variant viz figure files when present instead of regenerating.",
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
        help="Render IGV screenshots for selected variants using sniffcell igvviz.",
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
        help="Generate an alternate igv-reports HTML page for the selected variants.",
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

    discover_parser = subparsers.add_parser(
        "discover",
        help="Postprocess deconv split BAM outputs with SV, TR, and methylation workflows.",
    )
    discover_subparsers = discover_parser.add_subparsers(dest="discover_section")

    discover_tools_parser = discover_subparsers.add_parser(
        "tools",
        help="Run discover tool-facing commands.",
    )
    discover_tools_subparsers = discover_tools_parser.add_subparsers(dest="discover_tools_command")
    discover_tools_run_parser = discover_tools_subparsers.add_parser(
        "run",
        help="Run the main discover external-tool pipeline.",
        parents=[_build_discover_run_parent_parser()],
    )
    discover_tools_sv_parser = discover_tools_subparsers.add_parser(
        "sv",
        help="Run standalone SV discovery before discover/anno.",
        parents=[_build_discover_sv_discovery_parser(prog="sniffcell discover tools sv", add_help=False)],
    )
    discover_tools_check_parser = discover_tools_subparsers.add_parser(
        "check",
        help="Preflight-check discover external dependencies.",
        parents=[_build_discover_envcheck_parser(prog="sniffcell discover tools check", add_help=False)],
    )

    discover_ctprocessing_parser = discover_subparsers.add_parser(
        "ctprocessing",
        help="Run discover cell-type post-processing utilities.",
    )
    discover_ctprocessing_subparsers = discover_ctprocessing_parser.add_subparsers(dest="discover_ctprocessing_command")
    discover_ctprocessing_snv_parser = discover_ctprocessing_subparsers.add_parser(
        "snv",
        help="Run SNP post-processing on two Clair3 gVCF groups.",
        parents=[_build_discover_snv_post_parser(prog="sniffcell discover ctprocessing snv", add_help=False)],
    )
    discover_ctprocessing_sv_parser = discover_ctprocessing_subparsers.add_parser(
        "sv",
        help="Run SV post-processing on two split groups.",
        parents=[_build_discover_sv_post_parser(prog="sniffcell discover ctprocessing sv", add_help=False)],
    )
    discover_ctprocessing_tr_parser = discover_ctprocessing_subparsers.add_parser(
        "tr",
        help="Run tandem-repeat post-processing on two split groups.",
        parents=[_build_discover_tr_post_parser(prog="sniffcell discover ctprocessing tr", add_help=False)],
    )
    discover_ctprocessing_harmonize_parser = discover_ctprocessing_subparsers.add_parser(
        "harmonize",
        help="Merge TR and SV post-processing outputs into a harmonized variant TSV.",
        parents=[_build_discover_harmonize_parser(prog="sniffcell discover ctprocessing harmonize", add_help=False)],
    )


    top_level_passthrough = {"-h", "--help", "-v", "--version"}

    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif argv[0] in top_level_passthrough:
        return parser.parse_args(argv)
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
    elif len(argv) == 1 and argv[0] == "discover":
        discover_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 2 and argv[0] == "discover" and argv[1] == "tools":
        discover_tools_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 2 and argv[0] == "discover" and argv[1] == "ctprocessing":
        discover_ctprocessing_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "tools", "run"]:
        discover_tools_run_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "tools", "sv"]:
        discover_tools_sv_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "tools", "check"]:
        discover_tools_check_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "ctprocessing", "sv"]:
        discover_ctprocessing_sv_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "ctprocessing", "snv"]:
        discover_ctprocessing_snv_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "ctprocessing", "tr"]:
        discover_ctprocessing_tr_parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) == 3 and argv[:3] == ["discover", "ctprocessing", "harmonize"]:
        discover_ctprocessing_harmonize_parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(_normalize_discover_argv(list(argv)))
    return args
