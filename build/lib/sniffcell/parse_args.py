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
    anno_parser.add_argument("-b", "--bed", required=True, help="Input BED file with DMR indications")
    anno_parser.add_argument("-o", "--output", required=True, help="Output folder")
    anno_parser.add_argument( "-krn", "--kanpig_read_names", type=str, default=None, help="Read names TSV from kanpig output, will use Sniffles read names if not sepecified." )
    anno_parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use, default=1")
    anno_parser.add_argument("-w", "--window", type=int, default=5000, help="Window size for filtering BED based on variants, default=5000")
    
    svanno_parser = subparsers.add_parser("svanno", help="Use pre-annotated reads csv to annotate variants' cell types")
    svanno_parser.add_argument("-v", "--vcf", required=True, help="Input VCF file for variant annotation")
    svanno_parser.add_argument("-i", "--input", required=True, help="Input reads_classification.tsv file from anno step")
    svanno_parser.add_argument( "-krn", "--kanpig_read_names", type=str, default=None, help="Read names TSV from kanpig output, will use Sniffles read names if not sepecified." )
    svanno_parser.add_argument("-o", "--output", required=True, help="Output sv_assignment.tsv file")

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


    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(argv) >= 1 and argv[0] not in ["find", "deconv", "anno", "svanno", "dmsv"]:
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
    args = parser.parse_args(argv)
    return args
