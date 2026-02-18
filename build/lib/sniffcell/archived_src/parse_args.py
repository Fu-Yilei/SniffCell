#!/usr/bin/env python
import argparse
import sys
import os
from src.__init__ import __version__ as version


def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="sniffcell",
        description="Annotating mosaic structural variants (SVs) with cell type-specific methylation information.",
        epilog=f"Version {version}",
    )
    
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads, default 1.")
    

    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-b", "--bam", type=str, required=True, help="Input BAM file.")
    required_args.add_argument("-v", "--vcf", type=str, required=True, help="Input SV VCF file (supporting reads needed).")
    required_args.add_argument("-r", "--reference", type=str, required=True, help="Reference genome.")
    required_args.add_argument("-o", "--output", type=str, required=True, default=f"{os.path.abspath("./")}", help=f"Output directory. Default: {os.path.abspath('./')}")
    required_args.add_argument("-d", "--bed", type=str, required=True, help="Input cell type DMR BED file got from atlas.")
    optional_args = parser.add_argument_group("SniffCell deconvolution optional arguments")
    script_dir = os.path.abspath(os.path.dirname(__file__))
    atlas_celltypes_file_path = os.path.abspath(os.path.join(script_dir, "atlas", "39Bisulfite.tsv"))
    json_file_file_path = os.path.join(script_dir, "atlas", "tissue_celltypes.json")
    optional_args.add_argument("-a", "--atlas", type=str, default=atlas_celltypes_file_path, help=f"Cell type atlas file location, cell type JSON file need be at the same directory. Default: {atlas_celltypes_file_path} and {json_file_file_path}.")
    optional_args.add_argument("-c", "--tissue", type=str, required=True, help=f"Tissue type in {json_file_file_path}, need to update based on atlas. Required.")
    optional_args.add_argument("-n", "--region_number", type=int, default=300, help="Number of regions to be selected, default 300.")
    optional_args.add_argument("-me", "--method", type=str, default="std", help="Region selection method: std or diff, default std. Diff selects the regions with certain cell type as low methylation while all other cell types have high methylation. std selects regions with highest methylation value std in all cell types.")
    optional_args.add_argument("-vb", "--verbose", action="store_true", help="Enable verbose mode.")
    optional_args.add_argument("-ot", "--outlier_threshold", type=float, default=0.8, help="deviation threshold for filtering out wrongly deconvoluted regions. Default 0.8.")
    optional_args.add_argument("-confem", "--confidence", type=float, default=0.9, help="Minimum confidence threshold for EM algorithm. Default 0.9.")
    optional_args.add_argument("-confsnp", "--confidence_snp", type=float, default=0.9, help="Minimum confidence threshold for SNP filtering. Default 0.9.")
    optional_args.add_argument("-nc", "--n_closest", type=int, default=2, help="Number of closest methylation informative region for proportion estimation.")
    optional_args.add_argument("-wuoff",  "--wgbs_tools_uxm",  action="store_false", default=False, help="Use WGBS tools and UXM to set prior for EM algorithm. Default is False, but can be turned on.")
    optional_args.add_argument("-dis", "--hp_distance", type=float, default=0.3, help="Minimum distance of celltype probablity vector among hp1, hp2 and all hp. Default is 0.3. If the value is too large it means this region is confounded by allele specific methylation, hence discarded.")
    optional_args.add_argument("-rf", "--assigned_read_fraction", type=float, default=0.8, help="Minimum fraction of assigned reads inside a ctDMR (based on gemma value from EM algorithm). If this value is too low, it means the ctDMR is not reliable, hence discarded. Default is 0.8.")

    sniffcell_args = parser.add_argument_group("SniffCell optional arguments")
    sniffcell_args.add_argument(
        "-snp",
        "--snp",
        action="store_true",
        help="Indicate that the input VCF file contains SNPs. SniffCell will ignore SNPs with allele frequency between 0.4–0.6 and 0.8–1.0."
    )
    sniffcell_args.add_argument("-ob", "--output_bam", action="store_true", help="Output SV-related BAM file.")
    sniffcell_args.add_argument("-p", "--primary_only", action="store_true", help="Keep only primary alignments in BAM.")
    sniffcell_args.add_argument("-s", "--smoothing", type=int, default=0, help="Enable smoothing (0 = no smoothing). default 0.")
    sniffcell_args.add_argument("-i", "--interval", type=int, default=1000, help="Interval for checking methylation changes. default +-1000.")
    sniffcell_args.add_argument("-m", "--min_supporting", type=int, default=3, help="Minimum supporting reads for SV. default 3.")
    sniffcell_args.add_argument("-th", "--threshold", type=float, default=0.2, help="Threshold for differentially methylated regions. default 0.2.")
    sniffcell_args.add_argument("-tf", "--test_function", type=str, default="ttest", help="Statistical test function. default ttest.")

    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(argv)

    return args
