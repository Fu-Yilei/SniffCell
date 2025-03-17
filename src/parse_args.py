#!/usr/bin/env python
import argparse
import sys
import os

def parse_args(argv):
    parser = argparse.ArgumentParser(
        prog="SniffMeth",
        description="Sniffing CpG methylation changes around a (Mosaic) SV",
        epilog="Version 0.2",
    )
    
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads, default 1.")


    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-b", "--bam", type=str, required=True, help="Input BAM file.")
    required_args.add_argument("-v", "--vcf", type=str, required=True, help="Input SV VCF file (supporting reads needed).")
    required_args.add_argument("-r", "--reference", type=str, required=True, help="Reference genome.")
    required_args.add_argument("-o", "--output", type=str, required=True, default=f"{os.path.abspath("./")}", help=f"Output directory. Default: {os.path.abspath('./')}")
    required_args.add_argument("-d", "--deconv", action="store_true", help="Run deconvolution if specified.")

    optional_args = parser.add_argument_group("SniffMeth deconvolution optional arguments")
    script_dir = os.path.dirname(__file__)
    atlas_celltypes_file_path = os.path.join(script_dir, "atlas", "39Bisulfite.tsv")
    json_file_file_path = os.path.join(script_dir, "atlas", "tissue_celltypes.json")
    # optional_args.add_argument("-od", "--deconv_output", type=str, help="Output directory for deconvolution results. default: <output>/deconv_output")
    optional_args.add_argument("-a", "--atlas", type=str, default=atlas_celltypes_file_path, help=f"Cell type atlas location. Default: {atlas_celltypes_file_path}")
    optional_args.add_argument("-c", "--tissue", type=str, default="brain_cereb", help=f"JSON file with potential cell types. Default: {json_file_file_path}, brain_cereb.")
    optional_args.add_argument("-n", "--region_number", type=int, default=300, help="Number of regions to be selected, default 300.")
    optional_args.add_argument("-me", "--method", type=str, default="diff", help="Region selection method: std or diff, default diff. Diff selects the regions with certain cell type as low methylation while all other cell types have high methylation. std selects regions with highest methylation value std in all cell types.")
    optional_args.add_argument("-vb", "--verbose", action="store_true", help="Enable verbose mode.")
    optional_args.add_argument("-conf", "--confidence", type=float, default=0.9, help="Minimum confidence threshold for EM algorithm. Default 0.9.")

    sniffmeth_args = parser.add_argument_group("SniffMeth optional arguments")
    sniffmeth_args.add_argument("-ob", "--output_bam", action="store_true", help="Output SV-related BAM file.")
    sniffmeth_args.add_argument("-p", "--primary_only", action="store_true", help="Keep only primary alignments in BAM.")
    sniffmeth_args.add_argument("-s", "--smoothing", type=int, default=0, help="Enable smoothing (0 = no smoothing). default 0.")
    sniffmeth_args.add_argument("-i", "--interval", type=int, default=1000, help="Interval for checking methylation changes. default +-1000.")
    sniffmeth_args.add_argument("-m", "--min_supporting", type=int, default=3, help="Minimum supporting reads for SV. default 3.")
    sniffmeth_args.add_argument("-th", "--threshold", type=float, default=0.2, help="Threshold for differentially methylated regions. default 0.2.")
    sniffmeth_args.add_argument("-tf", "--test_function", type=str, default="ttest", help="Statistical test function. default ttest.")
    sniffmeth_args.add_argument("-b2", "--second_bam", type=str, help="Benchmark BAM file. If provided, benchmarking will be performed.")

    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args(argv)

    return args
