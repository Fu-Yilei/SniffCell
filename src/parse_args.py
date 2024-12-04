#!/usr/bin/env python
import argparse, sys

def parse_args(argv):
    parser = argparse.ArgumentParser(
                    prog='SniffMeth',
                    description='Sniffing CpG methylaiton changed around a (Mosaic) SV',
                    epilog='Version 0.1')
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument(
        "-b",
        "--bam",
        type=str,
        help="Input BAM file.",
        required=True,
    )
    required_args.add_argument(
        "-v",
        "--vcf",
        type=str,
        help="Input SV VCF file (supporting reads needed).",
        required=True,
    )
    required_args.add_argument(
        "-r",
        "--reference",
        type=str,
        help="Reference genome",
        required=True,
    )
    required_args.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output directory",
        required=True,
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads, default 1",
    )

    parser.add_argument(
        "-ob",
        "--output_bam",
        type=bool,
        default=True,
        help="Output SV-related BAM file for each processed SVs. Default: True",
    )

    parser.add_argument(
        "-s",
        "--smoothing",
        type=int,
        default=0,
        help="Enable smoothing, which consider neighboring [input] bps' CpG as an unit. Default: 0 (no smoothing)",
    )

    parser.add_argument(
        "-i",
        "--interval",
        type=int,
        default=1000,
        help="Inverval for checking methylation changes around a SV.",
    )
    
    parser.add_argument(
        "-m",
        "--min_supporting",
        type=int,
        default=3,
        help="Minimum supporting reads requirement for a SV.",
    )

    parser.add_argument(
        "-th",
        "--threshhold",
        type=float,
        default=0.2,
        help="Threshold for determining nearby region is differently methylated or not",
    )

    parser.add_argument(
        "-tf",
        "--test_function",
        type=str,
        default='ttest',
        help="Statistical test function for determining methylation difference. Default: ttest. Other options: ranksum, fisher, chi2",
    )

    parser.add_argument(
        "-b2",
        "--second_bam",
        type=str,
        default=None,
        help="Another BAM file that used for benchmarking.",
    )


    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)

    return args
