import argparse, sys, os

def parse_args(argv):
    parser = argparse.ArgumentParser(
                    prog='SniffMeth deconv',
                    description='Assing reads to cell types based on fine-graind atlas',
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

    optional_args = parser.add_argument_group("Optional arguments")
    script_dir = os.path.dirname(__file__)
    atlas_celltypes_file_path = os.path.join(script_dir, 'atlas', '39Bisulfite.tsv')

    optional_args.add_argument(
        "-a",
        "--atlas",
        type=str,
        help="Cell type atlas location",
        default=atlas_celltypes_file_path,
    )
    optional_args.add_argument(
        "-c",
        "--tissue",
        type=str,
        help="A file that contains potential cell types, could be: brain_cereb, brain_cortex, whole_blood, etc.",
        default="brain_cereb",
    )   
    optional_args.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads, default 1",
    )
    optional_args.add_argument(
        "-n",
        "--region_number",
        type=int,
        help="Number of regions to be selected",
        default=300,
    )  
    optional_args.add_argument(
        "-m",
        "--method",
        type=str,
        help="Method to select informative regions, could be: std, diff",
        default="diff",
    )
    optional_args.add_argument(
        "-v",
        "--verbose",
        action='store_true',
        help="Verbose mode",
        default=False,
    )
    optional_args.add_argument(
        "-th",
        "--threshold",
        type=float,
        help="Threshold for the minimum confidence of EM algorithm to assign reads to a cell type",
        default=0.9,
    )
    if len(argv) == 1 or 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        args = parser.parse_args(argv)
        return args
