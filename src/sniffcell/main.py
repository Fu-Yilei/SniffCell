import os
import sys
from sniffcell.anno import anno
from sniffcell.find import find
from sniffcell.parse_args import parse_args  # assuming you defined parse_args in args.py
from sniffcell.deconv import deconv  # assuming these modules exist
from sniffcell.dmsv import dmsv
from sniffcell.viz import viz
from sniffcell.viz import igvviz
from sniffcell.report import report_main

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)
    if args.command == "find":
        find.find_main(args)
    elif args.command == "anno":
        os.makedirs(args.output, exist_ok=True)
        anno.anno_main(args)
    elif args.command == "svanno":
        output_dir = os.path.dirname(os.path.abspath(args.output))
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        anno.sv_anno(args)
    elif args.command == "deconv":
        deconv.deconv_main(args)
    elif args.command == "dmsv":
        os.makedirs(args.output, exist_ok=True)
        dmsv.dmsv_main(args)
    elif args.command == "viz":
        viz.viz_main(args)
    elif args.command == "igvviz":
        igvviz.igvviz_main(args)
    elif args.command == "report":
        report_main(args)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
