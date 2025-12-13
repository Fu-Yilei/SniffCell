import os
import sys
from sniffcell.anno import anno
from sniffcell.find import find
from sniffcell.parse_args import parse_args  # assuming you defined parse_args in args.py
from sniffcell.deconv import deconv  # assuming these modules exist
from sniffcell.dmsv import dmsv

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
        os.makedirs(args.output, exist_ok=True)
        anno.sv_anno(args)
    elif args.command == "deconv":
        deconv.deconv_main(args)
    elif args.command == "dmsv":
        os.makedirs(args.output, exist_ok=True)
        dmsv.dmsv_main(args)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))