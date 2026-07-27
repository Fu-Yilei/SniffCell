import os
import sys

from sniffcell.anno import anno
from sniffcell.find import find
from sniffcell.parse_args import parse_args
from sniffcell.report import report


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)
    if args.command == "find":
        find.find_main(args)
    elif args.command == "anno":
        os.makedirs(args.output, exist_ok=True)
        anno.anno_main(args)
    elif args.command == "report":
        report.report_main(args)
    else:  # pragma: no cover
        raise ValueError(f"Unsupported command: {args.command}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
