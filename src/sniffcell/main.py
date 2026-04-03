import os
import sys
from sniffcell.anno import anno
from sniffcell.discover import envcheck as discover_envcheck
from sniffcell.discover import harmonize_variants as discover_harmonize_variants
from sniffcell.discover import snv_post_processing as discover_snv_post_processing
from sniffcell.discover import sv_post_processing as discover_sv_post_processing
from sniffcell.discover import tr_post_processing as discover_tr_post_processing
from sniffcell.find import find
from sniffcell.parse_args import parse_args  # assuming you defined parse_args in args.py
from sniffcell.deconv import deconv  # assuming these modules exist
from sniffcell.dmsv import dmsv
from sniffcell.discover import discover_main
from sniffcell import sv_discovery as discover_sv_discovery
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
        anno.svanno_main(args)
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
    elif args.command == "discover":
        discover_section = getattr(args, "discover_section", None)
        if discover_section is None:
            discover_main(args)
        elif discover_section == "tools":
            discover_tools_command = getattr(args, "discover_tools_command", None)
            if discover_tools_command == "run":
                discover_main(args)
            elif discover_tools_command == "sv":
                return discover_sv_discovery.main(argv[3:])
            elif discover_tools_command == "check":
                return discover_envcheck.main(argv[3:])
            else:  # pragma: no cover
                raise ValueError(f"Unsupported discover tools command: {discover_tools_command}")
        elif discover_section == "ctprocessing":
            discover_ct_command = getattr(args, "discover_ctprocessing_command", None)
            if discover_ct_command == "snv":
                return discover_snv_post_processing.main(argv[3:])
            if discover_ct_command == "sv":
                return discover_sv_post_processing.main(argv[3:])
            if discover_ct_command == "tr":
                return discover_tr_post_processing.main(argv[3:])
            if discover_ct_command == "harmonize":
                return discover_harmonize_variants.main(argv[3:])
            raise ValueError(f"Unsupported discover ctprocessing command: {discover_ct_command}")  # pragma: no cover
        else:  # pragma: no cover
            raise ValueError(f"Unsupported discover command group: {discover_section}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
