from __future__ import annotations

import argparse
import importlib
import json
import shutil
import sys
from pathlib import Path

from sniffcell.discover.discover import (
    DEFAULT_STAGE_ORDER,
    STAGE_ALIASES,
    _parse_stages,
    _required_tools_for_stages,
)


TOOL_DEFAULTS = {
    "bcftools": "bcftools",
    "bgzip": "bgzip",
    "kanpig": "kanpig",
    "medaka": "medaka",
    "modkit": "modkit",
    "sniffles": "sniffles",
    "tabix": "tabix",
    "tdb": "tdb",
    "truvari": "truvari",
    "clair3": "run_clair3.sh",
    "clairs": "run_clairs",
}


def _build_parser(*, prog: str = "sniffcell-check-discover", add_help: bool = True) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Preflight check for external binaries and Python modules used by sniffcell discover.",
        add_help=add_help,
    )
    parser.add_argument(
        "--stages",
        default="all",
        help=(
            "Comma-separated discover stages or aliases. "
            f"Aliases: {', '.join(sorted(STAGE_ALIASES))}. "
            f"Stages: {', '.join(DEFAULT_STAGE_ORDER)}."
        ),
    )
    parser.add_argument("--sniffles-bin", default=None)
    parser.add_argument("--bcftools-bin", default=None)
    parser.add_argument("--bgzip-bin", default=None)
    parser.add_argument("--kanpig-bin", default=None)
    parser.add_argument("--truvari-bin", default=None)
    parser.add_argument("--medaka-bin", default=None)
    parser.add_argument("--tdb-bin", default=None)
    parser.add_argument("--modkit-bin", default=None)
    parser.add_argument("--tabix-bin", default=None)
    parser.add_argument("--clair3-bin", default=None)
    parser.add_argument("--clair3-model-path", default=None)
    parser.add_argument("--clairs-bin", default=None)
    parser.add_argument("--json", action="store_true", default=False, help="Emit machine-readable JSON.")
    return parser


def _resolve_binary(candidate: str | None, default_name: str) -> tuple[bool, str]:
    if candidate:
        path = Path(candidate).expanduser().resolve()
        return path.exists(), str(path)
    found = shutil.which(default_name)
    return bool(found), found or default_name


def _check_python_module(name: str) -> tuple[bool, str]:
    try:
        importlib.import_module(name)
    except Exception as exc:
        return False, str(exc)
    return True, "ok"


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    stages = _parse_stages(args.stages)
    required_tools = _required_tools_for_stages(stages)
    explicit_paths = {
        "sniffles": args.sniffles_bin,
        "bcftools": args.bcftools_bin,
        "bgzip": args.bgzip_bin,
        "kanpig": args.kanpig_bin,
        "truvari": args.truvari_bin,
        "medaka": args.medaka_bin,
        "tdb": args.tdb_bin,
        "modkit": args.modkit_bin,
        "tabix": args.tabix_bin,
        "clair3": args.clair3_bin,
        "clairs": args.clairs_bin,
    }

    checks: list[dict[str, object]] = []
    failures = 0
    warnings = 0

    for tool in sorted(required_tools - {"python"}):
        ok, detail = _resolve_binary(explicit_paths.get(tool), TOOL_DEFAULTS[tool])
        checks.append(
            {
                "kind": "binary",
                "name": tool,
                "required": True,
                "status": "ok" if ok else "missing",
                "detail": detail,
            }
        )
        failures += 0 if ok else 1

    if "tdb_merge" in stages:
        ok, detail = _check_python_module("tdb")
        checks.append(
            {
                "kind": "python_module",
                "name": "tdb",
                "required": True,
                "status": "ok" if ok else "missing",
                "detail": detail,
            }
        )
        failures += 0 if ok else 1

        ok, detail = _check_python_module("seaborn")
        checks.append(
            {
                "kind": "python_module",
                "name": "seaborn",
                "required": False,
                "status": "ok" if ok else "warning",
                "detail": detail if ok else "optional: TR plots are skipped when seaborn is unavailable",
            }
        )
        warnings += 0 if ok else 1

    if "clair3" in stages:
        model_path = Path(args.clair3_model_path).expanduser().resolve() if args.clair3_model_path else None
        ok = bool(model_path and model_path.exists())
        checks.append(
            {
                "kind": "path",
                "name": "clair3_model_path",
                "required": True,
                "status": "ok" if ok else "missing",
                "detail": str(model_path) if model_path else "required when stage set includes clair3",
            }
        )
        failures += 0 if ok else 1

    payload = {
        "stages": list(stages),
        "failures": failures,
        "warnings": warnings,
        "checks": checks,
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"Stages: {', '.join(stages)}")
        for item in checks:
            marker = "OK" if item["status"] == "ok" else ("WARN" if item["status"] == "warning" else "MISSING")
            suffix = "required" if item["required"] else "optional"
            print(f"{marker:8} {item['kind']} {item['name']} ({suffix}) -> {item['detail']}")
        if failures:
            print(f"Result: FAIL ({failures} required dependency checks failed)")
        elif warnings:
            print(f"Result: OK with warnings ({warnings} optional checks failed)")
        else:
            print("Result: OK")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
