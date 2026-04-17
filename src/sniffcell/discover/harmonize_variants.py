from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any


_HARMONIZED_COLS: list[str] = [
    "chrom", "start", "end",
    "variant_class",       # TR | SV
    "variant_id",          # trid (TR) | sv_id (SV)
    "variant_subtype",     # expansion_hap1 / contraction_hap2 (TR) | DEL/INS/DUP/… (SV)
    "category",            # group_a_only | group_b_only | shared
    "change_size_bp",      # change_length_bp (TR) | abs(sv_len) (SV)
    "group_a_alt_reads",   # alt-supporting read count for group A
    "group_b_alt_reads",   # alt-supporting read count for group B
    "group_a_read_names",  # JSON list
    "group_b_read_names",  # JSON list
]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _truthy_or_none(value: object) -> bool | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def _load_tr_rows(
    tr_bed: Path | None,
    group_a_label: str,
    group_b_label: str,
    *,
    include_all_tr_candidates: bool,
) -> list[dict[str, Any]]:
    """Convert tr_changes.bed.tsv rows to the harmonized schema.

    The TR file records which group carries the changed allele (change_group)
    and which is the baseline (baseline_group).  The mapping to group_a /
    group_b is resolved via the supplied labels so read counts and read-name
    lists land in the right columns.
    """
    rows: list[dict[str, Any]] = []
    if tr_bed is None or not tr_bed.exists():
        return rows
    for r in _read_tsv(tr_bed):
        pass_flag = _truthy_or_none(r.get("tr_pass_for_harmonized"))
        if not include_all_tr_candidates and pass_flag is False:
            continue
        change_names = r.get("change_support_read_names", r.get("change_read_names", "[]"))
        change_alt = r.get("n_change_support_reads", r.get("n_change_reads", "."))
        change_group = r.get("change_group", "")
        if change_group == group_a_label:
            category = "group_a_only"
            a_alt = change_alt
            b_alt = 0
            a_names = change_names
            b_names = "[]"
        elif change_group == group_b_label:
            category = "group_b_only"
            a_alt = 0
            b_alt = change_alt
            a_names = "[]"
            b_names = change_names
        else:
            category = "unknown"
            a_alt = b_alt = 0
            a_names = b_names = "[]"

        change_type = r.get("change_type", ".")
        change_allele = r.get("change_allele", ".")
        subtype = f"{change_type}_{change_allele}"

        rows.append({
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "variant_class": "TR",
            "variant_id": r.get("trid", "."),
            "variant_subtype": subtype,
            "category": category,
            "change_size_bp": r.get("change_length_bp", "."),
            "group_a_alt_reads": a_alt,
            "group_b_alt_reads": b_alt,
            "group_a_read_names": a_names,
            "group_b_read_names": b_names,
        })
    return rows


def _load_sv_rows(sv_bed: Path) -> list[dict[str, Any]]:
    """Convert sv_changes.bed.tsv rows to the harmonized schema."""
    rows: list[dict[str, Any]] = []
    if not sv_bed.exists():
        return rows
    for r in _read_tsv(sv_bed):
        sv_len_str = r.get("sv_len", ".")
        change_size: int | str = "."
        if sv_len_str not in (".", ""):
            try:
                change_size = abs(int(sv_len_str))
            except ValueError:
                pass

        rows.append({
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "variant_class": "SV",
            "variant_id": r.get("sv_id", "."),
            "variant_subtype": r.get("sv_type", "."),
            "category": r.get("category", "."),
            "change_size_bp": change_size,
            "group_a_alt_reads": r.get("group_a_ad_alt", "."),
            "group_b_alt_reads": r.get("group_b_ad_alt", "."),
            "group_a_read_names": r.get("group_a_read_names", "[]"),
            "group_b_read_names": r.get("group_b_read_names", "[]"),
        })
    return rows


def _load_snv_rows(
    snv_tsv: Path | None,
    group_a_label: str,
    group_b_label: str,
    *,
    include_all_snv_candidates: bool,
) -> list[dict[str, Any]]:
    """Convert snv_changes.tsv rows to the harmonized schema.

    Read names are read from the group_a/b_read_names columns written by
    snv_post_processing (populated via BAM pileup there).  Falls back to
    empty ``[]`` if those columns are absent (e.g. older runs).
    """
    rows: list[dict[str, Any]] = []
    if snv_tsv is None or not snv_tsv.exists():
        return rows
    for r in _read_tsv(snv_tsv):
        if not include_all_snv_candidates:
            pass_flag = r.get("snv_pass_for_harmonized", "").strip().lower()
            if pass_flag not in {"true", "1", "yes"}:
                continue
        direction = r.get("direction", "")
        chrom = r["chrom"]
        pos = int(r["pos"])
        ref = r.get("ref", ".")
        alt = r.get("alt", ".")
        alt_ad = r.get("target_alt_ad", "0")
        if direction == "group_a_only":
            category = "group_a_only"
            a_alt, b_alt = alt_ad, "0"
        elif direction == "group_b_only":
            category = "group_b_only"
            a_alt, b_alt = "0", alt_ad
        else:
            category = direction
            a_alt = b_alt = "0"
        rows.append({
            "chrom": chrom,
            "start": pos - 1,
            "end": pos,
            "variant_class": "SNV",
            "variant_id": f"{chrom}:{pos}:{ref}>{alt}",
            "variant_subtype": f"{ref}>{alt}",
            "category": category,
            "change_size_bp": 1,
            "group_a_alt_reads": a_alt,
            "group_b_alt_reads": b_alt,
            "group_a_read_names": r.get("group_a_read_names", "[]") or "[]",
            "group_b_read_names": r.get("group_b_read_names", "[]") or "[]",
        })
    return rows


_CHROM_ORDER: dict[str, int] = {f"chr{i}": i for i in range(1, 23)}
_CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})


def _sort_key(row: dict[str, Any]) -> tuple[int, str, int, int]:
    chrom = str(row["chrom"])
    chrom_rank = _CHROM_ORDER.get(chrom, 99)
    return (chrom_rank, chrom, int(row["start"]), int(row["end"]))


def write_harmonized_variants(
    *,
    output: Path,
    group_a_label: str,
    group_b_label: str,
    tr_bed: Path | None = None,
    sv_bed: Path | None = None,
    snv_tsv: Path | None = None,
    include_all_tr_candidates: bool = True,
    include_all_snv_candidates: bool = False,
) -> dict[str, int]:
    tr_rows = _load_tr_rows(
        tr_bed,
        group_a_label,
        group_b_label,
        include_all_tr_candidates=include_all_tr_candidates,
    )
    sv_rows = _load_sv_rows(sv_bed) if sv_bed is not None else []
    snv_rows = _load_snv_rows(
        snv_tsv,
        group_a_label,
        group_b_label,
        include_all_snv_candidates=include_all_snv_candidates,
    )

    all_rows = tr_rows + sv_rows + snv_rows
    all_rows.sort(key=_sort_key)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=_HARMONIZED_COLS)
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    return {
        "total_rows": len(all_rows),
        "tr_rows": len(tr_rows),
        "sv_rows": len(sv_rows),
        "snv_rows": len(snv_rows),
    }


def _build_arg_parser(
    *,
    prog: str = "python -m sniffcell.discover.harmonize_variants",
    add_help: bool = True,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Combine TR (tr_changes.bed.tsv) and SV (sv_changes.bed.tsv) findings "
            "into a single unified BED-like TSV sorted by genomic position."
        ),
        add_help=add_help,
    )
    parser.add_argument("--tr-bed", required=True, help="tr_changes.bed.tsv from tr_post_processing")
    parser.add_argument("--sv-bed", required=True, help="sv_changes.bed.tsv from sv_post_processing")
    parser.add_argument("--snv-tsv", default=None, help="snv_changes.tsv from snv_post_processing (optional)")
    parser.add_argument(
        "--group-a-label", required=True,
        help="Sample label used for group A in the TR bed (values in the change_group / baseline_group columns)",
    )
    parser.add_argument(
        "--group-b-label", required=True,
        help="Sample label used for group B in the TR bed",
    )
    parser.add_argument("--output", required=True, help="Output harmonized TSV path")
    parser.add_argument(
        "--passing-tr-only",
        action="store_true",
        help="Keep only TR rows marked tr_pass_for_harmonized=true. Default keeps all TR candidates.",
    )
    return parser


def harmonize_main(cli_args=None) -> None:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = parser.parse_args(cli_args)

    tr_bed = Path(args.tr_bed)
    sv_bed = Path(args.sv_bed)
    snv_tsv = Path(args.snv_tsv) if args.snv_tsv else None
    output = Path(args.output)
    stats = write_harmonized_variants(
        output=output,
        group_a_label=args.group_a_label,
        group_b_label=args.group_b_label,
        tr_bed=tr_bed,
        sv_bed=sv_bed,
        snv_tsv=snv_tsv,
        include_all_tr_candidates=not args.passing_tr_only,
    )

    logging.info(
        "Wrote %d rows (%d TR, %d SV, %d SNV) → %s",
        stats["total_rows"], stats["tr_rows"], stats["sv_rows"], stats["snv_rows"], output,
    )


def main(argv: list[str] | None = None) -> int:
    harmonize_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
