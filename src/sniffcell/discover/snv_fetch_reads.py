"""Enrich harmonized_variants.tsv SNV rows with alt-supporting read names.

Reads harmonized_variants.tsv, and for every row where variant_class == "SNV"
performs a single indexed BAM pileup at the SNV position to collect the names
of reads carrying the alt allele.  Those names are written into the
group_a_read_names / group_b_read_names columns so the downstream anno step
can perform methylation-based cell-type assignment exactly as it does for
TR/SV variants.

Each SNV requires one indexed seek per BAM, so the total overhead is
negligible (seconds for a typical 50-100 SNV run).
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
from pathlib import Path
from typing import Any

import pysam


def _expand_path(value: str | "os.PathLike[str]") -> Path:
    return Path(value).expanduser().resolve()


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        import json as _json
        fh.write(_json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _fetch_alt_read_names(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    alt: str,
) -> list[str]:
    """Return read names at 1-based *pos* that carry the *alt* base.

    Uses stepper="nofilter" so that supplementary/secondary alignments are
    included — consistent with how Clair3 pileup generates its calls.
    """
    names: list[str] = []
    seen: set[str] = set()
    try:
        for col in bam.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            min_base_quality=0,
            min_mapping_quality=0,
            ignore_overlaps=False,
            stepper="nofilter",
        ):
            if col.reference_pos != pos - 1:
                continue
            for pread in col.pileups:
                if pread.is_del or pread.is_refskip:
                    continue
                aln = pread.alignment
                # skip secondary and supplementary alignments — Clair3 counts only
                # primary alignments; including secondary/supplementary would inflate
                # the read list with duplicate read names
                if aln.is_secondary or aln.is_supplementary:
                    continue
                qseq = aln.query_sequence
                if qseq is None:
                    continue
                base = qseq[pread.query_position]
                if base == alt:
                    name = aln.query_name
                    if name and name not in seen:
                        seen.add(name)
                        names.append(name)
    except (ValueError, KeyError):
        # chrom not in BAM — silently skip
        pass
    return names


def enrich_harmonized_with_snv_reads(
    *,
    harmonized_tsv: Path,
    group_a_bam: Path,
    group_b_bam: Path,
    output: Path,
) -> dict[str, Any]:
    """Overwrite group_a/b_read_names for SNV rows in *harmonized_tsv*.

    Non-SNV rows are copied through unchanged.  The output file is written
    atomically (write to final path directly; parent dir is created if absent).

    Returns a summary dict with enrichment counts.
    """
    logger = logging.getLogger("snv_fetch_reads")

    rows: list[dict[str, str]] = []
    fieldnames: list[str] = []
    with harmonized_tsv.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows.append(row)

    bam_a = pysam.AlignmentFile(str(group_a_bam), "rb")
    bam_b = pysam.AlignmentFile(str(group_b_bam), "rb")

    n_enriched = 0
    n_skipped = 0
    try:
        for row in rows:
            if row.get("variant_class") != "SNV":
                continue
            # variant_id format: {chrom}:{pos}:{ref}>{alt}
            vid = row.get("variant_id", "")
            try:
                chrom, pos_str, subst = vid.split(":")
                _ref, alt = subst.split(">")
                pos = int(pos_str)
            except ValueError:
                logger.warning("Cannot parse variant_id %r; skipping", vid)
                n_skipped += 1
                continue

            category = row.get("category", "")
            if category == "group_a_only":
                a_names = _fetch_alt_read_names(bam_a, chrom, pos, alt)
                b_names: list[str] = []
            elif category == "group_b_only":
                a_names = []
                b_names = _fetch_alt_read_names(bam_b, chrom, pos, alt)
            else:
                n_skipped += 1
                continue

            row["group_a_read_names"] = json.dumps(a_names)
            row["group_b_read_names"] = json.dumps(b_names)
            logger.debug(
                "%s category=%s alt=%s a_reads=%d b_reads=%d",
                vid, category, alt, len(a_names), len(b_names),
            )
            n_enriched += 1
    finally:
        bam_a.close()
        bam_b.close()

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    logger.info(
        "Enriched %d SNV rows (%d skipped), total rows=%d → %s",
        n_enriched, n_skipped, len(rows), output,
    )
    return {
        "snv_rows_enriched": n_enriched,
        "snv_rows_skipped": n_skipped,
        "total_rows": len(rows),
        "output": str(output),
    }


def _build_arg_parser(
    *,
    prog: str = "python -m sniffcell.discover.snv_fetch_reads",
    add_help: bool = True,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Enrich harmonized_variants.tsv SNV rows with alt-supporting read "
            "names extracted via BAM pileup.  Non-SNV rows are passed through "
            "unchanged."
        ),
        add_help=add_help,
    )
    parser.add_argument(
        "--harmonized-tsv", required=True,
        help="harmonized_variants.tsv produced by harmonize_variants",
    )
    parser.add_argument(
        "--group-a-bam", required=True,
        help="Deconvolved BAM for group A (e.g. Neuron.bam)",
    )
    parser.add_argument(
        "--group-b-bam", required=True,
        help="Deconvolved BAM for group B (e.g. Oligodendrocyte.bam)",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output path for the enriched harmonized_variants.tsv",
    )
    parser.add_argument(
        "--summary-json", default=None,
        help="Optional path to write a JSON summary of enrichment counts",
    )
    return parser


def snv_fetch_reads_main(cli_args: list[str] | None = None) -> dict[str, Any]:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = parser.parse_args(cli_args)

    harmonized_tsv = _expand_path(args.harmonized_tsv)
    group_a_bam = _expand_path(args.group_a_bam)
    group_b_bam = _expand_path(args.group_b_bam)
    output = _expand_path(args.output)

    for p in (harmonized_tsv, group_a_bam, group_b_bam):
        if not p.exists():
            raise FileNotFoundError(p)

    stats = enrich_harmonized_with_snv_reads(
        harmonized_tsv=harmonized_tsv,
        group_a_bam=group_a_bam,
        group_b_bam=group_b_bam,
        output=output,
    )
    if args.summary_json:
        _write_json(Path(args.summary_json), stats)
    return stats


def main(argv: list[str] | None = None) -> int:
    snv_fetch_reads_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
