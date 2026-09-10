from __future__ import annotations

import csv
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

import pysam


_FIELDS = {
    "SC_ID": ("String", "Original harmonized variant ID (percent-encoded)"),
    "SC_ROW": ("Integer", "One-based data row in harmonized_variants.tsv"),
    "SC_CLASS": ("String", "SniffCell variant class"),
    "SC_SUBTYPE": ("String", "SniffCell variant subtype (percent-encoded)"),
    "SC_CATEGORY": ("String", "Split-BAM comparison category; not methylation-confirmed assignment"),
    "SC_GROUP_A": ("String", "Comparison group A label (percent-encoded; composite groups remain composite)"),
    "SC_GROUP_B": ("String", "Comparison group B label (percent-encoded; composite groups remain composite)"),
    "SC_TARGET_GROUP": ("String", "Group carrying the candidate change (percent-encoded); absent for shared or unknown"),
    "SC_SUPPORT_A": ("Integer", "Group A support from harmonized TSV; zero does not establish reference genotype or coverage"),
    "SC_SUPPORT_B": ("Integer", "Group B support from harmonized TSV; zero does not establish reference genotype or coverage"),
    "SC_CHANGE_BP": ("Integer", "SniffCell change size; TR values compare groups, not the reference or a specific ALT allele"),
}


def _tr_catalog(path: Path | None) -> dict[str, tuple[str, int, int]]:
    loci = {}
    if path is None or not path.exists():
        return loci
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            locus = (fields[0], int(fields[1]), int(fields[2]))
            for token in fields[3].split(";"):
                if token.startswith("ID="):
                    identifier = token[3:]
                    if identifier in loci and loci[identifier] != locus:
                        raise ValueError(f"Ambiguous repeat catalog ID: {identifier}")
                    loci[identifier] = locus
    return loci


def _row_key(row: dict[str, str]) -> tuple:
    variant_class = row["variant_class"]
    if variant_class == "TR":
        return row["chrom"], int(row["start"]), int(row["end"])
    if variant_class == "SNV":
        ref, alt = row["variant_subtype"].split(">")
        return row["chrom"], int(row["start"]) + 1, ref, alt
    return row["chrom"], int(row["start"]) + 1, row["variant_id"], row["variant_subtype"]


def _record_keys(record, variant_class: str, loci: dict) -> list[tuple]:
    if variant_class == "TR":
        identifier = record.info.get("TRID")
        if identifier in loci:
            return [loci[identifier]]
        return []
    if variant_class == "SNV":
        return [(record.contig, record.pos, record.ref, alt) for alt in record.alts or ()]
    return [(record.contig, record.pos, record.id or f"{record.contig}:{record.pos}", record.info.get("SVTYPE", "."))]


def _annotate(record, row: dict[str, str], row_number: int, groups: tuple[str, str]) -> None:
    values = {
        "SC_ID": row["variant_id"], "SC_ROW": row_number,
        "SC_CLASS": row["variant_class"], "SC_SUBTYPE": row["variant_subtype"],
        "SC_CATEGORY": row["category"], "SC_GROUP_A": groups[0], "SC_GROUP_B": groups[1],
        "SC_SUPPORT_A": row["group_a_alt_reads"], "SC_SUPPORT_B": row["group_b_alt_reads"],
        "SC_CHANGE_BP": row["change_size_bp"],
    }
    target = {"group_a_only": groups[0], "group_b_only": groups[1]}.get(row["category"])
    if target is not None:
        values["SC_TARGET_GROUP"] = target
    for field, value in values.items():
        if value in ("", ".", None):
            continue
        record.info[field] = int(value) if _FIELDS[field][0] == "Integer" else quote(str(value), safe="-_.~")


def export_harmonized_vcfs(
    harmonized_tsv: Path,
    *,
    sources: dict[str, tuple[str, Path | None]],
    group_a: str,
    group_b: str,
    tr_bed: Path | None = None,
) -> dict:
    """Subset caller records and add discovery evidence without reconstructing alleles.

    Sources map output suffixes to (variant class, caller VCF). TR matching uses
    the catalog's TRID-to-interval mapping, independent of caller padding.
    Each exported record represents one harmonized row; SC_ROW disambiguates
    repeated locus IDs. Caller genotypes are retained, not somatic genotypes.
    """
    with harmonized_tsv.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    loci = _tr_catalog(tr_bed)
    results = {}
    for suffix, (variant_class, source) in sources.items():
        selected = [(number, row) for number, row in enumerate(rows, 1) if row["variant_class"] == variant_class]
        output = harmonized_tsv.with_name(f"{harmonized_tsv.stem}.{suffix}.vcf.gz")
        result = {"source": str(source) if source else None, "output": None,
                  "rows": len(selected), "exported_records": 0, "unmatched_rows": []}
        results[suffix] = result
        if source is None or not source.exists():
            result["status"] = "missing_source"
            result["unmatched_rows"] = [number for number, _ in selected]
            output.unlink(missing_ok=True)
            Path(str(output) + ".tbi").unlink(missing_ok=True)
            if selected:
                logging.warning("VCF export %s: source unavailable for %d harmonized rows", suffix, len(selected))
            continue
        lookup = defaultdict(list)
        for number, row in selected:
            lookup[_row_key(row)].append((number, row))
        records = []
        matched = set()
        with pysam.VariantFile(str(source)) as caller:
            header = caller.header.copy()
            for field, (value_type, description) in _FIELDS.items():
                if field in header.info:
                    raise ValueError(f"Source VCF already defines reserved export field {field}: {source}")
                header.info.add(field, 1, value_type, description)
            header.add_meta("SniffCellEvidence", "Split-BAM group-specific candidates, not anno-confirmed cell types; TR evidence is locus-level, not ALT-specific")
            for record in caller if lookup else ():
                for key in _record_keys(record, variant_class, loci):
                    for number, row in lookup.get(key, ()):
                        if number in matched:
                            raise ValueError(f"Ambiguous caller records for harmonized row {number}: {source}")
                        exported = record.copy()
                        exported.translate(header)
                        _annotate(exported, row, number, (group_a, group_b))
                        records.append(exported)
                        matched.add(number)
            contigs = {name: rank for rank, name in enumerate(header.contigs)}
            records.sort(key=lambda record: (contigs[record.contig], record.start, record.stop, record.info["SC_ROW"]))
            temporary = output.with_name(output.name + ".tmp.vcf.gz")
            try:
                with pysam.VariantFile(str(temporary), "wz", header=header) as writer:
                    for record in records:
                        writer.write(record)
                pysam.tabix_index(str(temporary), preset="vcf", force=True)
                temporary.replace(output)
                Path(str(temporary) + ".tbi").replace(Path(str(output) + ".tbi"))
            finally:
                temporary.unlink(missing_ok=True)
                Path(str(temporary) + ".tbi").unlink(missing_ok=True)
        result.update(output=str(output), exported_records=len(records), status="complete",
                      unmatched_rows=[number for number, _ in selected if number not in matched])
        if result["unmatched_rows"]:
            result["status"] = "partial"
            logging.warning("VCF export %s: %d/%d harmonized rows unmatched", suffix, len(result["unmatched_rows"]), len(selected))
    return results
