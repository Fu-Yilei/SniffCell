#!/usr/bin/env python3
"""Validate the expected native-GRCh38 SH3RF3 example results."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


TARGET_START = 109_199_301
TARGET_END = 109_199_876


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise AssertionError(f"Missing output: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {Path(sys.argv[0]).name} OUTPUT_DIR")
    output = Path(sys.argv[1])
    discover = (
        output
        / "deconv/deconv_requested_group_splits/discover/dual_sv_tr"
        / "harmonized_variants.tsv"
    )
    variants = read_tsv(discover)

    tr = next(
        row
        for row in variants
        if row["variant_class"] == "TR"
        and row["variant_id"] == "chr2_109199301_109199876"
    )
    assert tr["variant_subtype"] == "expansion_all"
    assert tr["category"] == "group_a_only"
    assert int(tr["change_size_bp"]) == 2502
    assert int(tr["group_a_alt_reads"]) == 5
    assert int(tr["group_b_alt_reads"]) == 0

    overlapping_insertions = [
        row
        for row in variants
        if row["variant_class"] == "SV"
        and row["variant_subtype"] == "INS"
        and row["chrom"] == "chr2"
        and int(row["start"]) < TARGET_END
        and int(row["end"]) > TARGET_START
    ]
    assert len(overlapping_insertions) == 1, "Expected one SH3RF3 insertion call"
    insertion = overlapping_insertions[0]
    assert insertion["category"] == "group_a_only"
    assert int(insertion["change_size_bp"]) == 726
    assert int(insertion["group_a_alt_reads"]) == 3
    assert int(insertion["group_b_alt_reads"]) == 0

    assignments = read_tsv(output / "anno/variant_assignment_readable.tsv")
    assignment = next(row for row in assignments if row["id"] == tr["variant_id"])
    assert assignment["classified_celltypes"] == "Neuron"
    assert int(assignment["n_supporting"]) == 5
    assert int(assignment["n_overlapped"]) == 5
    assert float(assignment["majority_pct"]) == 1.0

    insertion_assignment = next(
        row for row in assignments if row["id"] == insertion["variant_id"]
    )
    assert insertion_assignment["classified_celltypes"] == "Neuron"
    assert int(insertion_assignment["n_supporting"]) == 3
    assert int(insertion_assignment["n_overlapped"]) == 3
    assert float(insertion_assignment["majority_pct"]) == 1.0

    splits = read_tsv(output / "deconv/deconv_requested_group_splits/requested_group_splits.tsv")
    split_counts = {row["requested_group"]: int(row["n_reads"]) for row in splits}
    assert split_counts == {"Neuron": 182, "Oligodendrocyte": 356}

    assert (output / "anno/SH3RF3_TR_expansion.png").is_file()
    assert (output / "anno/report/index.html").is_file()
    print(
        "PASS: native GRCh38 SH3RF3 example produced a Neuron-only 2,502-bp "
        "TR expansion and a Neuron-only 726-bp insertion."
    )


if __name__ == "__main__":
    main()
