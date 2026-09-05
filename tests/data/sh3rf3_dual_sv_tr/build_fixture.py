#!/usr/bin/env python3
"""Build the public SH3RF3 regional SniffCell example fixture.

This is a maintainer utility. It requires the internal source BAM and legacy
Loyfer atlas files. The generated fixture stays on native GRCh38 coordinates
and retains only the nine atlas columns used by the example. Users provide
their own indexed GRCh38 FASTA.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pysam


SOURCE_CHROM = "chr2"
FETCH_START = 109_180_000
FETCH_END = 109_225_000
TARGET_START = 109_199_301
TARGET_END = 109_199_876

NEURON_PREFIXES = tuple(f"GSM56522{i}_Cortex-Neuron-" for i in range(25, 31))
OLIGO_PREFIXES = (
    "GSM5652220_Oligodendrocytes-",
    "GSM5652221_Oligodendrocytes-",
    "GSM5652222_Oligodendrocytes-",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_bam(source_bam: Path, output: Path) -> tuple[int, int]:
    kept = 0
    unique_reads: set[str] = set()
    with pysam.AlignmentFile(str(source_bam), "rb") as source:
        with pysam.AlignmentFile(str(output), "wb", header=source.header) as target:
            for record in source.fetch(SOURCE_CHROM, FETCH_START, FETCH_END):
                target.write(record)
                kept += 1
                unique_reads.add(record.query_name)
    pysam.index(str(output))
    return kept, len(unique_reads)


def build_atlas(
    source_npy: Path,
    source_index: Path,
    source_meta: Path,
    source_celltypes: Path,
    output_dir: Path,
) -> tuple[int, int]:
    samples = [line.strip() for line in source_meta.read_text().splitlines() if line.strip()]
    selected_columns = [
        index
        for index, sample in enumerate(samples)
        if sample.startswith(NEURON_PREFIXES) or sample.startswith(OLIGO_PREFIXES)
    ]
    selected_rows: list[int] = []
    shifted_index: list[str] = []
    open_index = gzip.open if source_index.suffix == ".gz" else open
    with open_index(source_index, "rt", encoding="utf-8") as handle:
        for row_index, line in enumerate(handle):
            chrom, start_text, end_text, start_cpg, end_cpg = line.rstrip().split("\t")
            start = int(start_text)
            end = int(end_text)
            if chrom != SOURCE_CHROM or start >= FETCH_END or end <= FETCH_START:
                continue
            selected_rows.append(row_index)
            shifted_index.append(
                f"{SOURCE_CHROM}\t{start}\t{end}"
                f"\t{start_cpg}\t{end_cpg}\n"
            )

    matrix = np.load(source_npy, mmap_mode="r")
    subset = np.asarray(matrix[selected_rows][:, selected_columns], dtype=np.float32)
    np.save(output_dir / "atlas.npy", subset)
    (output_dir / "atlas.index.tsv").write_text("".join(shifted_index), encoding="utf-8")
    (output_dir / "atlas.samples.txt").write_text(
        "".join(samples[index] + "\n" for index in selected_columns), encoding="utf-8"
    )

    source_mapping = json.loads(source_celltypes.read_text(encoding="utf-8"))
    mapping = {"brain_cereb": source_mapping["brain_cereb"]}
    (output_dir / "celltypes.json").write_text(
        json.dumps(mapping, indent=2) + "\n", encoding="utf-8"
    )
    return len(selected_rows), len(selected_columns)


def build_tr_catalog(source_tr_bed: Path, output: Path) -> int:
    rows: list[str] = []
    with source_tr_bed.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 4 or fields[0] != SOURCE_CHROM:
                continue
            start = int(fields[1])
            end = int(fields[2])
            if start >= FETCH_END or end <= FETCH_START:
                continue
            rows.append("\t".join(fields) + "\n")
    output.write_text("".join(rows), encoding="utf-8")
    return len(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-bam", type=Path, required=True)
    parser.add_argument("--source-atlas-npy", type=Path, required=True)
    parser.add_argument("--source-atlas-index", type=Path, required=True)
    parser.add_argument("--source-atlas-meta", type=Path, required=True)
    parser.add_argument("--source-celltypes", type=Path, required=True)
    parser.add_argument("--source-tr-bed", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "inputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bam = args.output_dir / "SH3RF3_example.bam"
    alignment_count, read_count = build_bam(args.source_bam, bam)
    atlas_rows, atlas_columns = build_atlas(
        args.source_atlas_npy,
        args.source_atlas_index,
        args.source_atlas_meta,
        args.source_celltypes,
        args.output_dir,
    )
    tr_count = build_tr_catalog(args.source_tr_bed, args.output_dir / "tr_catalog.bed")
    (args.output_dir / "target.bed").write_text(
        f"{SOURCE_CHROM}\t{TARGET_START}\t{TARGET_END}\tSH3RF3_AATGG\n",
        encoding="utf-8",
    )

    outputs = sorted(
        path
        for path in args.output_dir.iterdir()
        if path.is_file() and path.name != "fixture_manifest.json"
    )
    manifest = {
        "fixture": "SH3RF3 AATGG repeat expansion detected by SV and TR workflows",
        "source_coordinates_grch38": f"{SOURCE_CHROM}:{TARGET_START}-{TARGET_END}",
        "fixture_coordinates_grch38": f"{SOURCE_CHROM}:{TARGET_START}-{TARGET_END}",
        "coordinate_shift": 0,
        "fetch_window_grch38": f"{SOURCE_CHROM}:{FETCH_START}-{FETCH_END}",
        "alignment_records": alignment_count,
        "unique_read_names": read_count,
        "atlas_rows": atlas_rows,
        "atlas_samples": atlas_columns,
        "tr_catalog_rows": tr_count,
        "reference": "User-supplied indexed GRCh38 no-alt FASTA; no reference sequence is bundled.",
        "sha256": {path.name: sha256(path) for path in outputs},
    }
    (args.output_dir / "fixture_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
