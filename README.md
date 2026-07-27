# SniffCell Lite

SniffCell Lite keeps a compact workflow:

- `sniffcell-lite find`: call ctDMR catalogs from a methylation atlas.
- `sniffcell-lite anno`: annotate variants from supporting reads, a BAM, a reference FASTA, and a ctDMR catalog.
- `sniffcell-lite report`: build a lite-native HTML report, optionally with per-variant figures.

This branch intentionally keeps only the lite `find`, `anno`, and reporting workflow.

## Install

```bash
pip install .
```

## Find

`find` keeps the original atlas-driven flavor. The bundled tissue atlas supports tissue code or tissue name lookup through `-ck`.

```bash
sniffcell-lite find \
  -n atlas/all_celltypes_blocks.npy \
  -i atlas/all_celltypes_blocks.index.gz \
  -m atlas/all_celltypes.txt \
  -ck "Colon, Ascending" \
  -o colon_ascending.ctdmr.tsv
```

Equivalent code form:

```bash
sniffcell-lite find -ck 3E -o colon_ascending.ctdmr.tsv
```

The default cell-type JSON is the packaged `sniffcell/data/tissue_atlas.json`. You can still pass a custom atlas with `-cf`.

## Anno

Single-variant mode requires a BAM and reference FASTA. SniffCell Lite maps the variant-supporting reads in the BAM, selects ctDMRs from the catalog that overlap those supporting-read alignment spans, computes methylation from the BAM only for those reads at those ctDMRs, and assigns the variant from the targeted read-level methylation calls. ctDMR evidence is not capped to a fixed distance from the variant.

```bash
sniffcell-lite anno \
  -i sample.bam \
  -r ref.fa \
  --variant-name variant_001 \
  --variant-location chr1:100000-101000 \
  --supporting-reads readA,readB,readC \
  --catalog colon_ascending.ctdmr.tsv \
  -o anno_out
```

Batch mode uses a TSV or CSV with these columns:

```text
variant_name    variant_location    supporting_reads    catalog    bam    reference
```

Run:

```bash
sniffcell-lite anno --batch variants.tsv -o anno_out
```

`supporting_reads` accepts comma, pipe, semicolon, whitespace-delimited text, JSON list text, or `@path/to/read_names.txt`.

If every batch row uses the same reference, pass it once instead of adding a `reference` column:

```bash
sniffcell-lite anno --batch variants.tsv -r ref.fa -o anno_out
```

## Outputs

`sniffcell-lite find` writes:

- `*.tsv`: ctDMR catalog
- `*.tsv.igv.bed`: IGV BED companion
- `*.tsv.catalog.json`: catalog manifest when tissue metadata was used

`sniffcell-lite anno` writes:

- `variant_assignment.tsv`
- `variant_assignment_readable.tsv`
- `variant_assignment_readable_long.tsv`
- `reads_classification.tsv`
- `support_read_mappings.tsv`
- `anno_compact_manifest.json` or `anno_batch_manifest.json`

## Report

`report` reads a `sniffcell-lite anno` output folder directly. Batch reports keep
the per-row BAM and ctDMR catalog from the original lite batch manifest, so
multi-tissue/multi-BAM annotation outputs do not need to be coerced into a
single full-SniffCell manifest.

Table-only report:

```bash
sniffcell-lite report --anno_output anno_out -o anno_report
```

Report with lightweight per-variant figures:

```bash
sniffcell-lite report --anno_output anno_out --with_figures -o anno_report
```

Useful filters:

```bash
sniffcell-lite report \
  --anno_output anno_out \
  --min_overlap_pct 0 \
  --min_majority_pct 0 \
  --max_variants 100 \
  --with_figures \
  -o anno_report
```

To report only listed variants, or to remove listed variants, provide a TSV/CSV
with an `id` column or a one-ID-per-line file:

```bash
sniffcell-lite report \
  --anno_output anno_out \
  --include_variants selected_ids.tsv \
  --with_figures \
  -o selected_report

sniffcell-lite report \
  --anno_output anno_out \
  --exclude_variants rejected_ids.tsv \
  --with_figures \
  -o filtered_report
```

`--include_variants` and `--exclude_variants` are mutually exclusive. The report
writes `index.html`, `high_confidence_variants.tsv`, `report_manifest.json`,
normalized inclusion/exclusion provenance when used, and optional files under
`figures/`.

## Tests

```bash
pytest -q
```
