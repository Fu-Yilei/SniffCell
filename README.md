# SniffCell Lite

SniffCell Lite keeps only two commands:

- `sniffcell-lite find`: call ctDMR catalogs from a methylation atlas.
- `sniffcell-lite anno`: annotate variants from supporting reads, a BAM, a reference FASTA, and a ctDMR catalog.

This branch intentionally keeps only the lite `find` and `anno` workflow.

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

## Tests

```bash
pytest -q
```
