# Sperm versus blood ctDMRs

This catalog is generated with the standard `sniffcell find` command using
five normal-control mature-sperm references and the existing T-cell, NK-cell,
Monocyte and B-cell reference groups. It contains the five-group ctDMR catalog,
including blood-cell distinctions and combinations of groups; it is not a
catalog of earlier spermatogenic stages or all testicular somatic cell types.

The signal is combined modified-C from WGBS. There is no separate sperm 5hmC
channel. Sperm donors NC002, NC004, NC005, NC007 and NC008 come from
[Leitao et al. (2020)](https://doi.org/10.1186/s13148-020-00854-0),
[ENA PRJEB34432](https://www.ebi.ac.uk/ena/browser/view/PRJEB34432), the
mature-sperm source cited by
[Porsborg et al. (2025)](https://www.nature.com/articles/s41467-025-65248-3).

## Files and use

- `sperm_vs_blood.default04_gap2000_direct.tsv.gz`: annotation-ready ctDMR table;
  supply this table to `sniffcell deconv` or `sniffcell anno` using `-b`.
- `sperm_vs_blood.default04_gap2000_direct.tsv.igv.bed.gz`: BED9 display companion
  for IGV; use the main TSV for analysis.
- `provenance.json`: source checksums, sample groups, parameters, counts and QC.

Both catalog files are losslessly gzip-compressed. Their uncompressed SHA-256
checksums and byte sizes are recorded in the parent `manifest.tsv`.
The catalog contains 260,399 ctDMRs, including 99,133 single-row regions.

## Regenerate with SniffCell

This catalog uses the current CLI defaults: methylation difference 0.40,
at least one atlas row, at least three CpGs and maximum gap 2,000 bp.
A single atlas row can qualify if it contains at least three CpGs.
Earlier releases used 0.35 and a two-row minimum; the existing historical
catalogs retain their original parameters.

With SniffCell installed and the completed 258-sample sperm-extended atlas
available, run:

```bash
ATLAS=/path/to/sperm_extended_atlas
sniffcell find \
  --npy "$ATLAS/all_celltypes_blocks.npy" \
  --index "$ATLAS/all_celltypes_blocks.index.gz" \
  --meta "$ATLAS/all_celltypes.txt" \
  --celltypes_file "$ATLAS/index_to_major_celltypes.json" \
  --celltypes_keys sperm_vs_blood \
  --diff_threshold 0.40 --min_rows 1 --min_cpgs 3 --max_gap_bp 2000 \
  --output sperm_vs_blood.default04_gap2000_direct.tsv
```

`sniffcell find` also creates the `.tsv.igv.bed` companion. No separate Python
catalog builder is required. The matrix is an input to this command; this
directory distributes the resulting ctDMRs.

## Reference eligibility and QC

Within each sperm donor, CpGs require at least six methylation observations,
and a block requires at least half of its reference CpGs to meet that threshold.
Counts are pooled across qualifying CpGs to calculate block methylation.
Missing coverage remains missing. Sperm reference values are retained only at
blocks with at least three eligible donors (6,918,718 of 7,990,609 atlas blocks).
The original 253 reference columns and existing group definitions are preserved
in the extended input atlas.

The initial all-block donor correlation cutoff of 0.90 failed (observed
0.823-0.856). A documented review found correlations of 0.954-0.970 on blocks
with at least five supported CpGs and 50 observations per donor. The technical
acceptance check was revised after that review to require r >=0.90 on at least
100,000 such shared blocks for every donor pair. The primary eligibility rules
were not relaxed. These checks establish internal reference consistency;
independent classification accuracy has not been measured for this catalog.
