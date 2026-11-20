# Precomputed ctDMR atlases

This directory contains the supplied `sniffcell find` ctDMR tables and BED files
for brain universal, PBMC, breast, colon, kidney, liver, lung, pancreas, and
[sperm versus blood](sperm_vs_blood/README.md).
The [brain 5hmC catalog](brain_5hmc/README.md) provides separate modifiedC,
5mC and 5hmC scoring rows.
The [whole-blood 5hmC catalog](blood_5hmc/README.md) adds granulocytes and
two-donor FANS PBMC m/h markers in one shared five-class schema.
The [SMaHT production release](smaht_20261016/README.md) maps tissue codes to
the exact scoring catalogs used for the completed 407-sample run, including
expanded lung and segment-specific colon catalogs. It retains the existing
brain 5hmC file and documents the differences from the original tissue files.
The [experimental brain annotation catalog](brain_anno/README.md) supplies
nine-class modifiedC profiles for `anno`, while retaining the original catalogs
for deconvolution and read splitting.

Both `combined_dt04` and `default035_gap2000_direct` file sets are preserved
as supplied, including their tissue-specific BED filenames.

The original 32 source files and two sperm-versus-blood catalog files are stored
with lossless gzip compression. The sperm catalog uses the updated CLI defaults
(0.40 methylation difference, one-row minimum, three-CpG minimum, 2 kb gap); its directory documents the exact
`sniffcell find` command and reference provenance. To extract an
individual file while retaining the compressed copy:

```bash
gzip -dk pbmc/pbmc.combined_dt04.ctdmr.tsv.gz
```

`manifest.tsv` records compressed and uncompressed sizes and the SHA-256
checksum of uncompressed contents for the 34 standard tissue catalog files.
The brain 5hmC, blood 5hmC, SMaHT production, and brain annotation catalogs have
their own `manifest.json` files in their respective directories.
The brain universal and PBMC directories
contain copies of the source data resolved from directory symlinks.
