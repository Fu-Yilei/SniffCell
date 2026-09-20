# Precomputed ctDMR atlases

This directory contains the supplied `sniffcell_find` ctDMR tables and BED files
for brain/cerebellum, PBMC, breast, colon, kidney, liver, lung, pancreas, and
[sperm versus blood](sperm_vs_blood/README.md).
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

`manifest.tsv` records each file's compressed and uncompressed sizes and the
SHA-256 checksum of its uncompressed contents. The brain/cerebellum and PBMC
directories contain copies of the source data resolved from directory symlinks.
