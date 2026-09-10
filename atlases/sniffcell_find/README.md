# Precomputed ctDMR atlases

This directory contains the supplied `sniffcell_find` ctDMR tables and BED files
for brain/cerebellum, PBMC, breast, colon, kidney, liver, lung, and pancreas.
Both `combined_dt04` and `default035_gap2000_direct` file sets are preserved
as supplied, including their tissue-specific BED filenames.

All 32 source files are stored with lossless gzip compression. To extract an
individual file while retaining the compressed copy:

```bash
gzip -dk pbmc/pbmc.combined_dt04.ctdmr.tsv.gz
```

`manifest.tsv` records each file's compressed and uncompressed sizes and the
SHA-256 checksum of its uncompressed contents. The brain/cerebellum and PBMC
directories contain copies of the source data resolved from directory symlinks.
