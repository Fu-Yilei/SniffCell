# Brain cerebellum dual-modification ctDMR catalog

Catalog file: `brain_cereb.dual_5mc_5hmc.ctdmr.tsv.gz`

- GRCh38 coordinates.
- 478,514 channel-aware rows.
- All 138,094 historical `brain_cereb` modifiedC ctDMRs are preserved with
  identical coordinates.
- 197,393 5mC and 143,027 5hmC ctDMRs novel to the historical backbone are
  appended.
- The row identity is `(chr, start, end, modification)`.
- 5mC and 5hmC values must not be summed or averaged during downstream use.

ONT selection criteria:

- absolute median paired neuron-minus-oligodendrocyte effect at least 0.40;
- the same direction in all four donor pairs;
- at least three donor pairs with an absolute effect at least 0.30;
- no donor pair with an absolute effect below 0.15.

SHA-256:

```text
d947004d2f4a0fa06dfaa6316c5a02f6a03560c2fc6b72648ef65f53859a7cbf  brain_cereb.dual_5mc_5hmc.ctdmr.tsv.gz
```

See the [Wiki tutorial](https://github.com/Fu-Yilei/SniffCell/wiki/Build-a-Custom-MDB-Atlas-from-bedMethyl)
for instructions to generate a custom atlas from bedMethyl files.
