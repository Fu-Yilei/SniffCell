# Brain modifiedC + 5mC + 5hmC scoring atlas

[Download the compressed ctDMR catalog](brain_cereb.paper_backbone_plus_5mc_5hmc.four_donors.ctdmr.tsv.gz).
This GRCh38 atlas distinguishes Neuron and Oligodendrocyte and is intended
for the `5hmC_compatitible_atlas` branch of SniffCell.

| Component | Scoring rows |
|---|---:|
| Original paper modifiedC backbone | 152,961 |
| Added ONT 5mC | 194,662 |
| Added ONT 5hmC | 142,349 |
| Total | 489,972 |

The atlas contains 378,216 merged genomic components. Its four paired donors
are bcontrol1, bcontrol2, bmsa, and bparkinsons. The original paper backbone
fields are preserved exactly; new ONT markers overlap no backbone interval.
5mC and 5hmC rows can overlap each other and retain separate assay labels.
`manifest.json` records the source checksums, output checksum, and filter settings.
The ONT sources require four finite donor pairs, consistent direction,
median absolute effect >=0.40, at least three donor effects >=0.30, and all
four effects >=0.15. Discovery uses >=3 CpGs and a 3,000-bp merge gap.

```bash
sniffcell deconv \
  -i sample.bam -r GRCh38.fa \
  -b atlases/brain_5hmc/brain_cereb.paper_backbone_plus_5mc_5hmc.four_donors.ctdmr.tsv.gz \
  -o deconv_out --bam-modification auto \
  --read_assignment_mode closest_reference_mean \
  --per_read_min_agreement 0.66 -t 8
```

Use a BAM with separate m/h calls to supply both ONT channels. With an m-only
BAM, h-specific markers have no evidence. The m tag alone does not establish
whether the basecaller signal is separately measured 5mC or an older combined
modifiedC signal. Catalogs without modification labels still use the legacy
modifiedC path. The original precomputed tissue atlases from `main` remain
alongside this directory under [`atlases/`](../README.md).

The [held-out benchmark](https://github.com/Fu-Yilei/SniffCell-analysis/blob/5hmc/analyses/06_06_heldout_validation/BLOG_UPDATE.md)
uses a distinct three-donor atlas that excludes bcontrol1. Its results must
not be presented as an independent test of this all-four-donor artifact.
The catalog was structurally checked and the original backbone fields were
verified; it has not been independently benchmarked on new donors.

To reconstruct the integration with authorized source catalogs, run
`tools/build_dual_ctdmr_catalog.py --help`, then
`tools/export_5hmc_catalog.py --help` from the repository root. The export
retains only scoring fields and path-free provenance. Raw BAMs, MDBs,
individual read identifiers, and variant calls are not distributed here.
The catalog is a separate download, not part of the Python wheel.
