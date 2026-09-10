# Visualizing a 5mC/5hmC atlas

The development branch extends the existing `sniffcell viz` layout with
assay-aware ctDMR summaries. These changes are newer than the published
`v0.9.8a1` wheel; install the branch to use them:

```bash
pip install --upgrade "git+https://github.com/Fu-Yilei/SniffCell.git@5hmC_compatitible_atlas"
```

Read the [bmsa FGF14 blog and example figure](https://github.com/Fu-Yilei/SniffCell-analysis/blob/5hmc/analyses/06_07_assay_visualization/BLOG_UPDATE.md).

## Layout and measurements

The original read tracks, support arrows, haplotype filtering and CIGAR
insertion/deletion marks are retained. Read methylation is summarized over
ctDMR overlaps on supporting reads, rather than drawn as individual CpG marks.
The bottom panel retains genomic spacing and cell-type reference rows.

Exact-coordinate 5mC/5hmC marker pairs share a box. Red width is `m/(m+h)`;
blue width is `h/(m+h)`. Text gives absolute assay means on a 0–1 scale.
The interior split encodes the ratio, not a genomic boundary. Canonical C is
excluded from this ratio: independently summarized atlas means are not assumed
to form a measured three-state composition. If both means are zero, the ratio
is undefined and the box is neutral. A paired read box requires observations
of both channels at the same CpGs; missing evidence is not converted to zero.

Single-assay boxes retain a white-to-color mean scale. White means low signal
for that assay, not an independently established canonical-C fraction. The
paper modifiedC channel retains its historical extraction behavior. All read
ctDMR boxes have the same height; partial interval overlaps occupy separate
lanes. CIGAR insertion width is a display glyph, not a direct repeat-length
scale or an estimate of somatic excess over the sample baseline.

## Commands

```bash
sniffcell viz --anno_output sample_anno --sv_id VARIANT_ID \
  --sample_label "Sample name" -o sample.viz.png --export_tables

sniffcell report --anno_output sample_anno --with_figures -o sample_report
```

`report --with_figures` calls the same renderer, so its figures automatically
use this display when the atlas contains separate assay rows. Existing report
selection, review controls and cell-type filtering remain available.

`--linked_ctdmr_mode auto` is the default: separate-assay atlases show genomic
tracks, while legacy atlases keep the original distal-callout behavior.
`--flanking_ctdmrs 6` controls the number of nearest **distinct intervals**
per side of the variant in genomic mode; both assays are retained for a paired
interval. The requested base window is never reduced. For a fixed view use
`--linked_ctdmr_mode strict --window 8000 --exact_window`. In the assay-aware
view, `distal` uses genomic tracks too; `extend` retains its nearest-marker
window behavior. `--show_all_haplotypes` disables the default support-haplotype
filter when both haplotypes are wanted.

For saved regional deconvolution results, use `--read_assignment` with the
row-level classification table and `--read_summary` with the saved
`deconv_read_summary.tsv`. That filename is also detected beside the
classification table. Saved consensus labels take priority. Otherwise the
assay-aware view uses the classifier's intersection/conflict-threshold logic
on eligible variant-linked evidence, using the annotation manifest's threshold
when present (otherwise 0.66). This affects display labels, not scoring output.

## Explicit display exclusions

`--exclude_read_names excluded_read_ids.txt` omits listed IDs from the display
only. The figure prints the number omitted, and the supporting-read assignment
export records `display_excluded`. Original support membership, assignments,
and methylation exports are retained. Exclusions are not enabled by default.
Any selected illustration should state its selection in the caption and must
not be used to claim improved assignment accuracy.

The bmsa blog illustration excludes one borderline support read, whose existing
Oligodendrocyte assignment is retained. It shows 9/10 original support reads:
7 Neuron predictions and 2 ambiguous predictions. The unfiltered result is
7 Neuron, 1 Oligodendrocyte and 2 ambiguous. No read was relabeled to improve
this illustration, and no expansion threshold was changed.

## Compatibility

Legacy atlases with no modification column or only `modifiedC` continue through
the original extractor and renderer. A real-data regression check produced
exactly equal legacy extraction tables and a byte-identical PNG against the
pre-change branch. The scoring implementation and original atlas rows are
unchanged. Focused tests cover assay identity, matched CpGs, missing values,
proportional widths, equal heights, read consensus and legacy dispatch.
