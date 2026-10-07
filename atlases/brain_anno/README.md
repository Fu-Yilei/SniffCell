# Experimental brain annotation atlas

[Download the compressed ctDMR catalog](brain_universal.expanded9.modifiedC.ctdmr.tsv.gz).
This GRCh38 candidate contains **155,108 modifiedC regions** for `sniffcell anno`.
Its nine classes are `CorticalHippocampalNeuron`, `CerebellarNeuron`, `Astrocyte`,
`Oligodendrocyte`, `OPC`, `Microglia`, `Endothelial`, `Perivascular`, and `VLMC`.

Keep the original atlas for deconvolution and read splitting. To evaluate the
additional cell types during annotation, use the original unsplit BAM and this
catalog from the repository root:

```bash
sniffcell anno \
  -i sample.bam \
  -v harmonized_variants.tsv \
  -r GRCh38.fa \
  -b atlases/brain_anno/brain_universal.expanded9.modifiedC.ctdmr.tsv.gz \
  -o anno_brain_expanded9 \
  --read_assignment_mode closest_reference_mean \
  --evidence_mode per_read \
  -t 8
```

`anno` rescans the BAM with the supplied catalog. The former `--deconv-reads`
option has been removed. Annotation does not automatically combine the new
labels with original assignments or fall back to them when new evidence is
missing; retain the original results separately for comparison.

## Sources and construction

Original cortical, cerebellar and oligodendrocyte reference block means were
projected onto published [HBA major-cell-type DMRs](https://neomorph.salk.edu/hba/)
using CpG-overlap weighting. Partial blocks assume a uniform source block mean.
The cortical/hippocampal class uses equal weights for the original cortical
profile and the mean of hippocampal CA1, CA3, DG, HiF1 and HiF2 profiles.
Additional non-neuronal profiles come from the HBA matrix. These operations
do not reconstruct raw methylation counts.

The cerebellar source is [GSM5652231](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5652231),
titled Cerebellum-Neuron but described as "not sorted." Its purity is uncertain.
The [manifest](manifest.json) records the reference mappings, projection and
marker filters, build counts, and compressed/uncompressed SHA-256 checksums.

## Scope and validation

This candidate is restricted to selected published DMR intervals and does not
preserve all markers from the original atlas. It provides combined modifiedC
values, without separate 5mC/5hmC profiles for the added cell types.

Catalog structure and source-copy checks passed. The existing pilot used 901
molecule-by-locus observations from five frontal-cortex donors. Cell-of-origin
truth was unavailable, and unexpected cerebellar-reference-like assignments
require regional validation. Improved accuracy has not been established.
This catalog is a separate download, not part of the Python wheel.
