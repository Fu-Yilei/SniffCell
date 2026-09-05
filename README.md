# SniffCell — 5hmC-compatible branch

[![PyPI version](https://img.shields.io/pypi/v/sniffcell.svg)](https://pypi.org/project/sniffcell/)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Issues](https://img.shields.io/badge/Issues-GitHub-red?logo=github)](https://github.com/Fu-Yilei/SniffCell/issues)

SniffCell links somatic structural variants (SVs) and tandem-repeat (TR)
expansions/contractions to the cell type they occurred in, using long-read DNA
methylation. It calls cell-type-specific differentially methylated regions
(ctDMRs) from a reference methylation atlas, discovers SVs and TRs directly
from long-read BAMs, extracts read-level methylation, and assigns each
variant's supporting reads to the cell population whose methylation profile
they match.

This branch supports separate 5mC and 5hmC atlas signals for cell-type
deconvolution.

## Which Package Should I Use?

Use **`sniffcell`** when you are discovering variants from a BAM (no existing
callset yet):

- ctDMR discovery from an atlas (`find`)
- cell-type deconvolution and BAM splitting (`deconv`)
- cell-type-aware SV, tandem-repeat, and SNV discovery (`discover`)
- methylation-based annotation of the discovered variants (`anno`, `svanno`)
- visualization, IGV screenshots, differential methylation tests, and HTML reports

Use **`sniffcell-lite`** when you already have a variant callset with each
variant's supporting reads identified (e.g. from Sniffles/kanpig or another
caller) and just want fast ctDMR-based cell-type annotation, without running
the full deconvolution/discovery pipeline:

- `sniffcell-lite find`
- `sniffcell-lite anno`

`sniffcell-lite` is published as a separate PyPI package and does not replace
the full `sniffcell` package:

```bash
pip install sniffcell-lite
sniffcell-lite --help
```

The lite package is maintained in the separate
[`Fu-Yilei/SniffCell-lite`](https://github.com/Fu-Yilei/SniffCell-lite)
repository.

## Install SniffCell

```bash
pip install "git+https://github.com/Fu-Yilei/SniffCell.git@5hmC_compatitible_atlas"
```

For the full external-tool environment:

```bash
micromamba env create -f environment.yml
micromamba activate sniffcell
pip install "git+https://github.com/Fu-Yilei/SniffCell.git@5hmC_compatitible_atlas"
```

## Methylation Atlases

Download the processed methylation atlas inputs and precomputed cell-type-specific
differentially methylated regions (ctDMRs) from
[Zenodo](https://zenodo.org/records/22003085). For the legacy NumPy `find` interface, place
`all_celltypes_blocks.npy`, `all_celltypes_blocks.index.gz`,
`index_to_major_celltypes.json`, and `all_celltypes.txt` in the `atlas/` directory.
Precomputed ctDMR tables are available for brain/cerebellum, PBMC, lung, liver,
pancreas, kidney, breast, and colon, and can be supplied directly to `deconv`
or `anno` with `-b` instead of running `find`.

To discuss tissue-specific atlases, additional tissues or cell types, or custom
atlas support, please [open a GitHub issue](https://github.com/Fu-Yilei/SniffCell/issues)
or [email us](mailto:yilei.fu@bcm.edu).

**Restricted 5hmC data:** The 5hmC signal used in the paper is subject to
access restrictions. For access inquiries and guidance on achieving better
performance with the 5hmC-compatible workflow, please
[email us](mailto:yilei.fu@bcm.edu).

## Example Workflows

For a runnable native-GRCh38 regional example with a bundled subset BAM and
atlas, see the [SH3RF3 dual SV/TR wiki tutorial](https://github.com/Fu-Yilei/SniffCell/wiki/SH3RF3-Dual-SV-TR-Example).

### SniffCell: discover and annotate SVs/TRs from a BAM

1. Call ctDMRs from an MDB atlas:

```bash
sniffcell find \
  --mdb combined_loyfer_ont.mmdb \
  --assay dual \
  -cf atlas/celltypes.json \
  -ck brain_cereb_ont \
  -o brain_dual_ctdmr.tsv \
  --diff_threshold 0.40
```

`--assay dual` calls separate 5mC and 5hmC views and records the assay in the
`modification` column; it does not collapse the two signals. The older
`--npy/--index/--meta` interface remains available as a legacy compatibility
path.

See the [custom MDB atlas Wiki tutorial](https://github.com/Fu-Yilei/SniffCell/wiki/Build-a-Custom-MDB-Atlas-from-bedMethyl)
for the complete bedMethyl-to-ctDMR workflow. Use an atlas with separate 5mC
and 5hmC assays and matching cell-type metadata; replace the example MDB and
metadata paths with your own inputs. See the access note above for the 5hmC
signal used in the paper.

2. Split the BAM into cell-type groups using those ctDMRs:

```bash
sniffcell deconv \
  -i sample.bam \
  -r ref.fa \
  -b brain_dual_ctdmr.tsv \
  -o deconv_out \
  --bam-modification auto \
  --split_bam_groups "Neuron=Neuron;Oligodendrocyte=Oligodendrocyte" \
  -t 8
```

With `--bam-modification auto` (the default), each ctDMR uses the channel named
in its `modification` column: BAM `C+m` calls for `5mC`, `C+h` calls for
`5hmC`, and both for legacy `modifiedC` rows. Use an explicit `5mC`, `5hmC`,
or `modifiedC` value only to override every catalog row for a control analysis.
Catalogs without a `modification` column retain the previous combined-modified-C
behavior.

3. Discover SVs and tandem repeats from the two cell-type-split BAMs:

```bash
sniffcell discover tools run \
  --deconv-dir deconv_out \
  --reference ref.fa \
  --tr-bed tr_catalog.bed \
  --sex female \
  --run-id run1 \
  --threads 8
```

This runs Sniffles/Kanpig/Truvari (SVs), TRGT or Medaka (tandem repeats), and
optionally Clair3 (SNVs) on the two split BAMs, then harmonizes everything
into `deconv_out/deconv_requested_group_splits/discover/run1/harmonized_variants.tsv`.
Tool binaries are resolved from `PATH` by default, or pass explicit
`--sniffles-bin`/`--trgt-bin`/etc. flags; run `sniffcell discover tools check`
to preflight-check dependencies first.

4. Annotate the discovered SVs/TRs with cell-type methylation evidence:

```bash
sniffcell anno \
  -i sample.bam \
  -v deconv_out/deconv_requested_group_splits/discover/run1/harmonized_variants.tsv \
  -r ref.fa \
  -b brain_dual_ctdmr.tsv \
  -o anno_out \
  -t 8
```

5. Build an HTML review report:

```bash
sniffcell report --anno_output anno_out
```

Open `anno_out/report/index.html` to review high-confidence calls and figures.

### SniffCell Lite: annotate an existing callset

If you already have a variant and know its supporting read names, skip
deconvolution and discovery entirely:

```bash
sniffcell-lite find -ck "Colon, Ascending" -o colon_ascending.ctdmr.tsv

sniffcell-lite anno \
  -i sample.bam \
  -r ref.fa \
  --variant-name variant_001 \
  --variant-location chr1:100000-101000 \
  --supporting-reads readA,readB,readC \
  --catalog colon_ascending.ctdmr.tsv \
  -o anno_out
```

`sniffcell-lite anno` also supports a `--batch variants.tsv` mode for
annotating many variants at once. See the
[SniffCell Lite README](https://github.com/Fu-Yilei/SniffCell-lite)
for details.

## Main Commands

| Command | Purpose |
|---------|---------|
| `sniffcell find` | Call ctDMRs from a methylation atlas |
| `sniffcell anno` | Extract BAM methylation and assign SVs/TRs to cell types |
| `sniffcell svanno` | Re-score SVs from a saved read-classification table |
| `sniffcell deconv` | Classify all reads and optionally split BAMs by group |
| `sniffcell discover` | Run cell-type-aware SV / TR / SNV discovery |
| `sniffcell viz` | Render per-SV methylation figures |
| `sniffcell igvviz` | Generate IGV batch screenshots |
| `sniffcell report` | Filter calls and build an HTML review report |
| `sniffcell dmsv` | Test methylation differences around SVs |

## Inputs

| Input | Format | Used by |
|-------|--------|---------|
| Long-read alignment | BAM with `MM` / `ML` modification tags | `anno`, `deconv`, `discover`, `dmsv`, `viz` |
| Variants | VCF / VCF.GZ, harmonized TSV from `discover`, or a variant + supporting-read names (`sniffcell-lite`) | `anno`, `dmsv`, `viz`, `report` |
| Reference genome | FASTA plus index | `anno`, `deconv`, `discover`, `dmsv`, `viz` |
| ctDMR table | TSV from `sniffcell find` | `anno`, `deconv`, `viz` |
| Methylation atlas | MDB atlas with separate assay views, or legacy NumPy matrix, CpG index, and metadata | `find` |

## Outputs

A `find -> deconv -> discover -> anno -> report` run produces:

```text
brain_dual_ctdmr.tsv
deconv_out/
  deconv_summary.tsv
  deconv_reads_classification.tsv
  deconv_requested_group_splits/
    <group>.bam
    discover/run1/
      harmonized_variants.tsv
      run_summary.json
anno_out/
  reads_classification.tsv
  blocks_classification.tsv
  sv_assignment.tsv
  sv_assignment_readable.tsv
  sv_assignment_readable_long.tsv
  anno_run_manifest.json
  report/
    index.html
    high_confidence_sv.tsv
    figures/
```

A `sniffcell-lite find -> anno` run produces:

```text
colon_ascending.ctdmr.tsv
colon_ascending.ctdmr.tsv.igv.bed
anno_out/
  variant_assignment.tsv
  variant_assignment_readable.tsv
  reads_classification.tsv
  anno_compact_manifest.json
```

## Documentation

For command-line options, run `sniffcell --help` or
`sniffcell <command> --help` (for example, `sniffcell anno --help`).

Full documentation is in the [GitHub Wiki](https://github.com/Fu-Yilei/SniffCell/wiki):

| Page | Contents |
|------|----------|
| [Installation](https://github.com/Fu-Yilei/SniffCell/wiki/Installation) | PyPI, conda, Docker, and tool setup |
| [SniffCell Lite](https://github.com/Fu-Yilei/SniffCell/wiki/SniffCell-Lite) | Annotating an existing callset with `sniffcell-lite` |
| [Find Workflow](https://github.com/Fu-Yilei/SniffCell/wiki/Find-Workflow) | ctDMR discovery parameters |
| [Methods](https://github.com/Fu-Yilei/SniffCell/wiki/Methods-Deconv-Discover-Anno) | Technical methods for core commands |
| [CLI Reference](https://github.com/Fu-Yilei/SniffCell/wiki/CLI-Reference) | Command-line options |
| [Test Examples](https://github.com/Fu-Yilei/SniffCell/wiki/Test-Examples) | Validation and QA examples |
| [SH3RF3 dual SV/TR example](https://github.com/Fu-Yilei/SniffCell/wiki/SH3RF3-Dual-SV-TR-Example) | Native-GRCh38 regional BAM example from atlas through report |

## Citation

If you use SniffCell in your research, please cite:

> **[Cell-type-resolved somatic variant discovery from bulk long-read sequencing](https://www.medrxiv.org/content/10.64898/2026.09.01.26361966v1)**
> Yilei Fu et al. *medRxiv preprint (2026).* DOI: 10.64898/2026.09.01.26361966.

## License

MIT License. See [LICENSE](LICENSE) for details.
