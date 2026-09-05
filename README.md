# SniffCell

[![PyPI version](https://img.shields.io/pypi/v/sniffcell.svg)](https://pypi.org/project/sniffcell/)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/Docs-Wiki-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/wiki)
[![Issues](https://img.shields.io/badge/Issues-GitHub-red?logo=github)](https://github.com/Fu-Yilei/SniffCell/issues)

SniffCell links somatic structural variants (SVs) and tandem-repeat (TR)
expansions/contractions to the cell type they occurred in, using long-read DNA
methylation. It calls cell-type-specific differentially methylated regions
(ctDMRs) from a reference methylation atlas, discovers SVs and TRs directly
from long-read BAMs, extracts read-level methylation, and assigns each
variant's supporting reads to the cell population whose methylation profile
they match.

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

The lite package is maintained on the
[`sniffcell-lite`](https://github.com/Fu-Yilei/SniffCell/tree/sniffcell-lite)
branch.

## Install SniffCell

```bash
pip install sniffcell
```

For the full external-tool environment:

```bash
micromamba env create -f environment.yml
micromamba activate sniffcell
pip install sniffcell
```

See the [Installation wiki](https://github.com/Fu-Yilei/SniffCell/wiki/Installation)
for Docker, optional tools, and manual setup.

## Example Workflows

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
for the complete bedMethyl-to-ctDMR workflow. The brain modifiedC/5mC/5hmC
super-union catalog is distributed as
`atlas/brain_cereb.dual_5mc_5hmc.ctdmr.tsv.gz`.

2. Split the BAM into cell-type groups using those ctDMRs:

```bash
sniffcell deconv \
  -i sample.bam \
  -r ref.fa \
  -b atlas/brain_cereb.dual_5mc_5hmc.ctdmr.tsv.gz \
  -o deconv_out \
  --bam-modification auto \
  --split_bam_groups "Lymphoid=T-cell,NK-cell,B-cell;Myeloid=Monocyte" \
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
  -b pbmc_ctdmr.tsv \
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
[sniffcell-lite branch README](https://github.com/Fu-Yilei/SniffCell/tree/sniffcell-lite)
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
| Methylation atlas | MDB atlas (preferred), or legacy NumPy matrix plus CpG index and metadata | `find` |

## Outputs

A `find -> deconv -> discover -> anno -> report` run produces:

```text
pbmc_ctdmr.tsv
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

Full documentation is in the [GitHub Wiki](https://github.com/Fu-Yilei/SniffCell/wiki):

| Page | Contents |
|------|----------|
| [Installation](https://github.com/Fu-Yilei/SniffCell/wiki/Installation) | PyPI, conda, Docker, and tool setup |
| [SniffCell Lite](https://github.com/Fu-Yilei/SniffCell/wiki/SniffCell-Lite) | Annotating an existing callset with `sniffcell-lite` |
| [Find Workflow](https://github.com/Fu-Yilei/SniffCell/wiki/Find-Workflow) | ctDMR discovery parameters |
| [Methods](https://github.com/Fu-Yilei/SniffCell/wiki/Methods-Deconv-Discover-Anno) | Technical methods for core commands |
| [CLI Reference](https://github.com/Fu-Yilei/SniffCell/wiki/CLI-Reference) | Command-line options |
| [Test Examples](https://github.com/Fu-Yilei/SniffCell/wiki/Test-Examples) | Validation and QA examples |

## Citation

If you use SniffCell in your research, please cite:

> **SniffCell: cell-type annotation of somatic structural variants using long-read methylation**
> Yilei Fu et al. *(manuscript in preparation)*

## License

MIT License. See [LICENSE](LICENSE) for details.

Developed at Baylor College of Medicine by [Yilei Fu](mailto:yilei.fu@bcm.edu).
