# SniffCell

[![PyPI version](https://img.shields.io/pypi/v/sniffcell.svg)](https://pypi.org/project/sniffcell/)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/Docs-Wiki-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/wiki)
[![Issues](https://img.shields.io/badge/Issues-GitHub-red?logo=github)](https://github.com/Fu-Yilei/SniffCell/issues)

SniffCell annotates somatic structural variants with cell-type evidence from
long-read DNA methylation. It uses cell-type-specific differentially methylated
regions (ctDMRs) from a reference methylation atlas, extracts read methylation
from BAM files, and links SV-supporting reads to likely cell populations.

## Which Package Should I Use?

Use **`sniffcell`** when you want the full toolkit:

- ctDMR discovery from an atlas
- BAM-backed SV annotation
- read-level deconvolution and optional BAM splitting
- cell-type-aware SV / TR / SNV discovery workflows
- visualization, IGV screenshots, differential methylation tests, and reports

Use **`sniffcell-lite`** when you only need the lightweight workflow:

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

## Minimal Workflow

Call ctDMRs from an atlas:

```bash
sniffcell find \
  -n atlas/all_celltypes_blocks.npy \
  -i atlas/all_celltypes_blocks.index.gz \
  -cf atlas/index_to_major_celltypes.json \
  -m atlas/all_celltypes.txt \
  -ck pbmc \
  -o pbmc_ctdmr.tsv
```

Annotate variants using methylation from a BAM:

```bash
sniffcell anno \
  -i sample.bam \
  -v sample.vcf.gz \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -o anno_out \
  -t 8
```

Build an HTML review report:

```bash
sniffcell report --anno_output anno_out
```

Open `anno_out/report/index.html` to review high-confidence calls and figures.

## Main Commands

| Command | Purpose |
|---------|---------|
| `sniffcell find` | Call ctDMRs from a methylation atlas |
| `sniffcell anno` | Extract BAM methylation and assign SVs to cell types |
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
| Long-read alignment | BAM with `MM` / `ML` modification tags | `anno`, `deconv`, `dmsv`, `viz` |
| Structural variants | VCF / VCF.GZ or harmonized TSV | `anno`, `dmsv`, `viz`, `report` |
| Reference genome | FASTA plus index | `anno`, `deconv`, `dmsv`, `viz` |
| ctDMR table | TSV from `sniffcell find` | `anno`, `deconv`, `viz` |
| Methylation atlas | NumPy matrix, CpG index, and metadata | `find` |

## Outputs

A typical `find -> anno -> report` run produces:

```text
pbmc_ctdmr.tsv
anno_out/
  reads_classification.tsv
  sv_assignment.tsv
  sv_assignment_readable.tsv
  anno_run_manifest.json
  report/
    index.html
    high_confidence_sv.tsv
    figures/
```

## Documentation

Full documentation is in the [GitHub Wiki](https://github.com/Fu-Yilei/SniffCell/wiki):

| Page | Contents |
|------|----------|
| [Installation](https://github.com/Fu-Yilei/SniffCell/wiki/Installation) | PyPI, conda, Docker, and tool setup |
| [End-to-End Workflow](https://github.com/Fu-Yilei/SniffCell/wiki/End-to-End-Workflow) | Atlas-to-report walkthrough |
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
