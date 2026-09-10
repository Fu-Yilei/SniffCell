# SniffCell Wiki

Welcome to the SniffCell documentation. SniffCell links somatic structural variants (SVs) and tandem-repeat (TR) expansions/contractions to cell-type origin by combining long-read DNA methylation signals with a reference methylation atlas.

---

## Getting Started

| Page | What you will find |
|------|--------------------|
| [Installation](Installation) | PyPI, conda, Docker, manual tool setup, and verification |
| [SniffCell Lite](SniffCell-Lite) | Lightweight `find` + `anno` workflow for an existing variant callset with known supporting reads |

---

## Reference

| Page | What you will find |
|------|--------------------|
| [Find Workflow](Find-Workflow) | How `sniffcell find` turns an atlas into ctDMRs — scoring, bipartition logic, and region merging |
| [Methods](Methods-Deconv-Discover-Anno) | Technical methods text for `deconv`, `discover`, and `anno` |
| [Test Examples](Test-Examples) | Practical validation runs for each command; QA queries and common failure patterns |
| [Deconv and Discover Design](Deconv-Postprocess-Design) | Internal design notes for the deconvolution and discovery pipeline |

---

## Command Summary

| Command | Purpose |
|---------|---------|
| `sniffcell find` | Call cell-type-specific DMRs (ctDMRs) from a reference atlas |
| `sniffcell anno` | Annotate SVs/TRs with read-level methylation and cell-type codes |
| `sniffcell svanno` | Re-score SV assignments from a saved read table |
| `sniffcell deconv` | Deconvolve all reads in a BAM by cell type; split into per-group BAMs |
| `sniffcell discover` | Multi-stage SV / tandem-repeat / SNV pipeline on cell-type-split BAMs |
| `sniffcell viz` | Render a per-SV methylation figure (PNG or PDF) |
| `sniffcell igvviz` | Produce IGV batch-mode screenshots for one SV |
| `sniffcell report` | Build an interactive HTML review report with filtered high-confidence SVs |
| `sniffcell dmsv` | Test for differential methylation near each SV |

---

## Quick links

- [GitHub repository](https://github.com/Fu-Yilei/SniffCell)
- [PyPI package](https://pypi.org/project/sniffcell/)
- [Bug reports and feature requests](https://github.com/Fu-Yilei/SniffCell/issues)
