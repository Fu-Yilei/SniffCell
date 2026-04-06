# SniffCell

[![PyPI version](https://img.shields.io/pypi/v/sniffcell.svg)](https://pypi.org/project/sniffcell/)
[![Python](https://img.shields.io/badge/python-%3E%3D3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/Docs-Wiki-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/wiki)
[![Issues](https://img.shields.io/badge/Issues-GitHub-red?logo=github)](https://github.com/Fu-Yilei/SniffCell/issues)

**SniffCell** is a Python toolkit for annotating somatic structural variants (SVs) with cell-type origin using long-read DNA methylation. It integrates cell-type-specific differentially methylated regions (ctDMRs) derived from a reference methylation atlas with per-read methylation measurements from nanopore or PacBio long-read BAMs to assign each SV — or every read in a sample — to a cell population.

---

## Why SniffCell?

Somatic SVs identified from bulk long-read sequencing are a mixture of events from different cell types. Without knowing the cell of origin, it is difficult to interpret their functional significance or estimate their true variant allele fraction within a specific compartment. SniffCell solves this by reading the epigenetic "fingerprint" imprinted on each DNA molecule and matching it against a reference atlas of cell-type-specific methylation patterns.

**Core capabilities:**

- **ctDMR discovery** — Mine a reference methylation atlas to find genomic regions with distinct methylation in each cell population
- **Read-level deconvolution** — Assign every read in a BAM to a cell type using ctDMR methylation signals, with no single-cell data required
- **SV annotation** — Link cell-type identity to SV-supporting reads and produce a per-SV cell-of-origin call
- **Discovery pipeline** — Run a full multi-stage SV / tandem-repeat / SNV calling workflow on cell-type-split BAMs produced by deconvolution
- **Interactive reporting** — Filter high-confidence SVs and generate an HTML review report with clickable per-SV figures and IGV screenshots

---

## Overview

![SniffCell workflow](https://raw.githubusercontent.com/Fu-Yilei/SniffCell/main/img/workflow.png)

The typical workflow has three main stages:

```
Atlas (NPY + index + metadata)
        │
        ▼
  sniffcell find         ← Call cell-type-specific DMRs (ctDMRs)
        │
        ▼
  sniffcell anno         ← Extract methylation from BAM, classify reads, assign SVs
        │
        ▼
  sniffcell report       ← Filter high-confidence calls, build HTML review report
        │
        ├── sniffcell viz        ← Per-SV methylation figure (PNG / PDF)
        ├── sniffcell igvviz     ← IGV batch screenshots
        └── sniffcell dmsv       ← Differential methylation test near each SV
```

For multi-group analyses (e.g., comparing SVs enriched in one cell compartment vs. another):

```
  sniffcell deconv       ← Deconvolve all reads; split BAM by cell type
        │
        ▼
  sniffcell discover     ← Call SVs / TRs / SNVs independently per group
        │
        ▼
  sniffcell anno         ← Annotate harmonized variants
```

---

## Quick Start

### 1. Install

```bash
pip install sniffcell
```

For the full environment including bioinformatics tools (Sniffles, bcftools, samtools, Truvari …):

```bash
micromamba env create -f environment.yml
micromamba activate sniffcell
pip install sniffcell
```

See [Installation](https://github.com/Fu-Yilei/SniffCell/wiki/Installation) in the wiki for Docker instructions, optional extras, and manual tool setup.

### 2. Call ctDMRs from the reference atlas

```bash
sniffcell find \
  -n atlas/all_celltypes_blocks.npy \
  -i atlas/all_celltypes_blocks.index.gz \
  -cf atlas/index_to_major_celltypes.json \
  -m atlas/all_celltypes.txt \
  -ck pbmc \
  -o pbmc_ctdmr.tsv
```

### 3. Annotate SVs with cell-type evidence

```bash
sniffcell anno \
  -i sample.bam \
  -v sample.vcf.gz \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -o anno_out \
  -t 8
```

### 4. Build the review report

```bash
sniffcell report --anno_output anno_out
```

Open `anno_out/report/index.html` in a browser to review filtered high-confidence SVs with per-SV methylation evidence.

---

## Commands at a Glance

| Command | What it does |
|---------|-------------|
| `sniffcell find` | Mine a reference atlas to call cell-type-specific DMRs (ctDMRs) |
| `sniffcell anno` | Extract read-level methylation from a BAM and assign each SV a cell-type code |
| `sniffcell svanno` | Re-score SV assignments from a saved read table without re-processing the BAM |
| `sniffcell deconv` | Assign every read in a BAM to a cell type; optionally split into per-group BAMs |
| `sniffcell discover` | Multi-stage SV / tandem-repeat / SNV pipeline on cell-type-split BAMs |
| `sniffcell viz` | Render a per-SV methylation figure (PNG or PDF) |
| `sniffcell igvviz` | Produce IGV batch-mode screenshots for a single SV |
| `sniffcell report` | Filter high-confidence SVs and build an interactive HTML review report |
| `sniffcell dmsv` | Test for differential methylation between SV-supporting and non-supporting reads |

---

## Input Requirements

| Input | Format | Used by |
|-------|--------|---------|
| Long-read alignment | BAM with `MM`/`ML` base-modification tags | `anno`, `deconv`, `dmsv`, `viz` |
| Structural variants | VCF / VCF.GZ or harmonized TSV from `discover` | `anno`, `dmsv`, `viz`, `report` |
| Reference genome | FASTA + index | `anno`, `deconv`, `dmsv`, `viz` |
| ctDMR table | TSV from `sniffcell find` | `anno`, `deconv`, `viz` |
| Methylation atlas | NumPy matrix + CpG index + metadata | `find` |

---

## Key Outputs

After a complete `find → anno → report` run, the outputs include:

```
pbmc_ctdmr.tsv                      ← Cell-type-specific DMRs (input to anno)
anno_out/
  reads_classification.tsv          ← Per-read × ctDMR methylation and cell-type assignment
  sv_assignment.tsv                 ← Per-SV cell-type code and quality metrics
  sv_assignment_readable.tsv        ← Human-readable version with expanded cell-type labels
  anno_run_manifest.json            ← Full run manifest (paths, parameters, versions)
  report/
    index.html                      ← Interactive HTML review report
    high_confidence_sv.tsv          ← Filtered high-confidence SVs
    figures/                        ← Per-SV methylation panels (when --with_figures)
```

---

## Deconvolution and Discovery

For samples where you want to compare SVs across cell populations:

```bash
# Step 1: Deconvolve reads and split into cell-type-specific BAMs
sniffcell deconv \
  -i sample.bam \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -o deconv_out \
  --split_bam_groups "lymph=t_cell,b_cell,nk_cell;myeloid=monocyte" \
  -t 8

# Step 2: Call SVs, tandem repeats, and SNVs on each group independently
sniffcell discover tools run \
  --deconv-dir deconv_out \
  --reference ref.fa \
  --tr-bed atlas/adotto.v2.trgt.bed \
  --sex female \
  --stages sv,tdb \
  --threads 16

# Step 3: Annotate the harmonized variants
sniffcell anno \
  -i sample.bam \
  -v deconv_out/discover/harmonized_variants.tsv \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -o anno_out
```

Before running `discover`, validate your environment:

```bash
sniffcell-check-discover --stages all
```

---

## Visualizing Individual SVs

```bash
# Minimal — loads all inputs automatically from the anno manifest
sniffcell viz --anno_output anno_out -s sniffles.SV123

# With table exports
sniffcell viz --anno_output anno_out -s sniffles.SV123 --export_tables

# IGV batch screenshot
sniffcell igvviz --anno_output anno_out -s sniffles.SV123
```

---

## Documentation

Full documentation lives in the [GitHub Wiki](https://github.com/Fu-Yilei/SniffCell/wiki):

| Page | Contents |
|------|---------|
| [Installation](https://github.com/Fu-Yilei/SniffCell/wiki/Installation) | PyPI, conda, Docker, manual tool setup, verification |
| [End-to-End Workflow](https://github.com/Fu-Yilei/SniffCell/wiki/End-to-End-Workflow) | Step-by-step walkthrough from atlas to HTML report |
| [Find Workflow](https://github.com/Fu-Yilei/SniffCell/wiki/Find-Workflow) | ctDMR discovery internals and parameter guide |
| [Methods](https://github.com/Fu-Yilei/SniffCell/wiki/Methods-Deconv-Discover-Anno) | Technical methods for `deconv`, `discover`, and `anno` |
| [Test Examples](https://github.com/Fu-Yilei/SniffCell/wiki/Test-Examples) | Practical validation and QA queries |

---

## Citation

If you use SniffCell in your research, please cite:

> **SniffCell: cell-type annotation of somatic structural variants using long-read methylation**
> Yilei Fu et al. *(manuscript in preparation)*

---

## License

MIT License — see [LICENSE](LICENSE) for details.

Developed at Baylor College of Medicine by [Yilei Fu](mailto:yilei.fu@bcm.edu).
