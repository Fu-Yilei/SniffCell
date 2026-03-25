# SniffCell
[![PyPI version](https://img.shields.io/pypi/v/sniffcell.svg)](https://pypi.org/project/sniffcell/)
[![Install](https://img.shields.io/badge/Install-PyPI-3776AB?logo=pypi&logoColor=white)](https://pypi.org/project/sniffcell/)
[![Docs](https://img.shields.io/badge/Docs-GitHub-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/wiki)
[![Issues](https://img.shields.io/badge/Issues-GitHub-181717?logo=github)](https://github.com/Fu-Yilei/SniffCell/issues)

SniffCell annotates structural variants (SVs) using long-read methylation evidence and cell-type-specific ctDMR signals.

## Installation

```bash
pip install sniffcell
```

Requires Python `>=3.10`.

For a local checkout:

```bash
python -m pip install -e .
```

Optional Python extras:

```bash
python -m pip install ".[discover,igvreport]"
```

- `discover` installs the Python-side `tdb` package plus `seaborn` for TR summary plots in `sniffcell discover`.
- `igvreport` installs `igv-reports` for the alternate IGV.js HTML page in `sniffcell report`.

For a full `discover` environment with the conda-available tools preinstalled:

```bash
micromamba env create -f environment.yml
micromamba activate sniffcell
python -m pip install -e .
```

`environment.yml` covers the Python package plus conda/bioconda tools such as `sniffles`, `bcftools`, `samtools`, `truvari`, `severus`, and the Python `tdb` package. Some `discover` dependencies still require separate installation or dedicated envs: `kanpig`, `modkit`, `medaka`, `clair3`, and `ClairS`.

Before launching `sniffcell discover`, run a preflight:

```bash
sniffcell-check-discover --stages all
```

Examples:

```bash
sniffcell-check-discover --stages sv,mods
sniffcell-check-discover --stages clair3 --clair3-model-path /path/to/clair3_model
sniffcell-check-discover --stages medaka,tdb --medaka-bin /path/to/medaka --tdb-bin /path/to/tdb
```

The checker validates stage-specific binaries such as `sniffles`, `bcftools`, `bgzip`, `tabix`, `kanpig`, `truvari`, `medaka`, `tdb`, `modkit`, `run_clair3.sh`, and `run_clairs`. It also verifies the Python `tdb` package for TR postprocessing, warns when `seaborn` is missing for optional TR plots, and enforces `--clair3-model-path` when the `clair3` stage is requested.

Fresh-install smoke test helpers:

```bash
scripts/check_fresh_install.sh wheel
scripts/check_fresh_install.sh editable
```

## Docker

The repo now includes one unified `Dockerfile` for both the core `sniffcell` commands and the `discover` pipeline.

Build the base image:

```bash
docker build -t sniffcell:latest .
```

That image includes the Python package, `sniffles`, `bcftools`, `bgzip`, `tabix`, `samtools`, `truvari`, `severus`, the Python `tdb` package, and by default separate conda envs for `medaka` and `clair3`.

Optional build args for tools that are not reliably available from conda:

```bash
docker build -t sniffcell:latest \
  --build-arg KANPIG_URL=https://.../kanpig.tar.gz \
  --build-arg MODKIT_URL=https://.../modkit \
  --build-arg CLAIRS_URL=https://.../clairs.tar.gz \
  .
```

Optional archive layouts can be adjusted with:

```bash
--build-arg KANPIG_BIN_SUBPATH=kanpig
--build-arg MODKIT_BIN_SUBPATH=modkit
--build-arg CLAIRS_BIN_SUBPATH=run_clairs
```

For a true end-to-end `discover` image with no additional software installs at runtime, use the tracked builder helper:

```bash
KANPIG_URL=https://... \
MODKIT_URL=https://... \
CLAIRS_URL=https://... \
CLAIR3_MODEL_URL=https://... \
docker/build_full_image.sh sniffcell:full
```

That build runs in strict mode:
- it fails unless `kanpig`, `modkit`, `ClairS`, and a Clair3 model archive are supplied
- it keeps `medaka` and `clair3` installed in dedicated conda envs
- it runs `sniffcell-check-discover --stages all --clair3-model-path /opt/models/clair3` during the image build

If the full build succeeds, the container is self-contained for `sniffcell discover --stages all`.

Before trying the full build on a new machine, you can check whether the local OCI runtime is usable:

```bash
docker/check_builder_host.sh
```

Example container preflight:

```bash
docker run --rm sniffcell:latest sniffcell-check-discover --stages all
docker run --rm sniffcell:latest sniffcell-check-discover --stages clair3 --clair3-model-path /models/clair3
```

Example runtime:

```bash
docker run --rm -it \
  -v /path/to/data:/data \
  sniffcell:full \
  discover \
  --deconv-dir /data/sample/deconv \
  --reference /data/ref.fa \
  --tr-bed /data/tr.bed \
  --sex female \
  --clair3-model-path /opt/models/clair3
```

Notes:
- `sniffcell:latest` is still the flexible base image.
- `sniffcell:full` is the strict end-to-end discover image produced by `docker/build_full_image.sh`.
- The full image only becomes self-contained if the external build URLs stay valid.
- The container entrypoint accepts either full commands like `sniffcell anno ...` or short forms like `anno ...`.

## Commands

```
sniffcell {find, deconv, anno, svanno, dmsv, viz, igvviz, report}
```

## Typical Workflow

1. Call ctDMRs from an atlas with `find`.
3. Annotate SVs with ctDMR evidence using `anno`.
4. Re-run SV assignment from saved read tables with `svanno` (optional).
5. Generate an HTML review report with `report`.
6. Visualize individual SVs with `viz` or `igvviz`.
7. Test differential methylation near SVs with `dmsv` (optional).
8. Deconvolve cell-type composition from any BAM with `deconv` (optional).

---

## `find`: Call ctDMRs From an Atlas

Loads an atlas methylation matrix and calls cell-type-specific differentially methylated regions (ctDMRs).

```bash
sniffcell find \
  -n atlas/all_celltypes_blocks.npy \
  -i atlas/all_celltypes_blocks.index.gz \
  -cf atlas/index_to_major_celltypes.json \
  -m atlas/all_celltypes.txt \
  -ck pbmc \
  -o pbmc_ctdmr.tsv \
  --diff_threshold 0.40 \
  --min_rows 2 \
  --min_cpgs 3 \
  --max_gap_bp 500
```

> If `-n/-i/-cf/-m` are omitted, paths default to `./atlas/...` in your working directory.

**`-ck/--celltypes_keys`** selects a top-level JSON key mapping `{group_name: [sample_id, ...]}`.

**Outputs:**
- `<output>` — annotation-ready ctDMR TSV
- `<output>.igv.bed` — IGV BED9 companion file

---

## `anno`: Annotate SVs With ctDMRs

Classifies reads per ctDMR region, then assigns cell-type codes to each SV.

```bash
sniffcell anno \
  -i sample.bam \
  -v sample.vcf.gz \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -o anno_out \
  -w 10000 \
  --breakpoint_exclusion_frac 0.1 \
  -t 8 \
  --evidence_mode all_rows \
  --min_overlap_pct 0.0 \
  --min_agreement_pct 1.0
```

**Key options:**
- `--evidence_mode {all_rows,per_read}` — how ctDMR evidence is aggregated (default: `all_rows`)
- `--breakpoint_exclusion_frac` — excludes ctDMRs within `±frac × SVLEN` of the breakpoint (default: `0.0`)
- `--min_overlap_pct` / `--min_agreement_pct` — filtering thresholds

> `assigned_code` is suppressed when `has_hard_conflict=True`.

**Outputs in `<output>/`:**
- `reads_classification.tsv`
- `blocks_classification.tsv`
- `sv_assignment.tsv` / `sv_assignment_readable.tsv` / `sv_assignment_readable_long.tsv`
- `anno_run_manifest.json`

---

## `svanno`: Recompute SV Assignments

Re-runs only the SV assignment step from an existing `reads_classification.tsv`, useful for tuning thresholds without re-processing the BAM.

```bash
sniffcell svanno \
  -v sample.vcf.gz \
  -i anno_out/reads_classification.tsv \
  -w 10000 \
  --breakpoint_exclusion_frac 0.1 \
  --evidence_mode all_rows \
  --min_overlap_pct 0.0 \
  --min_agreement_pct 1.0 \
  -o anno_out
```

---

## `deconv`: Cell-Type Deconvolution

Assigns every read in a BAM a cell-type code using ctDMR methylation patterns, then produces per-read, per-group, and whole-sample summaries.

```bash
sniffcell deconv \
  -i sample.bam \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -o deconv_out \
  -t 8 \
  --read_assignment_mode closest_reference_mean
```

**Key options:**
- `--read_assignment_mode {closest_reference_mean,kmeans}` — assignment algorithm (default: `closest_reference_mean`)
- `--split_bam_groups` — after deconvolution, split reads into per-group BAMs. Use `;` between groups and `,` between labels within a group. Named splits use `=`. Example: `lymph=t_cell,b_cell,nk_cell;myeloid=monocyte`
- `--resume` — skip ctDMR classification and reload existing TSVs; useful for re-splitting without reprocessing
- `--skip_overall_summary` — skip writing `deconv_summary.tsv`; useful when you only need per-read outputs and split BAMs

**Outputs in `<output>/`:**
- `deconv_reads_classification.tsv` — one row per (read × ctDMR)
- `deconv_blocks_classification.tsv` — per-ctDMR block summary
- `deconv_read_summary.tsv` — one row per read with majority cell type and linked celltypes
- `deconv_summary.tsv` — whole-sample summary in `all_rows` and `per_read` modes
- `deconv_reads_by_group/` — per-group read tables (split by `best_group`)
- `deconv_requested_group_splits/` — user-defined BAM and TSV splits (when `--split_bam_groups` is used)
- `deconv_run_manifest.json`

---

## `viz`: Visualize One SV

Renders a figure (PNG or PDF) for a single SV with read-level methylation and ctDMR context.

```bash
# Minimal — loads inputs from anno manifest
sniffcell viz \
  --anno_output anno_out \
  -s sniffles.SV123

# Manual mode
sniffcell viz \
  -i sample.bam \
  -v sample.vcf.gz \
  -s sniffles.SV123 \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -a anno_out/reads_classification.tsv \
  -o figures/sniffles.SV123 \
  -f png
```

**Notable options:**
- `--indel_min_bp` — overlay read-level indels ≥ N bp on reads (default: `40`; set to `0` to disable)
- `--linked_ctdmr_mode {distal,extend,strict}` — controls how off-window winning ctDMRs are displayed (default: `distal`)
- `--export_tables` — also write `.summary.tsv`, `.supporting_reads_assignment.tsv`, and `.supporting_reads_ctdmr_methylation.tsv`

---

## `igvviz`: IGV Screenshots for One SV

Runs IGV batch mode and produces snapshots per BAM, with reads tagged and grouped by phase.

```bash
sniffcell igvviz \
  -i fans_a.bam fans_b.bam fans_c.bam \
  -v sample.vcf.gz \
  -s sniffles.SV123 \
  -r ref.fa \
  -b pbmc_ctdmr.tsv \
  -w 10000 \
  -o out/igvviz
```

**Notable options:**
- `--anno_output` — load inputs from anno manifest (manifest-driven mode)
- `--igv_cmd` — path to IGV executable (default: `igv.sh`)
- `--snapshot_width/--snapshot_height` — snapshot dimensions (default: `3600×1600`)
- `--batch_only` — write batch script only, don't run IGV

---

## `report`: HTML Review Report

Filters high-confidence SVs from `anno` output and builds an interactive HTML report.

```bash
# Basic report
sniffcell report \
  --anno_output anno_out \
  --min_overlap_pct 0.8 \
  --overlap_filter_mode gradient \
  --min_majority_pct 1.0

# With viz figures and IGV screenshots
sniffcell report \
  --anno_output anno_out \
  --with_figures \
  --with_igvviz \
  --igv_bams fans1.bam fans2.bam fans3.bam \
  --figure_threads 4

# With igv-reports alternate page (requires: pip install igv-reports)
sniffcell report \
  --anno_output anno_out \
  --with_igvreport \
  --igv_bams fans1.bam fans2.bam
```

**Default SV filters:**
- `assigned_code` must be non-empty
- `linked_celltypes` must be non-empty
- `has_hard_conflict` must be `False`
- `--overlap_filter_mode gradient` with `--min_overlap_pct 0.8`
- `--min_majority_pct` ≥ `1.0`

`gradient` uses `ceil(min_overlap_pct * n_supporting^exponent)` overlapped reads, with default `--overlap_gradient_exponent 0.5`. Use `--overlap_filter_mode hard_clip` to recover the previous fixed `overlap_pct` behavior.

**Outputs under `<anno_output>/report/`:**
- `index.html` — interactive report with genome-wide plots and per-SV panels
- `high_confidence_sv.tsv`
- `figures/` — viz panels (when `--with_figures`)
- `igvviz/` — IGV screenshots (when `--with_igvviz`)
- `igvreport/index.html` — alternate IGV.js page (when `--with_igvreport`)
- `report_manifest.json`

> Review labels (Real / Not real / Undecided) auto-save to browser `localStorage` and persist across sessions.

---

## `dmsv`: Differential Methylation Around SVs

Tests for methylation differences between SV-supporting and non-supporting reads near each SV.

```bash
sniffcell dmsv \
  -i sample.bam \
  -v sample.vcf.gz \
  -r ref.fa \
  -o dmsv_out \
  -m 3 \
  -f 1000 \
  -c 5 \
  -t 8
```

**Outputs:**
- `dmsv_out/significant_SVs.tsv`
- `dmsv_out/sv_details/<sv_id>.tsv.gz`

---

## Wiki
- End-to-end workflow: [`wiki/End-to-End-Workflow.md`](wiki/End-to-End-Workflow.md)
- Test examples: [`wiki/Test-Examples.md`](wiki/Test-Examples.md)
