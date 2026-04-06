# Installation

This page covers all installation paths for SniffCell — from a quick `pip install` for core annotation workflows to a full conda environment or Docker image for the multi-stage `discover` pipeline.

---

## Contents

- [Requirements](#requirements)
- [Option A — PyPI (core commands only)](#option-a--pypi-core-commands-only)
- [Option B — Conda environment (recommended for discover)](#option-b--conda-environment-recommended-for-discover)
- [Option C — Docker](#option-c--docker)
- [Optional Python extras](#optional-python-extras)
- [Manual tool installation](#manual-tool-installation)
- [Verifying your installation](#verifying-your-installation)
- [Troubleshooting](#troubleshooting)

---

## Requirements

| Requirement | Minimum version | Notes |
|-------------|----------------|-------|
| Python | 3.10 | Required |
| pysam | 0.21.0 | Installed automatically |
| numpy | 2.2.0 | Installed automatically |
| pandas | 2.3.0 | Installed automatically |
| scipy | any | Installed automatically |
| scikit-learn | any | Installed automatically |
| matplotlib | any | Installed automatically |
| tqdm | any | Installed automatically |

For the `discover` pipeline, additional bioinformatics tools are required. See [Manual tool installation](#manual-tool-installation) below.

---

## Option A — PyPI (core commands only)

This is the fastest path if you only need `find`, `anno`, `svanno`, `deconv`, `viz`, `igvviz`, `report`, and `dmsv`.

```bash
pip install sniffcell
```

Verify the install:

```bash
sniffcell --help
```

For optional features (`discover` post-processing or `igv-reports` integration):

```bash
pip install "sniffcell[discover]"        # adds tdb + seaborn for TR summary plots
pip install "sniffcell[igvreport]"       # adds igv-reports for HTML IGV.js pages
pip install "sniffcell[full]"            # all of the above
```

---

## Option B — Conda environment (recommended for discover)

The `environment.yml` at the root of the repository creates a conda environment with SniffCell and all tools available through conda/bioconda: `sniffles`, `bcftools`, `samtools`, `tabix`, `bgzip`, `truvari`, and the Python `tdb` package.

### Step 1 — Create the environment

Using micromamba (recommended):

```bash
micromamba env create -f environment.yml
micromamba activate sniffcell
```

Or with conda/mamba:

```bash
conda env create -f environment.yml
conda activate sniffcell
```

### Step 2 — Install the Python package

From PyPI:

```bash
pip install sniffcell
```

Or from a local checkout (editable install, useful for development):

```bash
git clone https://github.com/Fu-Yilei/SniffCell.git
cd SniffCell
pip install -e ".[full]"
```

### Step 3 — Install tools not available from conda

Some `discover` tools are not available through conda and require manual installation. See [Manual tool installation](#manual-tool-installation) for `kanpig`, `modkit`, `medaka`, and `clair3`.

### Step 4 — Run the environment preflight

Before running `sniffcell discover`, validate that all required binaries are available:

```bash
# Check all stages
sniffcell-check-discover --stages all

# Check only the SV and methylation stages
sniffcell-check-discover --stages sv,mods

# Check clair3 with a specific model path
sniffcell-check-discover --stages clair3 --clair3-model-path /path/to/clair3_model

# Check medaka and tdb with custom binary paths
sniffcell-check-discover \
  --stages medaka,tdb \
  --medaka-bin /path/to/medaka \
  --tdb-bin /path/to/tdb
```

The checker validates stage-specific binaries and will warn you of any missing dependencies before you invest compute time in a run.

---

## Option C — Docker

The repository ships with a unified `Dockerfile` that covers both the core `sniffcell` commands and the `discover` pipeline.

### Base image (flexible)

Includes: Python package + Sniffles + bcftools + samtools + tabix + bgzip + Truvari + Python `tdb` + dedicated conda envs for `medaka` and `clair3`.

```bash
docker build -t sniffcell:latest .
```

Tools not available from conda (`kanpig`, `modkit`) are included as optional build arguments:

```bash
docker build -t sniffcell:latest \
  --build-arg KANPIG_URL=https://example.com/kanpig.tar.gz \
  --build-arg MODKIT_URL=https://example.com/modkit \
  .
```

### Full end-to-end image (strict)

For a container that is self-contained for `sniffcell discover --stages all`, use the tracked builder helper. This build requires `kanpig`, `modkit`, and a Clair3 model to be supplied, and will fail if any are missing.

```bash
KANPIG_URL=https://example.com/kanpig.tar.gz \
MODKIT_URL=https://example.com/modkit \
CLAIR3_MODEL_URL=https://example.com/clair3_model.tar.gz \
docker/build_full_image.sh sniffcell:full
```

Before building on a new machine, check that your Docker runtime is usable:

```bash
docker/check_builder_host.sh
```

### Running the container

Preflight check:

```bash
docker run --rm sniffcell:latest sniffcell-check-discover --stages all
```

Run a discover pipeline on local data:

```bash
docker run --rm -it \
  -v /path/to/your/data:/data \
  sniffcell:full \
  discover tools run \
  --deconv-dir /data/sample/deconv \
  --reference /data/ref.fa \
  --tr-bed /data/tr.bed \
  --sex female \
  --clair3-model-path /opt/models/clair3
```

> The container entrypoint accepts both full commands (`sniffcell anno ...`) and short forms (`anno ...`).

---

## Optional Python extras

| Extra | Installs | When you need it |
|-------|----------|-----------------|
| `discover` | `tdb`, `seaborn` | TR summary plots in `sniffcell discover` |
| `igvreport` | `igv-reports` | Alternate IGV.js HTML page in `sniffcell report --with_igvreport` |
| `full` | All of the above | Full feature set |

```bash
pip install "sniffcell[discover]"
pip install "sniffcell[igvreport]"
pip install "sniffcell[full]"
```

---

## Manual tool installation

These tools are used by `sniffcell discover` but are not available from conda. Follow the instructions on each project's GitHub releases page.

### kanpig (SV re-genotyping)

Download the pre-compiled binary from the [kanpig releases page](https://github.com/ACEnglish/kanpig/releases) and place it somewhere on your `PATH`:

```bash
# Example: download and install to ~/bin
curl -L https://github.com/ACEnglish/kanpig/releases/latest/download/kanpig-linux-x86_64.tar.gz \
  | tar -xz -C ~/bin
chmod +x ~/bin/kanpig
```

### modkit (ONT methylation pileups)

Download the pre-compiled binary from the [modkit releases page](https://github.com/nanoporetech/modkit/releases):

```bash
curl -L https://github.com/nanoporetech/modkit/releases/latest/download/modkit -o ~/bin/modkit
chmod +x ~/bin/modkit
```

### medaka (tandem repeat calling)

medaka should be installed in its own dedicated conda environment to avoid dependency conflicts:

```bash
micromamba create -n medaka_env -c conda-forge -c bioconda medaka
```

Point `sniffcell discover` to the medaka binary:

```bash
sniffcell discover tools run \
  ... \
  --medaka-bin /path/to/medaka_env/bin/medaka
```

### clair3 (germline SNV/indel calling)

Install clair3 in a dedicated conda environment and download the appropriate model for your sequencing platform:

```bash
micromamba create -n clair3_env -c conda-forge -c bioconda clair3
```

Model download example (ONT R10 chemistry):

```bash
# Download from the Clair3 model repository
wget -qO- https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v430.tar.gz \
  | tar -xz -C /opt/clair3_models/
```

Pass the model path at runtime:

```bash
sniffcell discover tools run \
  ... \
  --clair3-model-path /opt/clair3_models/r1041_e82_400bps_sup_v430
```

---

## Verifying your installation

### Smoke test (fresh install)

The repository ships with helper scripts that exercise the package on synthetic data:

```bash
# Test a wheel install
scripts/check_fresh_install.sh wheel

# Test an editable install
scripts/check_fresh_install.sh editable
```

### Unit tests (development)

```bash
python -m unittest discover -s tests -v
```

### Quick CLI check

```bash
sniffcell --help
sniffcell find --help
sniffcell anno --help
sniffcell discover tools run --help
```

---

## Troubleshooting

**`sniffcell: command not found` after `pip install`**

Your Python `bin` directory may not be on `PATH`. Add it:

```bash
export PATH="$(python -m site --user-base)/bin:$PATH"
```

**`pysam` fails to import on macOS**

Install `htslib` via Homebrew first, then reinstall pysam:

```bash
brew install htslib
pip install --force-reinstall pysam
```

**`sniffcell discover` fails preflight with missing binaries**

Run `sniffcell-check-discover --stages <stage>` to identify which tools are missing and follow the [Manual tool installation](#manual-tool-installation) guide above.

**Import errors after a conda environment update**

Reinstall the package after updating the environment:

```bash
pip install --force-reinstall sniffcell
```

**Low memory when running `sniffcell find` on the full atlas**

The full atlas NPY matrix (~8 GB) is loaded entirely into RAM. On machines with less than 32 GB of memory, close other applications or use a subset of cell types via the `-ck` option.
