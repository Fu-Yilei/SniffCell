# Bioconda recipes

Recipes for the published `sniffcell` and `sniffcell-lite` Python packages live
under `recipes/`. Each uses a versioned release source archive and SHA-256 checksum,
builds as `noarch: python`, and tests imports and command-line entry points.

The packages share the `sniffcell` Python namespace. The reciprocal
`run_constrained` entries prevent installing both in the same Conda environment.
Create separate environments for the full and lite packages.

These recipes install the core Python dependencies, matching the default PyPI
installs. Optional discovery callers, plotting/report extras, reference genomes,
and methylation atlases are not bundled. See the upstream README and
`environment.yml` for the full discovery environment.

## Build locally

With `conda-build` installed, run from the repository root:

```bash
conda build --override-channels -c conda-forge -c bioconda bioconda/recipes/sniffcell
conda build --override-channels -c conda-forge -c bioconda bioconda/recipes/sniffcell-lite
```

## Submit to Bioconda

Copy the two recipe directories into `recipes/` in a checkout of
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes).
Follow the [Bioconda contribution workflow](https://bioconda.github.io/contributor/workflow.html)
to lint, build, and open a pull request. These files alone do not publish packages
to the Bioconda channel.

## Update a release

Update the recipe version and the SHA-256 checksum from the corresponding
source distribution, then reset the build number to zero. Check dependencies,
entry points, and the license in that release's archive; they may differ from
the current upstream branch.

Both SniffCell 0.9.7 and SniffCell-lite 0.9.7 use verified PyPI source
distributions.

Validation: source checksums, rendered YAML, source installation, Python imports,
and all listed CLI help tests were checked locally. A full Conda build and
Bioconda lint run remain to be performed with the Bioconda build toolchain.
