# Publishing a 5hmC prerelease

The 5hmC branch can be installed directly today:

```bash
pip install "git+https://github.com/Fu-Yilei/SniffCell.git@5hmC_compatitible_atlas"
```

The `v0.9.8a1` alpha provides a versioned test release without replacing the
stable version selected by normal pip installs. See the
[release notes](releases/v0.9.8a1.md) and their benchmark blog link. Testers can
explicitly select the alpha:

```bash
pip install 'sniffcell==0.9.8a1'
```

Maintainer procedure:

1. Set `setup.cfg` to `version = 0.9.8a1` and `src/sniffcell/__init__.py` to
   `__version__ = "v0.9.8a1"` on the tested 5hmC branch.
2. Run the test suite, build distributions, and check package metadata.
3. Commit and push the version change, then create and push tag `v0.9.8a1`
   at that commit. Tag publication triggers the release workflow.

The workflow checks that the tag and both package versions agree. It publishes
alpha/beta/RC distributions to PyPI, marks the GitHub release as a prerelease,
and pushes a versioned Docker image without moving Docker `latest`. Stable
version tags retain the normal latest-release behavior. PyPI environment
protection and Trusted Publishing configuration still apply.

The atlas is a separate download under `atlases/5hmc/`; it is not bundled into
the Python wheel. Existing modifiedC atlases remain usable. The published
benchmark is a single-donor held-out test; the four-donor atlas is a distinct
artifact and has not been independently evaluated on those same four donors.

References: [Python versioning](https://packaging.python.org/en/latest/discussions/versioning/)
and [pip prerelease selection](https://pip.pypa.io/en/stable/cli/pip_install/).
