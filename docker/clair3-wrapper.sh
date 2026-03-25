#!/usr/bin/env bash
set -euo pipefail

tool="/opt/conda/envs/clair3/bin/run_clair3.sh"
if [[ ! -x "${tool}" ]]; then
  echo "Clair3 is not installed in this image. Rebuild with --build-arg INSTALL_CLAIR3=1." >&2
  exit 127
fi

exec "${tool}" "$@"
