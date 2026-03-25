#!/usr/bin/env bash
set -euo pipefail

tool="/opt/conda/envs/medaka/bin/medaka"
if [[ ! -x "${tool}" ]]; then
  echo "medaka is not installed in this image. Rebuild with --build-arg INSTALL_MEDAKA=1." >&2
  exit 127
fi

exec "${tool}" "$@"
