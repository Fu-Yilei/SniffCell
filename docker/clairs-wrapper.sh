#!/usr/bin/env bash
set -euo pipefail

tool="/usr/local/bin/run_clairs.real"
if [[ -x "${tool}" ]]; then
  exec "${tool}" "$@"
fi

fallback="$(find /opt/sniffcell-tools/clairs -type f -name run_clairs | head -n1 || true)"
if [[ -n "${fallback}" && -x "${fallback}" ]]; then
  exec "${fallback}" "$@"
fi

echo "ClairS is not installed in this image. Rebuild with --build-arg CLAIRS_URL=<archive-or-binary-url>." >&2
exit 127
