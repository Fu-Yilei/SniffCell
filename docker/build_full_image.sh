#!/usr/bin/env bash
set -euo pipefail

image_name="${1:-sniffcell:full}"

required_vars=(
  KANPIG_URL
  MODKIT_URL
  CLAIRS_URL
  CLAIR3_MODEL_URL
)

for name in "${required_vars[@]}"; do
  if [[ -z "${!name:-}" ]]; then
    echo "Missing required environment variable: ${name}" >&2
    exit 2
  fi
done

exec docker build \
  -t "${image_name}" \
  --build-arg FULL_DISCOVER=1 \
  --build-arg INSTALL_MEDAKA=1 \
  --build-arg INSTALL_CLAIR3=1 \
  --build-arg KANPIG_URL="${KANPIG_URL}" \
  --build-arg KANPIG_BIN_SUBPATH="${KANPIG_BIN_SUBPATH:-kanpig}" \
  --build-arg MODKIT_URL="${MODKIT_URL}" \
  --build-arg MODKIT_BIN_SUBPATH="${MODKIT_BIN_SUBPATH:-modkit}" \
  --build-arg CLAIRS_URL="${CLAIRS_URL}" \
  --build-arg CLAIRS_BIN_SUBPATH="${CLAIRS_BIN_SUBPATH:-run_clairs}" \
  --build-arg CLAIR3_MODEL_URL="${CLAIR3_MODEL_URL}" \
  --build-arg CLAIR3_MODEL_SUBDIR="${CLAIR3_MODEL_SUBDIR:-}" \
  .
