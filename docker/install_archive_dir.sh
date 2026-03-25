#!/usr/bin/env bash
set -euo pipefail

url="${1:-}"
install_dir="${2:?install dir required}"
subdir="${3:-}"

if [[ -z "${url}" ]]; then
  exit 0
fi

rm -rf "${install_dir}"
mkdir -p "${install_dir}"
tmpdir="$(mktemp -d)"
archive="${tmpdir}/payload"
extract_dir="${tmpdir}/extract"

cleanup() {
  rm -rf "${tmpdir}"
}
trap cleanup EXIT

mkdir -p "${extract_dir}"
curl -fsSL "${url}" -o "${archive}"

case "${url}" in
  *.tar.gz|*.tgz)
    tar -xzf "${archive}" -C "${extract_dir}"
    ;;
  *.tar.xz)
    tar -xJf "${archive}" -C "${extract_dir}"
    ;;
  *.zip)
    unzip -q "${archive}" -d "${extract_dir}"
    ;;
  *)
    echo "Unsupported model archive type for ${url}" >&2
    exit 1
    ;;
esac

source_dir="${extract_dir}"
if [[ -n "${subdir}" ]]; then
  source_dir="${extract_dir}/${subdir}"
fi

if [[ ! -d "${source_dir}" ]]; then
  first_dir="$(find "${extract_dir}" -mindepth 1 -maxdepth 1 -type d | head -n1 || true)"
  if [[ -n "${first_dir}" ]]; then
    source_dir="${first_dir}"
  fi
fi

if [[ ! -d "${source_dir}" ]]; then
  echo "Could not resolve extracted model directory from ${url}" >&2
  exit 1
fi

cp -a "${source_dir}/." "${install_dir}/"
