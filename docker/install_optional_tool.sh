#!/usr/bin/env bash
set -euo pipefail

url="${1:-}"
install_dir="${2:?install dir required}"
binary_name="${3:?binary name required}"
subpath="${4:-$binary_name}"

if [[ -z "${url}" ]]; then
  exit 0
fi

mkdir -p "${install_dir}"
tmpdir="$(mktemp -d)"
archive="${tmpdir}/payload"

cleanup() {
  rm -rf "${tmpdir}"
}
trap cleanup EXIT

curl -fsSL "${url}" -o "${archive}"

case "${url}" in
  *.tar.gz|*.tgz)
    tar -xzf "${archive}" -C "${install_dir}"
    ;;
  *.tar.xz)
    tar -xJf "${archive}" -C "${install_dir}"
    ;;
  *.zip)
    unzip -q "${archive}" -d "${install_dir}"
    ;;
  *)
    install -m 0755 "${archive}" "${install_dir}/${subpath}"
    ;;
esac

resolved=""
if [[ -f "${install_dir}/${subpath}" ]]; then
  resolved="${install_dir}/${subpath}"
else
  resolved="$(find "${install_dir}" -type f -path "*/${subpath}" | head -n1 || true)"
fi

if [[ -z "${resolved}" ]]; then
  echo "Could not find ${subpath} after downloading ${url}" >&2
  exit 1
fi

chmod +x "${resolved}"
ln -sf "${resolved}" "/usr/local/bin/${binary_name}"
