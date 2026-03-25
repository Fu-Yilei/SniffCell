#!/usr/bin/env bash
set -euo pipefail

runtime=""
if command -v docker >/dev/null 2>&1; then
  runtime="docker"
elif command -v podman >/dev/null 2>&1; then
  runtime="podman"
fi

checks_total=0
checks_passed=0
checks_failed=0

record_check() {
  local name="$1"
  local code="$2"
  local detail="$3"
  checks_total=$((checks_total + 1))
  if [[ "${code}" -eq 0 ]]; then
    checks_passed=$((checks_passed + 1))
    printf 'OK      %s -> %s\n' "${name}" "${detail}"
  else
    checks_failed=$((checks_failed + 1))
    printf 'FAIL    %s -> %s\n' "${name}" "${detail}"
  fi
}

if [[ -n "${runtime}" ]]; then
  record_check "runtime_present" 0 "${runtime}"
else
  record_check "runtime_present" 1 "neither docker nor podman is installed"
  printf 'SUMMARY passed=%d failed=%d total=%d\n' "${checks_passed}" "${checks_failed}" "${checks_total}"
  exit 1
fi

if "${runtime}" info >/dev/null 2>&1; then
  record_check "runtime_info" 0 "${runtime} info succeeded"
else
  record_check "runtime_info" 1 "${runtime} info failed"
fi

if [[ "${runtime}" == "podman" ]]; then
  if grep -q "^${USER}:" /etc/subuid 2>/dev/null && grep -q "^${USER}:" /etc/subgid 2>/dev/null; then
    record_check "rootless_idmap" 0 "subuid/subgid entries found for ${USER}"
  else
    record_check "rootless_idmap" 1 "missing subuid/subgid entries for ${USER}"
  fi
fi

base_image="docker.io/mambaorg/micromamba:latest"
if "${runtime}" pull "${base_image}" >/tmp/sniffcell_builder_pull.log 2>&1; then
  record_check "base_pull" 0 "${base_image}"
else
  detail="$(tail -n 3 /tmp/sniffcell_builder_pull.log | tr '\n' ' ' | sed 's/[[:space:]]\\+/ /g')"
  record_check "base_pull" 1 "${detail}"
fi
rm -f /tmp/sniffcell_builder_pull.log

printf 'SUMMARY passed=%d failed=%d total=%d\n' "${checks_passed}" "${checks_failed}" "${checks_total}"
if [[ "${checks_failed}" -ne 0 ]]; then
  exit 1
fi
