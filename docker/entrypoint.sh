#!/usr/bin/env bash
set -euo pipefail

if [[ $# -eq 0 ]]; then
  exec sniffcell -h
fi

case "$1" in
  -h|--help|-v|--version|find|deconv|anno|svanno|dmsv|viz|igvviz|report|discover)
    exec sniffcell "$@"
    ;;
  sniffcell|sniffcell-check-discover|sniffcell-discover-sv|python|bash|sh|bcftools|bgzip|tabix|sniffles|kanpig|modkit|medaka|run_clair3.sh|run_clairs)
    exec "$@"
    ;;
  *)
    exec "$@"
    ;;
esac
