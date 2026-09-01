#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: $0 /path/to/GRCh38_no_alt.fa [output_dir] [threads]" >&2
    exit 2
fi

REFERENCE=$(realpath "$1")
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUTPUT_DIR=${2:-"${SCRIPT_DIR}/output"}
THREADS=${3:-4}
INPUT_DIR="${SCRIPT_DIR}/inputs"

if [[ ! -f "${REFERENCE}" || ! -f "${REFERENCE}.fai" ]]; then
    echo "Reference and index are required: ${REFERENCE} and ${REFERENCE}.fai" >&2
    exit 2
fi

CHR2_LENGTH=$(awk '$1 == "chr2" {print $2; exit}' "${REFERENCE}.fai")
if [[ "${CHR2_LENGTH}" != "242193529" ]]; then
    echo "Expected GRCh38 no-alt chr2 length 242193529; found ${CHR2_LENGTH:-missing}." >&2
    exit 2
fi

mkdir -p "${OUTPUT_DIR}"

sniffcell find \
    -n "${INPUT_DIR}/atlas.npy" \
    -i "${INPUT_DIR}/atlas.index.tsv" \
    -m "${INPUT_DIR}/atlas.samples.txt" \
    -cf "${INPUT_DIR}/celltypes.json" \
    -ck brain_cereb \
    -o "${OUTPUT_DIR}/brain_cereb.ctdmr.tsv" \
    --diff_threshold 0.35 \
    --min_rows 1

sniffcell deconv \
    -i "${INPUT_DIR}/SH3RF3_example.bam" \
    -r "${REFERENCE}" \
    -b "${OUTPUT_DIR}/brain_cereb.ctdmr.tsv" \
    -o "${OUTPUT_DIR}/deconv" \
    --regions "${INPUT_DIR}/target.bed" \
    --regions-ctdmrs 4 \
    --split_bam_groups 'Neuron=Neuron;Oligodendrocyte=Oligodendrocyte' \
    --skip_overall_summary \
    -t "${THREADS}"

sniffcell discover tools run \
    --deconv-dir "${OUTPUT_DIR}/deconv" \
    --reference "${REFERENCE}" \
    --tr-bed "${INPUT_DIR}/tr_catalog.bed" \
    --sex female \
    --sample-id SH3RF3_example \
    --run-id dual_sv_tr \
    --stages all \
    --sniffles-cluster-merge-len 0.31 \
    --threads "${THREADS}"

HARMONIZED="${OUTPUT_DIR}/deconv/deconv_requested_group_splits/discover/dual_sv_tr/harmonized_variants.tsv"

sniffcell anno \
    -i "${INPUT_DIR}/SH3RF3_example.bam" \
    -v "${HARMONIZED}" \
    -r "${REFERENCE}" \
    -b "${OUTPUT_DIR}/brain_cereb.ctdmr.tsv" \
    -o "${OUTPUT_DIR}/anno" \
    --deconv-reads "${OUTPUT_DIR}/deconv/deconv_reads_classification.tsv" \
    --evidence_mode per_read \
    -w 10000 \
    -t "${THREADS}"

sniffcell viz \
    --anno_output "${OUTPUT_DIR}/anno" \
    -s chr2_109199301_109199876 \
    -f png \
    --export_tables \
    -o "${OUTPUT_DIR}/anno/SH3RF3_TR_expansion.png"

sniffcell report \
    --anno_output "${OUTPUT_DIR}/anno" \
    --include_unassigned \
    --with_figures \
    --figure_threads "${THREADS}" \
    -o "${OUTPUT_DIR}/anno/report"

python "${SCRIPT_DIR}/validate_outputs.py" "${OUTPUT_DIR}"
