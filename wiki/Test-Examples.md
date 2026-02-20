# Test Examples

This page focuses on practical validation runs for `sniffcell` in real workflows.

## 1. Quick Regression Tests (Developer)

Run the full unit test suite:

```bash
python -m unittest discover -s tests -v
```

Run only SV assignment tests:

```bash
python -m unittest -v tests/test_variant_assignment.py
```

What this covers:
- code schema decoding (`code_order` + binary codes)
- vote aggregation from read-level evidence
- hard conflict handling (`has_hard_conflict`)
- readable report generation

## 2. End-to-End Smoke Test (CLI)

This sequence validates `find -> anno -> svanno -> report -> viz -> dmsv`.

### 2.1 Generate ctDMRs (`find`)

```bash
sniffcell find \
  -n atlas/all_celltypes_blocks.npy \
  -i atlas/all_celltypes_blocks.index.gz \
  -cf atlas/index_to_major_celltypes.json \
  -m atlas/all_celltypes.txt \
  -ck pbmc \
  -o out/pbmc_ctdmr.tsv \
  --diff_threshold 0.40 \
  --min_rows 2 \
  --min_cpgs 3 \
  --max_gap_bp 500
```

Checks:
- `out/pbmc_ctdmr.tsv` exists and has rows
- `out/pbmc_ctdmr.tsv.igv.bed` exists
- output has columns like `best_group_leaves`, `other_group_leaves`, `code_order`

### 2.2 Annotate SV-supporting reads (`anno`)

```bash
sniffcell anno \
  -i sample.bam \
  -v sample.vcf.gz \
  -r ref.fa \
  -b out/pbmc_ctdmr.tsv \
  -o out/anno \
  -w 10000 \
  -t 8 \
  --evidence_mode all_rows \
  --min_overlap_pct 0.0 \
  --min_agreement_pct 1.0
```

Checks:
- `out/anno/reads_classification.tsv`
- `out/anno/blocks_classification.tsv`
- `out/anno/sv_assignment.tsv`
- `out/anno/sv_assignment_readable.tsv`
- `out/anno/anno_run_manifest.json`

### 2.3 Recompute assignment from existing read table (`svanno`)

Use this to test threshold sensitivity without rerunning `anno` clustering:

```bash
sniffcell svanno \
  -v sample.vcf.gz \
  -i out/anno/reads_classification.tsv \
  -o out/anno_relaxed \
  -w 10000 \
  --evidence_mode all_rows \
  --min_overlap_pct 0.0 \
  --min_agreement_pct 0.6
```

Checks:
- compare `out/anno/sv_assignment.tsv` vs `out/anno_relaxed/sv_assignment.tsv`
- ensure lower `min_agreement_pct` increases assigned SV count (expected in noisy datasets)

### 2.4 Render one SV panel (`viz`)

```bash
sniffcell viz \
  --anno_output out/anno \
  -s sniffles.SV123 \
  -f png \
  --export_tables
```

Checks:
- figure exists: `out/anno/sniffles.SV123.viz.png`
- table exports exist:
  - `*.summary.tsv`
  - `*.supporting_reads_assignment.tsv`
  - `*.supporting_reads_ctdmr_methylation.tsv`

### 2.5 Build high-confidence SV HTML report (`report`)

```bash
sniffcell report \
  --anno_output out/anno \
  --min_overlap_pct 0.5 \
  --min_majority_pct 0.95 \
  -f png
```

Checks:
- `out/anno/report/index.html`
- `out/anno/report/high_confidence_sv.tsv`
- `out/anno/report/figures/*.viz.png`
- if many SVs pass, expect report generation to take longer (one `viz` render per SV)

### 2.6 Differential methylation near SVs (`dmsv`)

```bash
sniffcell dmsv \
  -i sample.bam \
  -v sample.vcf.gz \
  -r ref.fa \
  -o out/dmsv \
  -m 3 \
  -f 1000 \
  -c 5 \
  -t 8
```

Checks:
- `out/dmsv/significant_SVs.tsv`
- `out/dmsv/sv_details/<sv_id>.tsv.gz`

## 3. Practical QA Queries

Use these after a run:

```bash
# How many SVs received any assignment?
awk -F'\t' 'NR>1 && $7!="" {c++} END{print c+0}' out/anno/sv_assignment.tsv

# How many SVs are marked hard conflict?
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="has_hard_conflict") c=i}
            NR>1 && $c=="True"{n++} END{print n+0}' out/anno/sv_assignment.tsv

# Top linked cell types in readable report
cut -f11 out/anno/sv_assignment_readable.tsv | tail -n +2 | tr '|' '\n' | sort | uniq -c | sort -nr | head
```

## 4. Common Failure Patterns

- Empty `reads_classification.tsv`:
  - likely no ctDMRs passed window filtering around SVs
  - increase `-w/--window` or inspect ctDMR genomic distribution
- Many SVs with empty `assigned_code`:
  - strict default is `--min_agreement_pct 1.0`
  - try `svanno` sweep at `0.6` to inspect sensitivity
- `viz` missing reference-driven methylation tables:
  - make sure `reference` exists in `anno_run_manifest.json` or pass `-r`
