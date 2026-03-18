# Deconv Postprocess Design

## Goal

Add a new `sniffcell postprocess` command that operates on finished `sniffcell deconv` output, discovers the two split BAMs in `deconv_requested_group_splits/`, and runs downstream SV, tandem-repeat, and methylation postprocessing in either:

- local sequential mode
- Slurm script generation / submission mode

This command is specifically for the "two split BAMs from deconv" workflow.

## Why `postprocess` Instead Of `deconv-postprocess`

The current CLI uses one flat layer of subcommands in [parse_args.py](/users/u254106/130/SniffMeth/src/sniffcell/parse_args.py) and [main.py](/users/u254106/130/SniffMeth/src/sniffcell/main.py). The cleanest implementation is therefore a new top-level command:

```bash
sniffcell postprocess ...
```

instead of a nested `sniffcell deconv postprocess ...` CLI.

## Scope

The command should support four branches:

1. Sniffles on two split BAMs
2. Medaka tandem on two split BAMs
3. Automatic variant postprocessing
4. CpG 5mC / 5hmC extraction on the split BAMs

Branch 3 means:

- for Sniffles outputs: Kanpig mosaic + Truvari collapse
- for Medaka outputs: `tdb create` + `tdb merge` + optional custom follow-up script

## Command Shape

```bash
sniffcell postprocess \
  --deconv-dir /path/to/sample/deconv \
  --reference /path/to/ref.fa \
  --tr-bed /path/to/adotto.v2.trgt.lite.bed \
  --sex male \
  --scheduler local|slurm \
  [stage and tool arguments...]
```

### Required Arguments

- `--deconv-dir`
- `--reference`
- `--tr-bed`
- `--sex`
- `--scheduler`

### Core Optional Arguments

- `--split-dir`
  Default: `<deconv-dir>/deconv_requested_group_splits`
- `--sample-id`
  Optional override; otherwise infer from directory layout
- `--groups`
  Optional comma-separated override for the two requested groups
- `--stages`
  Comma-separated subset, for example:
  - `sv`
  - `sv,medaka`
  - `medaka,tdb,mods`
  - `sniffles,kanpig,collapse`
- `--dry-run`
- `--force`
- `--rerun-failed`
- `--keep-going`

## Stage Model

The implementation should use task instances, not only coarse stages.

For a 2-group deconv output, build this DAG:

- `discover`
- `sniffles[groupA]`
- `sniffles[groupB]`
- `sniffles_filter[groupA]`
- `sniffles_filter[groupB]`
- `kanpig[groupA]`
- `kanpig[groupB]`
- `truvari_collapse[groupA,groupB]`
- `medaka[groupA]`
- `medaka[groupB]`
- `tdb_create[groupA]`
- `tdb_create[groupB]`
- `tdb_merge[sample]`
- `custom_post_tdb[sample]` optional
- `modkit[groupA]`
- `modkit[groupB]`
- `finalize`

Dependencies:

- `discover -> all downstream tasks`
- `sniffles[groupX] -> sniffles_filter[groupX] -> kanpig[groupX]`
- `kanpig[groupA] + kanpig[groupB] -> truvari_collapse`
- `medaka[groupX] -> tdb_create[groupX]`
- `tdb_create[groupA] + tdb_create[groupB] -> tdb_merge`
- `tdb_merge -> custom_post_tdb` when enabled
- `modkit[groupX]` depends only on `discover`
- `finalize` depends on:
  - `truvari_collapse`
  - `tdb_merge` or `custom_post_tdb`
  - `modkit[groupA]`
  - `modkit[groupB]`

SV and methylation branches should remain independent after discovery.

## Discovery Contract

Discovery should read:

- `<split-dir>/requested_group_splits.tsv`

and resolve:

- exactly two requested groups
- each BAM path
- each BAM index
- output root
- sample ID

The command should fail early if:

- the split manifest is missing
- fewer or more than two groups are found
- a BAM or `.bai` file is missing

## Sniffles Branch

### Required Behavior

Sniffles should run once per split BAM and should default to:

- `--mosaic`
- `--mosaic-include-germline`
- `--output-rnames`

These should be on by default and should not require the user to pass them each time.

### Exposed Arguments

- `--sniffles-bin`
- `--sniffles-threads`
- `--sniffles-output-rnames`
  Default: on
- `--sniffles-mosaic`
  Default: on
- `--sniffles-include-germline`
  Default: on

### Command Template

```bash
sniffles \
  --input <group.bam> \
  --reference <ref.fa> \
  --vcf <group>.sniffles.vcf.gz \
  --snf <group>.snf \
  --threads <N> \
  --mosaic \
  --mosaic-include-germline \
  --output-rnames
```

### Outputs

Under:

`<split-dir>/postprocess/<run_id>/sv/sniffles/<group>/`

write:

- `sniffles.raw.vcf.gz`
- `sniffles.raw.vcf.gz.tbi`
- `sniffles.raw.snf`
- `sniffles.command.txt`
- `sniffles.done.json`

## Sniffles Mosaic Filter

Kanpig expects a filtered mosaic-oriented input VCF. This should be an explicit stage.

### Exposed Arguments

- `--sniffles-mosaic-filter-expression`
- `--bcftools-bin`

### Behavior

For each group, derive:

- `sniffles.mosaic_only.vcf.gz`
- `sniffles.mosaic_only.vcf.gz.tbi`

This stage is the only deliberately open design point right now because the exact filter expression has not been frozen yet.

## Kanpig Branch

### Exposed Arguments

- `--kanpig-bin`
- `--kanpig-threads`
- `--kanpig-seqsim`
  Default: `0.8`
- `--kanpig-sizesim`
  Default: `0.85`
- `--kanpig-passonly`
  Default: on
- `--kanpig-sample-name-template`
  Default: `{sample_id}_{group}`

### Command Template

```bash
kanpig mosaic \
  --input <sniffles.mosaic_only.vcf.gz> \
  --reference <ref.fa> \
  --reads <group.bam> \
  --sample <sample_name> \
  --passonly \
  --seqsim <seqsim> \
  --sizesim <sizesim> \
  -t <threads> \
  --rnames <kanpig.rnames.tsv> \
  -o <kanpig.mosaic.vcf>
```

### Outputs

Under:

`<split-dir>/postprocess/<run_id>/sv/kanpig/<group>/`

write:

- `kanpig.mosaic.vcf.gz`
- `kanpig.mosaic.vcf.gz.tbi`
- `kanpig.rnames.tsv`
- `kanpig.command.txt`
- `kanpig.done.json`

## Truvari Collapse

Truvari collapse should run after both group-level Kanpig outputs are complete.

### Exposed Arguments

- `--truvari-bin`
- `--truvari-refdist`
- `--truvari-pctseq`
- `--truvari-pctsize`
- `--collapse-use`
  Allowed values:
  - `kanpig`
  - `sniffles`
  Default: `kanpig`

### Behavior

By default, collapse the two Kanpig-refined VCFs rather than the raw Sniffles VCFs.

### Outputs

Under:

`<split-dir>/postprocess/<run_id>/sv/truvari_collapse/<groupA>_vs_<groupB>/`

write:

- `collapsed.vcf.gz`
- `collapse.report.json`
- `collapse.command.txt`
- `collapse.done.json`

## Medaka Tandem Branch

### Important CLI Fact

In the installed Medaka version used here, tandem-repeat calling uses:

- `--workers`

not `--threads`.

### Exposed Arguments

- `--medaka-bin`
- `--medaka-workers`
- `--medaka-model`
  Default: `dna_r10.4.1_e8.2_400bps_sup@v4.3.0:consensus`
- `--medaka-padding`
  Default: `250`
- `--medaka-sample-name-template`
  Default: `{sample_id}.{group}`

### Command Template

```bash
medaka tandem \
  --workers <N> \
  --model <model> \
  --sample_name <sample_name> \
  --padding <padding> \
  <group.bam> \
  <ref.fna> \
  <tr.bed> \
  <sex> \
  <group_output_dir>
```

### Outputs

Under:

`<split-dir>/postprocess/<run_id>/medaka_tandem/<group>.medaka/`

expect at least:

- `medaka_to_ref.TR.vcf`
- `consensus.fasta`
- `medaka_to_ref.bam`
- `medaka_to_ref.bam.bai`
- `trimmed_reads_to_poa.bam`
- `trimmed_reads_to_poa.bam.bai`
- `medaka.command.txt`
- `medaka.done.json`

## TDB Create And Merge

TDB creation and merging should be explicit stages, even if a Medaka run already emits a `.tdb`.

### Exposed Arguments

- `--tdb-bin`
- `--tdb-create-mem`
- `--tdb-create-force`
- `--tdb-merge-threads`
- `--skip-tdb-create-if-present`
- `--skip-tdb-merge-if-present`

There should be no `--tdb-merge-mem` parameter.

### TDB Create Template

```bash
tdb create \
  -o <group>.tdb \
  --mem <GB> \
  [--force] \
  <group>.medaka/medaka_to_ref.TR.vcf
```

### TDB Merge Template

```bash
tdb merge \
  -o <sample>.merged.tdb \
  --threads <N> \
  <groupA>.tdb \
  <groupB>.tdb
```

### Outputs

Under:

`<split-dir>/postprocess/<run_id>/medaka_tandem/tdb/`

write:

- `<group>.tdb/`
- `<sample>.merged.tdb/`
- `tdb_create.<group>.command.txt`
- `tdb_merge.command.txt`
- `tdb_create.<group>.done.json`
- `tdb_merge.done.json`

## Custom Post-TDB Hook

The design should support a user-provided follow-up script.

### Exposed Arguments

- `--post-tdb-script`
- `--post-tdb-args`

### Invocation Contract

Pass at least:

- merged TDB path
- per-group TDB paths
- sample ID
- group names
- output root

This keeps custom logic out of the main module.

## 5mC / 5hmC CpG Branch

Use `modkit` for this branch.

### Exposed Arguments

- `--modkit-bin`
- `--modkit-threads`
- `--emit-read-level-mods`
- `--mods-mode`
  Allowed values:
  - `separate`
  - `combined`
  Default: `separate`

### Site-Level Command Template

```bash
modkit pileup \
  <group.bam> \
  <group>.cpg.bedmethyl.gz \
  --cpg \
  --ref <ref.fa> \
  --modified-bases 5mC 5hmC \
  --bgzf \
  -t <threads>
```

### Optional Read-Level Command Template

```bash
modkit extract calls \
  <group.bam> \
  <group>.cpg.calls.tsv.gz \
  --reference <ref.fa> \
  --cpg \
  --bgzf \
  --pass-only \
  -t <threads>
```

### Output Policy

Keep 5mC and 5hmC separate first. If a combined track is needed later, derive it in Python from the separate outputs.

### Outputs

Under:

`<split-dir>/postprocess/<run_id>/modkit/<group>/`

write:

- `<group>.cpg.bedmethyl.gz`
- `<group>.cpg.bedmethyl.gz.tbi` if generated
- optionally `<group>.cpg.calls.tsv.gz`
- `modkit.command.txt`
- `modkit.done.json`

## Scheduler Modes

## Local Mode

Local mode should execute the DAG topologically in-process.

Behavior:

- stop on first failure by default
- optionally continue independent branches with `--keep-going`
- support `--stages`
- support `--groups`
- support `--rerun-failed`
- support `--force`

## Slurm Mode

Slurm mode should support:

- `--slurm-mode render`
- `--slurm-mode submit`

Recommended job splitting:

- `discover`: one sample-level job
- `sniffles`: one task per group
- `sniffles_filter`: one task per group
- `kanpig`: one task per group
- `truvari_collapse`: one sample-level job
- `medaka`: one task per group
- `tdb_create`: one task per group
- `tdb_merge`: one sample-level job
- `custom_post_tdb`: one sample-level job
- `modkit`: one task per group
- `finalize`: one sample-level job

Recommended scripts for a two-group sample:

- `01_sniffles.array.sbatch`
- `02_sniffles_filter.array.sbatch`
- `03_kanpig.array.sbatch`
- `04_truvari_collapse.sbatch`
- `05_medaka.array.sbatch`
- `06_tdb_create.array.sbatch`
- `07_tdb_merge.sbatch`
- `08_custom_post_tdb.sbatch`
- `09_modkit.array.sbatch`
- `10_finalize.sbatch`

Use `afterok:<jobid>` dependencies between stage families. Do not overcomplicate this with per-index dependency management unless later support for many-group splitting requires it.

## Run Directory And Manifests

Each run should write to:

`<split-dir>/postprocess/<run_id>/`

Inside:

- `manifest/`
- `logs/`
- `commands/`
- `status/`
- `slurm/`
- `sv/`
- `medaka_tandem/`
- `modkit/`

### Required Manifests

1. `manifest/postprocess_run_manifest.json`

Contains:

- sample ID
- deconv dir
- split dir
- split manifest path
- resolved groups
- BAM paths
- tool paths
- scheduler mode
- enabled stages
- parameters
- run ID
- timestamps

2. `manifest/postprocess_task_manifest.tsv`

One row per task with:

- `task_id`
- `stage`
- `scope`
- `sample_id`
- `group_name`
- `deps`
- `inputs`
- `outputs`
- `resources_json`
- `command_json`
- `status`
- `job_script`
- `slurm_job_id`

3. `status/postprocess_status.json`

Live task states:

- `pending`
- `running`
- `completed`
- `failed`
- `skipped`

with timestamps, exit code, stdout, stderr, output validation, and command hash.

## Resource Model

Each task should carry scheduler-agnostic resources:

- `threads`
- `memory_gb` optional
- `time_h`
- `partition` optional
- `account` optional

Local mode only uses `threads`.
Slurm mode maps these to `#SBATCH` directives.

## Restart And Idempotency Rules

### Completion

A task counts as complete only if:

- exit code is `0`
- expected outputs exist
- outputs are non-empty where appropriate
- done marker JSON has been written

### Default Skip Logic

If a task has a valid done marker and matching outputs, skip it.

### Rerun Cases

- `--rerun-failed`
  rerun failed tasks only
- `--force`
  rebuild selected tasks even if complete
- `--stages`
  rerun only a subset of branches

### Command Drift

If outputs exist but the effective command hash changes, require `--force` or write to a new run directory.

## Proposed Python Module Layout

Add:

- `src/sniffcell/postprocess/__init__.py`
- `src/sniffcell/postprocess/postprocess.py`
- `src/sniffcell/postprocess/config.py`
- `src/sniffcell/postprocess/discovery.py`
- `src/sniffcell/postprocess/dag.py`
- `src/sniffcell/postprocess/tasks.py`
- `src/sniffcell/postprocess/runner_local.py`
- `src/sniffcell/postprocess/runner_slurm.py`
- `src/sniffcell/postprocess/manifest.py`
- `src/sniffcell/postprocess/logging_utils.py`
- `src/sniffcell/postprocess/stages/sv.py`
- `src/sniffcell/postprocess/stages/tr.py`
- `src/sniffcell/postprocess/stages/mods.py`

Wire the command in:

- [parse_args.py](/users/u254106/130/SniffMeth/src/sniffcell/parse_args.py)
- [main.py](/users/u254106/130/SniffMeth/src/sniffcell/main.py)

## Recommended First Implementation Order

1. CLI parser and config normalization
2. discovery + run manifest
3. local runner
4. Sniffles stage
5. Sniffles filter stage
6. Kanpig stage
7. Truvari collapse stage
8. Medaka stage
9. TDB create / merge
10. Modkit stage
11. Slurm renderer / submitter
12. custom post-TDB hook

## Open Item

The one design item that still needs a final rule is the exact Sniffles mosaic-only filter expression used before Kanpig. Everything else in this design has a concrete precedent in the current workspace or installed tools.
