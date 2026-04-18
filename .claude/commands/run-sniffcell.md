# SniffCell End-to-End Analysis

You are guiding a wet-lab researcher (no command-line experience) through a complete SniffCell analysis on a long-read BAM file. Your job is to handle all technical details invisibly and communicate in plain biology language. Never show raw commands or file contents to the user — always translate into plain language.

---

## PHASE 0 — VERIFY SNIFFCELL INSTALLATION

Run `sniffcell --version` to confirm the tool is installed. If it fails, stop and tell the user:
> "SniffCell doesn't appear to be installed or isn't in your PATH. Please install it from https://github.com/Fu-Yilei/SniffCell and make sure your conda/mamba environment is activated."

Then locate the atlas directory. The atlas is a folder named `atlas/` that ships with the SniffCell repository. Find it by:

```bash
# Try to find the sniffcell package location
python3 -c "import sniffcell, os; print(os.path.dirname(sniffcell.__file__))" 2>/dev/null

# Or find it via the installed entry point
python3 -c "import importlib.util; s=importlib.util.find_spec('sniffcell'); print(s.submodule_search_locations[0] if s else 'not found')" 2>/dev/null
```

Look for an `atlas/` directory alongside the sniffcell source code. The atlas must contain:
- `all_celltypes_blocks.npy`
- `all_celltypes_blocks.index.gz`
- `all_celltypes.txt`
- `index_to_major_celltypes.json`

If you find the sniffcell source location (e.g., `/path/to/SniffCell/src/sniffcell`), the atlas is typically at `/path/to/SniffCell/atlas/`. If not found automatically, ask the user:
> "Where is your SniffCell installation folder? It should contain an `atlas/` subfolder."

Set `ATLAS_DIR=<found path>` for all subsequent steps.

Read the available tissue types from the atlas:
```bash
python3 -c "
import json
with open('${ATLAS_DIR}/index_to_major_celltypes.json') as f:
    d = json.load(f)
for k, v in d.items():
    if k != '__hierarchy__':
        print(f'{k}: {list(v.keys()) if isinstance(v, dict) else v}')
"
```

Store this tissue→cell-types mapping for later use.

---

## PHASE 1 — COLLECT USER INPUTS

Ask the user the following questions conversationally (not as a numbered list). Wait for all answers before proceeding.

1. **BAM file**: "Where is your sequencing file (BAM file)? Please paste the full path."
   - Verify the file exists: `ls -lh <BAM>`
   - Check for BAM index: `ls <BAM>.bai <BAM>.csi 2>/dev/null`. If missing, offer to create it: "I'll need to index your BAM file first — this takes about 5–15 minutes."

2. **Sample name**: Suggest the filename without extension. Let them confirm or rename.

3. **Tissue type**: Show the available tissue types from the atlas JSON in plain English. Ask them to pick. Present as a simple numbered list like:
   > "What tissue type is this sample from? Options are:
   > 1. PBMC / blood cells
   > 2. Brain (cerebellum)
   > 3. Brain (frontal cortex)
   > ... etc."
   Map their answer to the correct atlas key (`ctdmr_key`).

4. **Biological sex**: "Is this sample from a male or female donor?"

5. **Reference genome FASTA**: "Please provide the path to your reference genome FASTA file (human GRCh38/hg38 recommended for human samples)."
   - Verify it exists: `ls <REF>.fai 2>/dev/null || ls <REF>.fa.fai 2>/dev/null`. If missing, offer: "I'll need to index the reference — this takes a few minutes."

6. **Output directory**: Suggest `./sniffcell_results/<sample_name>`. Let them change it.

7. **Run label**: Auto-generate as `run_YYYYMMDD` from today's date. Show it but don't require changes.

After collecting all answers, show a confirmation:
> "Ready to analyze **[sample_name]** ([tissue]) from a [sex] donor.
> I'll look for structural variants in these cell types: [list in plain language].
> Results will go to: [output_dir]
> Estimated time: 2–5 hours depending on file size.
> Shall I start?"

Wait for "yes" before running anything.

---

## PHASE 2 — DETECT EXTERNAL TOOLS

Run this detection silently. For each tool, try `which <tool>` first. Do not bother the user unless a tool is missing.

Required tools for the full pipeline:
- `samtools` — for BAM handling
- `sniffles` — structural variant caller
- `bcftools` — VCF processing
- `bgzip` — compression
- `tabix` — indexing
- `kanpig` — SV genotyping
- `truvari` — SV merging/comparison
- `medaka` — consensus polishing
- `modkit` — methylation extraction
- `tdb` — tandem repeat database tool (may be called `tdb` or `truvari-db`)

For clair3 (SNV caller), check for `run_clair3.sh` or `clair3`:
```bash
which run_clair3.sh 2>/dev/null || which clair3 2>/dev/null
```
If clair3 is found, also look for its model directory:
```bash
# Check common model locations near the clair3 binary
CLAIR3_BIN=$(which run_clair3.sh 2>/dev/null || which clair3 2>/dev/null)
CLAIR3_MODEL_DIR=$(dirname $CLAIR3_BIN)/models 2>/dev/null
ls $CLAIR3_MODEL_DIR 2>/dev/null | head -5
```
Ask the user which clair3 model they want to use if multiple are present (show plain names like "R10.4.1 super-accuracy"), or default to the most recent `r1041_e82_400bps_sup` variant if found.

For any tool missing from PATH, ask the user:
> "I couldn't find [tool name] automatically. Do you know where it's installed? (You can type `which sniffles` in your terminal to check.)"

Build a `TOOL_PATHS` map with all discovered binaries. Store clair3 model path as `CLAIR3_MODEL`.

Also find the tandem repeat BED file for discover:
- Check if there is a `*.bed` file in the atlas directory that looks like a TR bed
- Or ask: "Do you have a tandem repeat BED file? It's used to improve structural variant calling near repetitive regions. If you're unsure, I can skip this (results may be less accurate in TR regions)."
  - If they provide one, use it. If not, omit `--tr-bed` from the discover command.

---

## PHASE 3 — BAM QUALITY CHECK

Run these checks and report results in plain language.

```bash
# Mapping statistics
samtools flagstat <BAM>

# Read length distribution (sample 10,000 reads)
samtools view <BAM> | head -10000 | awk '{print length($10)}' \
  | awk 'BEGIN{n=0;s=0;min=999999;max=0} {n++;s+=$1; if($1<min)min=$1; if($1>max)max=$1} \
         END{print "n="n, "min="min, "max="max, "mean="int(s/(n+1))}'

# Soft-clip fraction (proxy for alignment quality)
samtools view <BAM> | head -20000 | awk '{
  n=split($6, c, /[0-9]+/); m=split($6, l, /[MIDNSHPX=]/);
  sc=0; al=0;
  for(i=1;i<n;i++) { if(c[i]=="S") sc+=l[i]; else if(c[i]=="M"||c[i]=="X"||c[i]=="=") al+=l[i] }
  total_sc+=sc; total_al+=al
} END{printf "softclip_fraction=%.3f\n", total_sc/(total_sc+total_al+0.001)}'

# Verify methylation tags (MM/ML required for SniffCell)
samtools view <BAM> | head -500 | awk '{for(i=12;i<=NF;i++) if($i~/^MM:/ || $i~/^ML:/) {found=1; exit}} END{print "has_methylation_tags=" (found ? "yes" : "no")}'
```

Interpret and tell the user:
- **Mapping rate**: >90% = good; 80–90% = acceptable; <80% = warn ("Less than 80% of reads mapped — this could indicate a genome build mismatch or low data quality.")
- **Read length**: mean >5,000 bp = long-read data (expected); mean <2,000 bp = warn ("Your reads appear shorter than expected for long-read sequencing. SniffCell is designed for Oxford Nanopore or PacBio HiFi data.")
- **Soft-clip fraction**: <15% = good; 15–30% = acceptable; >30% = warn ("A high fraction of your reads have unaligned ends, which can affect variant calling accuracy.")
- **Methylation tags**: if `has_methylation_tags=no`, **STOP** and tell the user: "This BAM file doesn't contain methylation information (MM/ML tags), which SniffCell requires. Please re-run basecalling with a methylation-aware model (e.g., Dorado with a `sup` model that includes 5mC calling)." Do not proceed.

Summarize as a single plain-English sentence like:
> "Your sequencing file looks good: 94% of reads mapped, average read length 12,500 bp, methylation tags are present."

---

## PHASE 4 — FIND CELL-TYPE MARKER REGIONS (ctDMRs)

Tell the user: "First I'll identify genomic regions that are methylated differently in each cell type — these are the 'fingerprints' I'll use to sort your reads."

Generate the `split_bam_groups` string automatically from the atlas JSON for the chosen tissue key. Read the cell types:
```bash
python3 -c "
import json
with open('${ATLAS_DIR}/index_to_major_celltypes.json') as f:
    d = json.load(f)
celltypes = list(d['${CTDMR_KEY}'].keys()) if isinstance(d.get('${CTDMR_KEY}'), dict) else d.get('${CTDMR_KEY}', [])
print(','.join(celltypes))
"
```

For the `split_bam_groups` argument: default to one group per cell type (`CellType=CellType`), joined by semicolons. Exception: if the tissue key is `pbmc`, use the tested grouping `T-cell_NK-cell_B-cell=T-cell,NK-cell,B-cell;Monocyte=Monocyte` (groups lymphocytes together for better sensitivity).

Run find:
```bash
sniffcell find \
  -n ${ATLAS_DIR}/all_celltypes_blocks.npy \
  -i ${ATLAS_DIR}/all_celltypes_blocks.index.gz \
  -m ${ATLAS_DIR}/all_celltypes.txt \
  -cf ${ATLAS_DIR}/index_to_major_celltypes.json \
  -ck <ctdmr_key> \
  -o ${OUTPUT_DIR}/find/<ctdmr_key>.ctdmr.tsv \
  --diff_threshold 0.40 \
  --min_rows 2 \
  --min_cpgs 3 \
  --max_gap_bp 500
```

After completion, count the output TSV lines (`wc -l`). Set `CTDMR_TSV=${OUTPUT_DIR}/find/<ctdmr_key>.ctdmr.tsv`.

Tell the user: "Found [N] cell-type marker regions for [tissue]. [if N < 500: This is on the lower end — results may be sparser than usual. Consider whether the tissue type matches your sample exactly.]"

---

## PHASE 5 — CELL-TYPE DECONVOLUTION

Tell the user: "Now I'm sorting your sequencing reads into cell types based on their methylation patterns. This typically takes 30–90 minutes."

```bash
sniffcell deconv \
  -i <BAM> \
  -r <REFERENCE> \
  -b ${CTDMR_TSV} \
  -o ${OUTPUT_DIR}/deconv \
  -t 8 \
  --split_bam_groups "<split_bam_groups>" \
  --read_assignment_mode closest_reference_mean \
  --per_read_min_agreement 0.66
```

After completion, read `${OUTPUT_DIR}/deconv/deconv_summary.tsv`.

Get total mapped read count from flagstat (run earlier) for context.

**Interpret and report:**
- Compute total assigned reads (sum of all `n_unique_reads` in the summary).
- Assignment rate = assigned / total mapped reads.
- Good: >60%. Warn if <40%.
- Report top 3 cell types by fraction in plain English.
- Check that cell-type-specific BAM files exist in `deconv/deconv_requested_group_splits/` with non-zero size.

Tell the user:
> "Cell sorting complete. [N] ([X]%) of your sequencing reads were assigned to a cell type. Here's the breakdown:
> - [Cell type A]: [X]% of reads
> - [Cell type B]: [X]% of reads
> - [etc.]
> This looks [consistent / unusual] for a [tissue] sample."

Set `DECONV_DIR=${OUTPUT_DIR}/deconv`.

---

## PHASE 6 — STRUCTURAL VARIANT DISCOVERY

Tell the user: "Now I'm calling structural variants (insertions, deletions, etc.) in each cell-type subset. This typically takes 1–3 hours."

Build the discover command with all detected tool paths. Only include `--tr-bed` if a TR BED was provided. Only include clair3 arguments if clair3 was found.

```bash
sniffcell discover tools run \
  --deconv-dir ${DECONV_DIR} \
  --reference <REFERENCE> \
  [--tr-bed <TR_BED> if available] \
  --sex <sex> \
  --scheduler local \
  --sample-id <sample_id> \
  --run-id <run_id> \
  --stages all \
  --rerun-failed \
  --threads 8 \
  --sniffles-bin <SNIFFLES_BIN> \
  --bcftools-bin <BCFTOOLS_BIN> \
  --bgzip-bin <BGZIP_BIN> \
  --kanpig-bin <KANPIG_BIN> \
  --truvari-bin <TRUVARI_BIN> \
  --medaka-bin <MEDAKA_BIN> \
  --tdb-bin <TDB_BIN> \
  --modkit-bin <MODKIT_BIN> \
  --tabix-bin <TABIX_BIN> \
  --sniffles-mosaic-filter-expression 'INFO/MOSAIC=1' \
  --sniffles-cluster-merge-len 0.2 \
  --kanpig-seqsim 0.8 \
  --kanpig-sizesim 0.85 \
  --kanpig-passonly \
  --kanpig-sample-name-template '{sample_id}_{group}' \
  --truvari-refdist 500 \
  --truvari-pctseq 0.95 \
  --truvari-pctsize 0.95 \
  --truvari-passonly \
  --medaka-model 'dna_r10.4.1_e8.2_400bps_sup@v4.3.0:consensus' \
  --medaka-padding 250 \
  --medaka-sample-name-template '{sample_id}.{group}' \
  --tdb-create-mem 4 \
  --mods-mode separate \
  [--clair3-bin <CLAIR3_BIN> --clair3-platform ont --clair3-model-path <CLAIR3_MODEL> if clair3 available]
```

After completion, set `DISCOVER_DIR=${DECONV_DIR}/deconv_requested_group_splits/discover/<run_id>`.

Check `run_summary.json` for completed stages and `harmonized_variants.tsv` line count.

Tell the user: "Variant discovery complete. Found [N] candidate structural variants across your cell-type-specific data."

Set `HARMONIZED_VCF=${DISCOVER_DIR}/harmonized_variants.tsv`.

---

## PHASE 7 — CELL-TYPE ANNOTATION OF VARIANTS

Tell the user: "Now I'm linking each structural variant to the cell type whose methylation pattern best supports it."

```bash
sniffcell anno \
  --input <BAM> \
  --bed ${CTDMR_TSV} \
  --vcf ${HARMONIZED_VCF} \
  --reference <REFERENCE> \
  --output ${OUTPUT_DIR}/anno \
  --threads 8
```

After completion, check `${OUTPUT_DIR}/anno/sv_assignment.tsv`:

```bash
# Total SVs
TOTAL=$(tail -n +2 ${OUTPUT_DIR}/anno/sv_assignment.tsv | wc -l)

# Assigned (assigned_code column not empty — column 8)
ASSIGNED=$(awk -F'\t' 'NR>1 && $8!="" {c++} END{print c+0}' ${OUTPUT_DIR}/anno/sv_assignment.tsv)

# High-confidence: overlap_pct == 1.0 and majority_pct >= 0.8
HIGH_CONF=$(awk -F'\t' 'NR>1 && $5=="1.0" && $7+0>=0.8 {c++} END{print c+0}' ${OUTPUT_DIR}/anno/sv_assignment.tsv)

# Hard conflicts
CONFLICTS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="has_hard_conflict") col=i}
                         NR>1 && $col=="True" {n++} END{print n+0}' ${OUTPUT_DIR}/anno/sv_assignment.tsv)

echo "total=$TOTAL assigned=$ASSIGNED high_conf=$HIGH_CONF conflicts=$CONFLICTS"
```

Tell the user:
> "Annotation complete. Out of [total] structural variants, [N] ([%]) were assigned to a specific cell type. [H] of those are high-confidence. [C] showed conflicting signals."

Warn if assignment rate <20%: "Fewer variants than expected were assigned. This can happen with low sequencing depth or if the sample doesn't cleanly match the reference atlas. You can try re-running annotation with more relaxed thresholds — just ask me."

---

## PHASE 8 — GENERATE REPORT

Tell the user: "Generating your results report..."

```bash
sniffcell report \
  --anno_output ${OUTPUT_DIR}/anno \
  --min_overlap_pct 0.8 \
  --min_majority_pct 0.8 \
  --with_figures \
  --figure_threads 4
```

Check that `${OUTPUT_DIR}/anno/report/index.html` exists and `high_confidence_sv.tsv` has data rows.

---

## PHASE 9 — INTERPRET RESULTS FOR THE USER

This is the most important step. Read the key output files and explain everything in plain biology language.

### 9a. BAM evidence inspection for top SVs

Read `${OUTPUT_DIR}/anno/report/high_confidence_sv.tsv`. For the top 5 high-confidence SVs (or all if fewer than 5), extract coordinates from `sv_assignment_readable.tsv` and run CIGAR-level inspection:

```bash
# For each top SV at chrom:start-end:
samtools view <BAM> <chrom:start-end> | awk '{
  cigar=$6
  ins=0; del=0; sc=0; al=0; reads++
  while(match(cigar, /([0-9]+)([MIDSH])/, arr)) {
    len=arr[1]+0; op=arr[2]
    if(op=="I") ins+=len
    else if(op=="D") del+=len
    else if(op=="S") sc+=len
    else if(op=="M") al+=len
    cigar=substr(cigar, RSTART+RLENGTH)
  }
  total_ins+=ins; total_del+=del; total_sc+=sc; total_al+=al
} END{
  print "reads="reads
  print "avg_ins_bp=" (reads>0 ? int(total_ins/reads) : 0)
  print "avg_del_bp=" (reads>0 ? int(total_del/reads) : 0)
  print "softclip_frac=" (reads>0 ? int(100*total_sc/(total_sc+total_al+1)) "%" : "N/A")
}'

# Also check how many reads support the variant (carry the insertion/deletion signature)
samtools view <BAM> <chrom:start-end> | awk -v size=<SV_SIZE_BP> '
  {cigar=$6; has_event=0
   while(match(cigar, /([0-9]+)([ID])/, arr)) {
     if(arr[1]+0 >= size*0.7) has_event=1
     cigar=substr(cigar, RSTART+RLENGTH)
   }
   total++; if(has_event) support++
  }
  END{printf "support_fraction=%d/%d (%.0f%%)\n", support, total, 100*support/(total+0.001)}'
```

Interpret for the user:
- If avg_ins_bp or avg_del_bp is close to the SV size: reads carrying the variant are clearly visible — "strong CIGAR support."
- If softclip_frac at the locus is >40%: possible alignment artifact — "some reads had trouble aligning here, which is common near complex variants."
- If support_fraction <10%: "This variant appears in a small minority of reads at this location (mosaic or low-frequency variant)."

### 9b. Final plain-language summary

Produce a structured plain-English summary covering these sections. Use simple headers so the user can scan it.

**Your sequencing data quality**
Summarize the BAM QC results: mapped reads, read length, methylation tag status.

**Cell type composition**
Describe what the deconvolution found: which cell types are present, their proportions, whether this matches expectations for the tissue.

**Structural variants found**
- Total variants called
- How many are assigned to specific cell types and which ones
- How many are shared across cell types (these may be germline variants)
- How many are high-confidence

**Your most reliable findings** (high-confidence SVs)
List the top 5 with:
- Variant ID (as it will appear in the report)
- Location in plain terms (chromosome, rough position)
- Cell type it was assigned to
- Size and type (e.g., "a 3,500 bp deletion in chromosome 7")
- CIGAR support quality

**What good output looks like** — always include this section
> - **Assignment rate**: 50–80% of variants assigned is typical. Above 80% is excellent. Below 20% suggests a tissue type mismatch or low coverage.
> - **High-confidence variants**: these have 100% of their supporting reads from one cell type and strong methylation signal. These are your most actionable findings.
> - **Shared/germline variants**: variants seen equally across all cell types are likely inherited (present in all your cells) rather than somatic.
> - **Hard conflicts**: variants where different reads point to different cell types — these should be treated with caution. A small number (<10% of assigned) is normal.
> - **Coverage note**: SniffCell works best with >30x sequencing depth. Below 20x, variant calls may be incomplete.

**Where to find your results**
- Full interactive report: `[path]/anno/report/index.html` — open in a web browser
- High-confidence variant table: `[path]/anno/report/high_confidence_sv.tsv`
- All variant assignments: `[path]/anno/sv_assignment_readable.tsv`

**What to do next**
- "To see a visualization of any specific variant, just ask me and provide the variant ID from the report."
- "To check if any variants overlap known disease genes, ask me to cross-reference with a gene list."
- "If you want to try different assignment thresholds (more or fewer variants), I can re-run the annotation step without starting from scratch."

---

## ERROR HANDLING RULES

- **If any step fails**: read the error output, diagnose, and explain in plain language. Common causes:
  - Missing BAM index → create with `samtools index <BAM>`
  - Reference FASTA not indexed → index with `samtools faidx <REF>` and `bwa index <REF>` if needed
  - Tool not found → check PATH, ask user for binary location
  - Disk full → check with `df -h <OUTPUT_DIR>` and report remaining space
  - Low per-cell-type read count → explain that a cell type BAM may have too few reads for variant calling (typical minimum ~5x)
  - Out of memory → suggest reducing `--threads` or running on a machine with more RAM

- **Never silently skip a failed step** — always report what failed and offer options.
- **Between long steps**, give a progress update using elapsed time from the previous step.
