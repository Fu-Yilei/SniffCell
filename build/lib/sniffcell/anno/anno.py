import os
import pandas as pd
from sniffcell.anno.kmeans import kmeans_cluster_cells
from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam
from sniffcell.anno.filter_bed_based_on_variants import filter_bed_based_on_variants
from sniffcell.anno.vcf_to_df import read_vcf_to_df
from sniffcell.anno.variant_assignment import assign_sv_celltypes
from tqdm import tqdm
import multiprocessing as mp
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(processName)s] %(levelname)s: %(message)s")

def _one_dmr(args):
    # args: (row_dict, input_bam, reference_fasta)
    logger = logging.getLogger("anno._one_dmr")

    row, input_file, reference = args
    chrom = str(row["chr"])
    start = int(row["start"])
    end   = int(row["end"])

    best_group = str(row["best_group"])
    best_dir   = row.get("best_dir", None)

    # collect arbitrary cell types from mean_* columns in THIS row
    cell_types = []
    for k, v in row.items():
        if isinstance(k, str) and k.startswith("mean_") and k not in ("mean_best_value", "mean_rest_value", "mean_margin"):
            cell_types.append(k[len("mean_"):])
    if best_group not in cell_types:
        cell_types.append(best_group)
    cell_types = sorted(dict.fromkeys(cell_types))  # stable order

    # logger.info(f"[{chrom}:{start}-{end}] best_group={best_group} best_dir={best_dir} "
    #             f"cell_types={cell_types} (n={len(cell_types)})")

    try:
        # load methylation matrix + CpG positions
        mm, cpgs = methyl_matrix_from_bam(
            input_file, reference, chrom=chrom, start=start, end=end, return_positions=True
        )
        n_reads_raw = 0 if mm is None else mm.shape[0]
        n_cpgs = len(cpgs)
        if n_cpgs == 0:
            logger.warning(f"[{chrom}:{start}-{end}] no CpGs found; skipping")
            return None

        # drop rows that are entirely NaN across CpGs
        mm = mm.dropna(how="all")
        if mm.empty or mm.shape[0] < 2:
            logger.warning(f"[{chrom}:{start}-{end}] usable_reads={mm.shape[0] if not mm.empty else 0} "
                           f"(raw={n_reads_raw}) < 2; skipping")
            return None

        # call your untouched kmeans wrapper to assign target vs Other
        dmr_row = {
            "best_group": best_group,
            "best_dir": best_dir,
            "mean_best_value": row.get("mean_best_value", np.nan),
            "mean_rest_value": row.get("mean_rest_value", np.nan),
        }
        out = kmeans_cluster_cells(mm, dmr_row=dmr_row)
        
        # CpG bounds from cpgs
        cpgstart = int(cpgs[0])
        cpgend   = int(cpgs[-1])

        # read names from MultiIndex level 0 if present
        if isinstance(mm.index, pd.MultiIndex) and "read_name" in mm.index.names:
            readnames = mm.index.get_level_values("read_name").astype(str).values
        else:
            readnames = (mm.index.astype(str).values if mm.index.dtype == object
                         else np.array([f"read_{i}" for i in range(len(mm))], dtype=str))

        # target mask from your output column
        mask_target = (out["celltype_or_other"].astype(str).str.lower()
                       == best_group.strip().lower()).values
        # build variable-length code strings
        pos = {ct: i for i, ct in enumerate(cell_types)}
        t_idx = pos[best_group]
        target_bits = ["0"] * len(cell_types); target_bits[t_idx] = "1"
        other_bits  = ["1"] * len(cell_types); other_bits[t_idx]  = "0"
        target_code = "".join(target_bits)
        other_code  = "".join(other_bits)
        code_col = np.where(mask_target, target_code, other_code)

        # --- per-read assignments (each read = one row / index) ---
        assign_df = pd.DataFrame({
            "chr": chrom,
            "start": start,
            "end": end,
            "cpgstart": cpgstart,
            "cpgend": cpgend,
            "code_order": "|".join(cell_types),
            "code": code_col,
        }, index=pd.Index(readnames, name="readname"))

        # per-block means (per-read mean methylation, then avg by target vs other)
        X_imp = mm.astype(float).copy().fillna(mm.astype(float).mean())
        read_mean = X_imp.mean(axis=1).values
        tgt_mean = float(np.nanmean(read_mean[mask_target])) if mask_target.any() else np.nan
        oth_mean = float(np.nanmean(read_mean[~mask_target])) if (~mask_target).any() else np.nan

        # logger.info(f"[{chrom}:{start}-{end}] target_mean={tgt_mean:.4f} other_mean={oth_mean:.4f} "
        #             f"cpg_bounds={cpgstart}-{cpgend}")

        state_payload = {
            "chr": chrom, "start": start, "end": end,
            "cpgstart": cpgstart, "cpgend": cpgend,
        }
        for ct in cell_types:
            state_payload[f"{ct}_methylation"] = tgt_mean if ct == best_group else oth_mean
        state_df = pd.DataFrame([state_payload])

        return assign_df, state_df

    except Exception as e:
        logger.exception(f"[{chrom}:{start}-{end}] failed with error")
        return None

def sv_anno(args):
    logger = logging.getLogger("anno.sv_anno")
    logger.info("Starting SV annotation from pre-annotated reads")
    if args.command == "svanno":    
        input_file = args.input
    else:
        input_file = os.path.join(args.output, "reads_classification.tsv")
    if args.kanpig_read_names is not None:
        logger.info(f"Using kanpig read names from: {args.kanpig_read_names}")
    else:
        logger.info("No kanpig read names provided; using Sniffles read names from VCF")
    sv_assignment_df = assign_sv_celltypes(read_vcf_to_df(args.vcf, kanpig_read_names=args.kanpig_read_names), pd.read_csv(input_file, sep="\t", index_col=0))
    sv_assignment_df.to_csv(os.path.join(args.output, "sv_assignment.tsv"), sep="\t", index=False)




def anno_main(args):
    # print(args)
    # return
    logger = logging.getLogger("anno.main")

    bed_file   = args.bed
    base_out   = args.output    # writes <output>.reads.tsv and <output>.blocks.tsv
    input_file = args.input
    reference  = args.reference
    threads    = int(args.threads)
    window     = int(args.window)
    logger.info(f"Starting annotation: bed={bed_file} bam={input_file} ref={reference} "
                f"threads={threads} out_base={base_out}")

    # Output paths
    reads_out  = os.path.join(base_out, "reads_classification.tsv")
    blocks_out = os.path.join(base_out, "blocks_classification.tsv")

    # Load and (optionally) filter BED
    bed = pd.read_csv(bed_file, sep="\t")
    logger.info(f"Loaded BED with {len(bed)} DMR rows")

    sv_df = read_vcf_to_df(args.vcf)
    filtered_bed = filter_bed_based_on_variants(bed, sv_df=sv_df, window=window)

    for col in ["chr", "start", "end", "best_group", "best_dir"]:
        if col not in filtered_bed.columns:
            logger.error(f"BED missing required column: {col}")
            raise ValueError(f"BED missing required column: {col}")

    n_tasks = len(filtered_bed)
    logger.info(f"Filtered BED to {n_tasks} DMRs after variant overlap filtering, window size = {window}")

    tasks = [(dict(row), input_file, reference) for _, row in filtered_bed.iterrows()]

    # --- Prepare outputs: truncate files and reset header flags ---
    # We'll only write headers on the first real chunk for each file.
    open(reads_out,  "w").close()
    open(blocks_out, "w").close()
    reads_header_written  = False
    blocks_header_written = False
    blocks_cols_locked: list[str] | None = None  # we lock schema to the first block we see

    # Stream results and append immediately
    with mp.Pool(threads) as pool:
        for res in tqdm(pool.imap(_one_dmr, tasks, chunksize=1),
                        total=n_tasks, desc="Processing DMRs"):
            if res is None:
                continue

            a_df, s_df = res

            # --- APPEND READS ---
            if a_df is not None and not a_df.empty:
                if not reads_header_written:
                    # first write: include header and index (readname)
                    a_df.to_csv(reads_out, sep="\t", index=True, mode="a", header=True)
                    reads_header_written = True
                else:
                    a_df.to_csv(reads_out, sep="\t", index=True, mode="a", header=False)

            # --- APPEND BLOCKS (variable columns across DMRs) ---
            if s_df is not None and not s_df.empty:
                if not blocks_header_written:
                    # lock the schema to the first encountered block columns
                    blocks_cols_locked = list(s_df.columns)
                    s_df.to_csv(blocks_out, sep="\t", index=False, mode="a", header=True)
                    blocks_header_written = True
                else:
                    # align columns to locked header; drop extras, add missing as NaN
                    assert blocks_cols_locked is not None
                    s_df_aligned = s_df.reindex(columns=blocks_cols_locked)
                    s_df_aligned.to_csv(blocks_out, sep="\t", index=False, mode="a", header=False)

    # If nothing was written, emit empty files with headers to be friendly downstream
    if not reads_header_written:
        empty_reads = pd.DataFrame(
            columns=["chr","start","end","cpgstart","cpgend","code_order","code"]
        )
        empty_reads.index.name = "readname"
        empty_reads.to_csv(reads_out, sep="\t", index=True, header=True)
        logger.warning("No per-read assignments generated; wrote empty reads header only")

    if not blocks_header_written:
        pd.DataFrame(columns=["chr","start","end","cpgstart","cpgend"]).to_csv(
            blocks_out, sep="\t", index=False, header=True
        )
        logger.warning("No block states generated; wrote empty blocks header only")
    sv_anno(args)
    logger.info("Annotation complete")


