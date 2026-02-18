import os
import json
import pandas as pd
from sniffcell.dmsv.sv_methylation import get_methylation_around_sv
from sniffcell.dmsv.statistical_test_around_sv import get_statistical_tests_around_sv
from sniffcell.anno.vcf_to_df import read_vcf_to_df
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(processName)s] %(levelname)s: %(message)s")

def dmsv_main(args):
    logger = logging.getLogger("dmsv.main")
    input_bam = args.input
    input_vcf = args.vcf
    min_supporting = args.min_supporting
    output_folder = args.output
    reference_fasta = args.reference
    flank_size = args.flank_size
    test_type = args.test_type
    haplotype_majority_threshold = args.haplotype_majority_threshold
    threads = args.threads
    min_cpgs = args.min_cpgs
    os.makedirs(output_folder, exist_ok=True)
    logger.info(f"Starting DMSV analysis: min_supporting={min_supporting} output={output_folder}")  

    sv_df = read_vcf_to_df(input_vcf)  # limit to first 10 SVs for testing
    sv_methylation_data = get_methylation_around_sv(
        input_bam, reference_fasta, sv_df,
        min_supporting_reads=min_supporting,
        flank_size=flank_size, n_threads=threads
    )
    # print(sv_methylation_data)
    sv_statistical_tests = get_statistical_tests_around_sv(
        sv_methylation_data,
        haplotype_majority_threshold=haplotype_majority_threshold
    )

    # --- Output section (safe against empty slices) ---
    logger.info("Writing results to output folder...")

    alpha = 0.05
    min_hits = int(min_cpgs) if isinstance(min_cpgs, (int, float)) else 5

    details_dir = os.path.join(output_folder, "sv_details")
    os.makedirs(details_dir, exist_ok=True)

    def _safe_median(x: pd.Series) -> float:
        if x is None:
            return float("nan")
        x = pd.Series(x).dropna()
        return float(x.median()) if len(x) else float("nan")

    def _safe_max_abs(x: pd.Series) -> float:
        if x is None:
            return float("nan")
        x = pd.Series(x).dropna()
        return float(x.abs().max()) if len(x) else float("nan")

    rows = []
    for sv_id, stats in sv_statistical_tests.items():
        # Case A: modern return -> DataFrame
        if isinstance(stats, pd.DataFrame):
            df_stats = stats.copy()

            # If this SV ended up with zero CpGs, write an empty detail file and continue
            total_cpgs = int(df_stats.shape[0])
            out_detail = os.path.join(details_dir, f"{sv_id}.tsv.gz")
            # Write detail file (even if empty) to make outputs predictable
            df_stats.to_csv(out_detail, sep="\t", compression="infer")

            if total_cpgs == 0:
                rows.append({
                    "sv_id": sv_id,
                    "num_significant_CpGs": 0,
                    "total_CpGs": 0,
                    "median_abs_delta_sig": float("nan"),
                    "max_abs_delta_sig": float("nan"),
                    "n_hyper_sig": 0,
                    "n_hypo_sig": 0,
                    "median_agreement_sup": float("nan"),
                })
                continue

            # significance mask
            if "mwu_q" in df_stats.columns:
                sig_mask = (df_stats["mwu_q"] <= alpha)
            elif "mwu_p" in df_stats.columns:
                sig_mask = (df_stats["mwu_p"] <= alpha)
            else:
                sig_mask = pd.Series(False, index=df_stats.index)

            # require 'keep' if available
            if "keep" in df_stats.columns:
                sig_mask = sig_mask & df_stats["keep"].fillna(False)

            num_sig = int(sig_mask.sum())

            # effect summaries on sig subset (guard for empties)
            if "delta_mean" in df_stats.columns and num_sig > 0:
                deltas = df_stats.loc[sig_mask, "delta_mean"].dropna()
                median_abs_delta = _safe_median(deltas.abs())
                max_abs_delta    = _safe_max_abs(deltas)
                n_hyper = int((deltas > 0).sum())
                n_hypo  = int((deltas < 0).sum())
            else:
                median_abs_delta = float("nan")
                max_abs_delta    = float("nan")
                n_hyper = 0
                n_hypo  = 0

            # agreement summary (guard for empties)
            if "agreement_sup" in df_stats.columns:
                agreement_median = _safe_median(df_stats["agreement_sup"])
            else:
                agreement_median = float("nan")

            rows.append({
                "sv_id": sv_id,
                "num_significant_CpGs": num_sig,
                "total_CpGs": total_cpgs,
                "median_abs_delta_sig": median_abs_delta,
                "max_abs_delta_sig": max_abs_delta,
                "n_hyper_sig": n_hyper,
                "n_hypo_sig": n_hypo,
                "median_agreement_sup": agreement_median,
            })

        # Case B: legacy return -> dict {cpg: p}
        else:
            try:
                s = pd.Series(stats, dtype="float64")
            except Exception:
                s = pd.Series(dtype="float64")

            s_nonan = s.dropna()
            total_cpgs = int(s_nonan.size)
            num_sig = int((s_nonan <= alpha).sum())

            rows.append({
                "sv_id": sv_id,
                "num_significant_CpGs": num_sig,
                "total_CpGs": total_cpgs,
                "median_abs_delta_sig": float("nan"),
                "max_abs_delta_sig": float("nan"),
                "n_hyper_sig": 0,
                "n_hypo_sig": 0,
                "median_agreement_sup": float("nan"),
            })

    # Build summary
    summary_df = pd.DataFrame(rows).set_index("sv_id")

    # Join meta from sv_df (guard missing columns)
    meta_cols = [c for c in ["chr", "location", "vaf", "sv_len"] if c in sv_df.columns]
    if "id" in sv_df.columns and meta_cols:
        summary_df = summary_df.join(sv_df.set_index("id")[meta_cols], how="left")

    # Flag significant SVs using CLI min_cpgs
    summary_df["is_significant"] = summary_df["num_significant_CpGs"] >= min_hits

    # Sort robustly even if some columns are NaN/empty
    sort_cols = [c for c in ["num_significant_CpGs", "max_abs_delta_sig"] if c in summary_df.columns]
    if sort_cols:
        summary_df = summary_df.sort_values(sort_cols, ascending=[False]*len(sort_cols))

    # Write summary
    out_summary = os.path.join(output_folder, "significant_SVs.tsv")
    summary_df.to_csv(out_summary, sep="\t")
    logger.info(f"Summary written to {out_summary} (details in {details_dir}/)")
