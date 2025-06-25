import os
import sys
import json
import logging
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
from multiprocessing import Pool, current_process
from src import parse_args, calc_methylation_diff_regions
from src.deconv import (
    get_base_modification_dictionary_basic_supporting_reads,
    filter_cell_type_regions,
    assign_read_to_readgroup,
)
from src.deconv_em import em_haplotype_and_combined
from src.deconv_estimate import estimate_celltype_assignment
from src.figures import assign_variant_with_cell_type_names, plot_cell_type_box_distributions
from src.__init__ import __version__ as version

os.environ["HTS_LOG_LEVEL"] = "error"


def process_individual_region(
    filtered_regions, celltypes, bam_file_path, output, ref_seq_path, threshold, celltype_prior, verbose=False
):
    current_process().name = "sniffcell"
    chromosome, phase_region = filtered_regions.chr, [filtered_regions.start, filtered_regions.end]

    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    ref_seq = pysam.FastaFile(ref_seq_path)
    reads_methylation_df = get_base_modification_dictionary_basic_supporting_reads(
        bam_file, ref_seq, chromosome, phase_region
    )

    if reads_methylation_df.empty:
        logging.error(f"No reads found in the region {chromosome}:{phase_region[0]}-{phase_region[1]}")
        return None

    p_init_region = filtered_regions[celltypes]
    alpha_init = celltype_prior if celltype_prior else np.ones(len(celltypes)) / len(celltypes)

    alphas, _, gamma_df = em_haplotype_and_combined(
        df=reads_methylation_df, alpha_init=alpha_init, p_init=p_init_region, max_iter=100, tol=1e-5, random_state=42
    )

    max_distance, alpha_final, gamma_final = None, None, None
    if alphas[0] is None and alphas[1] is None:
        alpha_final = alphas[2]
        gamma_final = gamma_df.copy()
        gamma_final["gamma"] = gamma_df["gamma_list"].apply(
            lambda x: x[-1] if isinstance(x, list) and len(x) > 0 else None
        )
    elif len(alphas) == 3:
        alpha_final = np.mean([i for i in alphas if i is not None], axis=0)
        gamma_final = gamma_df.copy()
        gamma_final["gamma"] = gamma_df["gamma_list"].apply(
            lambda x: np.mean([i for i in x if i is not None], axis=0)
        )
        if alphas[0] is not None and alphas[1] is not None:
            max_distance = max(
                np.linalg.norm(alphas[0] - alphas[1]),
                np.linalg.norm(alphas[0] - alphas[2]),
                np.linalg.norm(alphas[1] - alphas[2]),
            )

    if verbose:
        gamma_final.to_csv(f"{output}/{chromosome}_{phase_region[0]}_{phase_region[1]}.tsv", sep="\t")

    assigned_reads_count, cell_types_list = assign_read_to_readgroup(
        filtered_regions, celltypes, gamma_final, bam_file, output, gamma_max_confidence=threshold
    )

    return pd.DataFrame(
        [
            {
                "chr": chromosome,
                "start": phase_region[0],
                "end": phase_region[1],
                "total_reads": reads_methylation_df.shape[0],
                "assigned_reads": assigned_reads_count,
                "cell_type_reads_counts": cell_types_list,
                "cell_type_prob_em": alpha_final,
                "max_distance": max_distance,
            }
        ]
    )


def process_individual_region_wrapper(args):
    return process_individual_region(*args)


def main(argv):
    args = parse_args.parse_args(argv)
    input_bam, input_vcf, reference_genome, output = args.bam, args.vcf, args.reference, args.output
    threads, output_bam, smoothing = args.threads, args.output_bam or False, args.smoothing
    sv_discovery_interval, region_decide_threshold = args.interval, args.threshold
    benchmark_second_bam, min_supporting_reads = args.second_bam, args.min_supporting
    test_function, primary_only, outlier_thresold = args.test_function, args.primary_only, args.outlier_thresold

    os.makedirs(output, exist_ok=True)

    if benchmark_second_bam:
        benchmark_file = os.path.join(output, "benchmarking_sv_df.csv")
        if os.path.exists(benchmark_file):
            print("Benchmarking results already exist. Skipping benchmarking.")
        else:
            benchmarking_sv_df = calc_methylation_diff_regions.calculate_methylation_diff_region_benchmark(
                input_bam, input_vcf, reference_genome, output, sv_discovery_interval, benchmark_second_bam,
                min_supporting_reads, test_function, smoothing, threads, output_bam
            )
            benchmarking_sv_df.to_csv(benchmark_file)
    else:
        sv_file = os.path.join(output, "sv_methylation_df.csv")
        if os.path.exists(sv_file):
            logging.info("SV methylation data already exists. Skipping calculation.")
            sv_methylation_df = pd.read_csv(sv_file)
        else:
            sv_methylation_df = calc_methylation_diff_regions.calculate_methylation_diff_region_bam(
                sv_vcf=input_vcf, input_bam=input_bam, reference_genome=reference_genome, output_bam_folder=output,
                smoothing=smoothing, output_bam=output_bam, min_supporting_read_num=min_supporting_reads,
                sv_discovery_range=sv_discovery_interval, test_function=test_function, threads=threads
            )
            sv_methylation_df.to_csv(sv_file)

    deconv_output_location = os.path.join(output, "deconv_bam_output")
    os.makedirs(deconv_output_location, exist_ok=True)

    deconv_atlas, deconv_tissue = os.path.abspath(args.atlas), args.tissue
    deconv_region_number, deconv_method = args.region_number, args.method
    deconv_confidence, deconv_verbose = args.confidence, args.verbose
    min_hp_distance, read_fraction_threshold = args.hp_distance, args.assigned_read_fraction

    with open(os.path.join(os.path.dirname(deconv_atlas), "tissue_celltypes.json"), "r", encoding="utf-8") as f:
        celltype_dict = json.load(f)
    celltypes = celltype_dict[deconv_tissue]["cell_type"]

    celltypes_prior = None
    if args.wgbs_tools_uxm:
        from src.wgbs_uxm import run_wgbstools, run_uxm, get_uxm_prior

        run_wgbstools(args.wgbs_path, input_bam, output, args.uxm_atlas, threads)
        uxm_output = run_uxm(args.uxm_path, output, threads, args.uxm_atlas)
        with open(os.path.join(os.path.dirname(deconv_atlas), "column_mapping.json"), "r", encoding="utf-8") as f:
            column_mapping = json.load(f)
        celltypes_prior = get_uxm_prior(uxm_output, column_mapping, selected_cell_types=celltypes)
    else:
        celltypes_prior = celltype_dict[deconv_tissue].get("prior", None)

    filtered_regions_df = filter_cell_type_regions(celltypes, deconv_atlas, deconv_region_number, by=deconv_method)
    summary_classification_df = pd.DataFrame(
        columns=["chr", "start", "end", "total_reads", "assigned_reads", "cell_type_reads_counts", "cell_type_prob_em"]
    )

    args_list = [
        (
            filtered_regions, celltypes, input_bam, deconv_output_location, reference_genome,
            deconv_confidence, celltypes_prior, deconv_verbose
        )
        for _, filtered_regions in filtered_regions_df.iterrows()
    ]

    with Pool(threads) as p:
        summary_classification_series = list(
            tqdm(p.imap(process_individual_region_wrapper, args_list), total=len(args_list),
                    desc="Processing methylation informative regions", unit="region")
        )

    for summary_classification in summary_classification_series:
        summary_classification_df = pd.concat([summary_classification_df, summary_classification], ignore_index=True)

    summary_classification_df = summary_classification_df[
        (summary_classification_df.max_distance <= min_hp_distance) | (summary_classification_df.max_distance.isna())
    ]
    summary_classification_df = summary_classification_df[
        (summary_classification_df.assigned_reads >= read_fraction_threshold * summary_classification_df.total_reads)
    ]

    probs_array = np.stack(summary_classification_df["cell_type_prob_em"].values)
    nan_mask = np.isnan(probs_array).any(axis=1)
    summary_classification_df_clean = summary_classification_df[~nan_mask].copy()
    probs_array_clean = probs_array[~nan_mask]
    mean_profile_clean = probs_array_clean.mean(axis=0)
    summary_classification_df_clean["deviation_from_mean"] = np.linalg.norm(
        probs_array_clean - mean_profile_clean, axis=1
    )
    threshold = summary_classification_df_clean["deviation_from_mean"].quantile(outlier_thresold)
    summary_classification_df_filtered = summary_classification_df_clean[
        summary_classification_df_clean["deviation_from_mean"] <= threshold
    ].copy()

    summary_classification_df = summary_classification_df_filtered
    summary_classification_df.to_csv(f"{output}/summary_classification.csv", sep="\t", index=False)

    cell_type_dicts = assign_variant_with_cell_type_names(summary_classification_df, celltypes)
    plot_cell_type_box_distributions(cell_type_dicts, f"{output}/celltype_prediction_distributions.png")
    logging.info(f"Cell type proportion estimation figure saved to {output}/celltype_prediction_distributions.png")

    sniffcell_version = version
    sniffcell_command = " ".join(sys.argv)

    estimate_celltype_assignment(
        input_vcf, sv_methylation_df, summary_classification_df, celltypes,
        os.path.join(output, f"{os.path.basename(input_vcf).split('.')[0]}.celltype_annotated.vcf"), cmd_info=(sniffcell_version, sniffcell_command),
        assignment_method=deconv_method
    )
    logging.info("VCF file annotated with cell type assignment.")
