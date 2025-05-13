import src.get_base_modification_dictionary
import sys
import src.parse_args
import src.calc_methylation_diff_regions
import pysam
import os
import pandas as pd
from tqdm import tqdm
# from src import main, parse_args, vcf_to_df
import logging
import json
import numpy as np
from multiprocessing import Pool
import multiprocessing
from src.deconv import get_base_modification_dictionary_basic_supporting_reads, filter_cell_type_regions, assign_read_to_readgroup
from src.deconv_em import em_k_cluster_methylation
from src.deconv_estimate import estimate_celltype_assignment
os.environ["HTS_LOG_LEVEL"] = "error"

def process_individual_region(filtered_regions, celltypes, bam_file_path, output, ref_seq_path, threshold, celltype_prior, verbose=False):
    process = multiprocessing.current_process()
    process.name = f"SniffMeth"
    chromosome = filtered_regions.chr
    phase_region = [filtered_regions.start, filtered_regions.end]

    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    ref_seq = pysam.FastaFile(ref_seq_path)
    reads_methylation_df = get_base_modification_dictionary_basic_supporting_reads(bam_file, ref_seq, chromosome, phase_region)
    if reads_methylation_df.empty:
        exception = f"No reads found in the region {chromosome}:{phase_region[0]}-{phase_region[1]}"
        logging.error(exception)
        return None

    p_init_region = filtered_regions[celltypes]
    if celltype_prior:
        alpha_init = celltype_prior
    else:
        alpha_init = 1 / len(celltypes) * np.ones(len(celltypes))
    alpha_final, _, gamma_final = em_k_cluster_methylation(
        df=reads_methylation_df, alpha_init=alpha_init, p_init=p_init_region,
        max_iter=100, tol=1e-5, random_state=42
    )
    # print(alpha_final)
    if verbose:
        gamma_final.to_csv(f"{output}/{chromosome}_{phase_region[0]}_{phase_region[1]}.tsv", sep='\t')
    assigned_reads_count, cell_types_list = assign_read_to_readgroup(
        filtered_regions, celltypes, gamma_final, bam_file, output, gamma_max_confidence=threshold
    )
    return pd.DataFrame([{
        "chr": chromosome, "start": phase_region[0], "end": phase_region[1],
        "total_reads": reads_methylation_df.shape[0], "assigned_reads": assigned_reads_count,
        "cell_type_reads_counts": cell_types_list, "cell_type_prob_em": alpha_final
    }])

def process_individual_region_wrapper(args):
    return process_individual_region(*args)

def main(argv):
    # Parse command line arguments
    args = src.parse_args.parse_args(argv)
    input_bam = args.bam
    input_vcf = args.vcf
    reference_genome = args.reference
    output = args.output
    threads = args.threads
    output_bam = args.output_bam
    smoothing = args.smoothing  # int
    sv_discovery_interval = args.interval
    region_decide_threshold = args.threshold
    benchmark_second_bam = args.second_bam
    min_supporting_reads = args.min_supporting
    test_function = args.test_function
    primary_only = args.primary_only

    # Create output directory if it doesn't exist
    os.makedirs(output, exist_ok=True)
    # sniffmeth_output_bam_folder = os.path.join(output, "methylation_sv_bam")
    # Benchmarking mode
    if benchmark_second_bam:
        benchmarking_sv_df = src.calc_methylation_diff_regions.calculate_methylation_diff_region_benchmark(
            input_bam, input_vcf, reference_genome, output, sv_discovery_interval=sv_discovery_interval,
            benchmark_second_bam=benchmark_second_bam, min_supporting_reads=min_supporting_reads,
            test_function=test_function, smoothing=smoothing, threads=threads, output_bam=output_bam
        )
        benchmarking_sv_df.to_csv(os.path.join(output, "benchmarking_sv_df.csv"))
    else:
        # Normal mode
        sv_methylation_df = src.calc_methylation_diff_regions.calculate_methylation_diff_region_bam(
            sv_vcf=input_vcf, input_bam=input_bam, reference_genome=reference_genome, output_bam_folder=output,
            smoothing=smoothing, output_bam=output_bam, min_supporting_read_num=min_supporting_reads,
            sv_discovery_range=sv_discovery_interval, test_function=test_function, threads=threads
        )
        sv_methylation_df.to_csv(os.path.join(output, "sv_methylation_df.csv"))

    # Deconvolution mode
    if args.deconv:
        deconv_atlas = args.atlas
        deconv_tissue = args.tissue
        deconv_region_number = args.region_number
        deconv_method = args.method
        deconv_confidence = args.confidence
        deconv_verbose = args.verbose

        deconv_output_location = os.path.join(output, "deconv_bam_output")
        os.makedirs(deconv_output_location, exist_ok=True)
        input_vcf_filename = os.path.basename(input_vcf).split(".")[0]


        deconv_folder = os.path.dirname(os.path.abspath(deconv_atlas))
        with open(os.path.join(deconv_folder, "tissue_celltypes.json"), 'r', encoding='utf-8') as tissue_celltypes:
            celltype_dict = json.load(tissue_celltypes)
        celltypes = celltype_dict[deconv_tissue]["cell_type"]
        celltypes_prior = celltype_dict[deconv_tissue]["prior"]
        
        logging.info("Parameters - BAM: %s, Reference: %s, Tissue: %s, Atlas: %s, Output: %s, Threads: %s, Verbose: %s, Region number: %s, Threshold: %s, Methods: %s, Using prior: %s", 
                     input_bam, reference_genome, deconv_tissue, deconv_atlas, deconv_output_location, threads, deconv_verbose, deconv_region_number, deconv_confidence, deconv_method, celltypes_prior)
        filtered_regions_df = filter_cell_type_regions(celltypes, deconv_atlas, deconv_region_number, by=deconv_method)
        summary_classification_df = pd.DataFrame(columns=["chr", "start", "end", "total_reads", "assigned_reads", "cell_type_reads_counts", "cell_type_prob_em"])
        args_list = [(filtered_regions, celltypes, input_bam, deconv_output_location, reference_genome, deconv_confidence, celltypes_prior, deconv_verbose) for _, filtered_regions in filtered_regions_df.iterrows()]
        with Pool(threads) as p:
            summary_classification_series = list(tqdm(p.imap(process_individual_region_wrapper, args_list), total=len(args_list), desc="Processing SVs"))
        for summary_classification in summary_classification_series:
            summary_classification_df = pd.concat([summary_classification_df, summary_classification], ignore_index=True)
        summary_classification_df.to_csv(f"{output}/summary_classification.csv", sep='\t', index=False)
        logging.info("Processing complete. Summary classification saved.")
        deconv_output_vcf = os.path.join(output, f"{input_vcf_filename}.celltype_annotated.vcf")
        
        estimate_celltype_assignment(input_vcf, sv_methylation_df, summary_classification_df, celltypes, deconv_output_vcf)
        logging.info("VCF file annotated with cell type assignment.")