import pysam, sys, json, os
from src.deconv.deconv import get_base_modification_dictionary_basic_supporting_reads, filter_cell_type_regions, assign_read_to_readgroup
from src.deconv.deconv_argparse import parse_args
import pandas as pd
import numpy as np
from src.deconv.deconv_em import em_k_cluster_methylation
import logging
from tqdm import tqdm
from multiprocessing import Pool

def process_individual_region(filtered_regions, celltypes, bam_file_path, output, ref_seq_path, threshold, verbose=False):
    chromosome = filtered_regions.chr
    phase_region = [filtered_regions.start, filtered_regions.end]
    # logging.info(f"Processing region {chromosome}:{phase_region[0]}-{phase_region[1]}")

    bam_file = pysam.AlignmentFile(bam_file_path)
    ref_seq = pysam.FastaFile(ref_seq_path)
    reads_methylation_df = get_base_modification_dictionary_basic_supporting_reads(bam_file, ref_seq, chromosome, phase_region)
    if reads_methylation_df.empty:
        exception = f"No reads found in the region {chromosome}:{phase_region[0]}-{phase_region[1]}"
        logging.error(exception)
        return None

    p_init_region = filtered_regions[celltypes]
    alpha_init = 1 / len(celltypes) * np.ones(len(celltypes))
    alpha_final, p_final, gamma_final = em_k_cluster_methylation(
        df=reads_methylation_df,
        alpha_init=alpha_init,
        p_init=p_init_region,  # same single p for all sites initially
        max_iter=50,
        tol=1e-5,
        random_state=42
    )
    if verbose:
        gamma_final.to_csv(f"{output}/{chromosome}_{phase_region[0]}_{phase_region[1]}.tsv", sep='\t')
    assigned_reads_count, cell_types_list = assign_read_to_readgroup(filtered_regions, celltypes, gamma_final, bam_file, output, gamma_max_confidence=threshold)
    return pd.DataFrame([{"chr": chromosome, "start": phase_region[0], "end": phase_region[1], "total_reads": reads_methylation_df.shape[0], "assigned_reads": assigned_reads_count, "cell_type_reads_counts": cell_types_list, "cell_type_prob_em": alpha_final}])

def process_individual_region_wrapper(args):
    return process_individual_region(*args)


def main(argv):    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = parse_args(argv)
    # print(args)
    bam_file = args.bam
    ref_seq = args.reference
    tissue = args.tissue
    atlas = args.atlas
    output = args.output
    os.makedirs(output, exist_ok=True)
    threads = args.threads
    verbose = args.verbose
    region_number = args.region_number
    threshold = args.threshold
    methods = args.method
    logging.info(f"Parameters - BAM: {bam_file}, Reference: {ref_seq}, Tissue: {tissue}, Atlas: {atlas}, Output: {output}, Threads: {threads}, Verbose: {verbose}, Region number: {region_number}, Threshold: {threshold}, Methods: {methods}")
    # Get the directory of the current script
    script_dir = os.path.dirname(__file__)
    tissue_celltypes_file_path = os.path.join(script_dir, 'atlas', 'tissue_celltypes.json')
    atlas_celltypes_file_path = atlas #os.path.join(script_dir, 'atlas', '39Bisulfite.tsv')

    with open(tissue_celltypes_file_path, 'r') as tissue_celltypes:
        celltype_dict = json.load(tissue_celltypes)

    # region_number = 300
    celltypes = celltype_dict[tissue]["cell_type"]

    filtered_regions_df = filter_cell_type_regions(celltypes, atlas_celltypes_file_path, region_number, by=methods)
    summary_classification_df = pd.DataFrame(columns=["chr", "start", "end", "total_reads", "assigned_reads", "cell_type_reads_counts", "cell_type_prob_em"])
    args_list = [(filtered_regions, celltypes, bam_file, output, ref_seq, threshold, verbose) for _, filtered_regions in filtered_regions_df.iterrows()]

    with Pool(threads) as pool:
        summary_classification_series = list(tqdm(pool.imap(process_individual_region_wrapper, args_list), total=len(args_list), desc="Processing SVs"))
    for summary_classification in summary_classification_series:
        summary_classification_df = pd.concat([summary_classification_df, summary_classification], ignore_index=True)   
    summary_classification_df.to_csv(f"{output}/summary_classification.tsv", sep='\t', index=False)
    logging.info("Processing complete. Summary classification saved.")




