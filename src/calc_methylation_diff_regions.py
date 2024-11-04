import pysam, os, re
import src
import src.smoothing
from src.vcf_to_df import read_vcf_to_df
from src.get_base_modification_dictionary import get_base_modification_dictionary_basic_supporting_reads, get_base_modification_dictionary_new_bam
from src.statistic_tests import calculate_fisher_basic, calculate_ranksum_basic, calculate_ttest, calculate_ranksum, calculate_fisher, calculate_ttest_basic
from tqdm import tqdm
from multiprocessing import Pool


def filter_dict_with_sv(input_dict, sv_start, sv_end):
    """
    Filters a dictionary by excluding keys that fall within a specified range.

    Args:
        input_dict (dict): The dictionary to be filtered.
        sv_start (int): The start of the range.
        sv_end (int): The end of the range.

    Returns:
        dict: A new dictionary with keys outside the specified range.
    """
    output_dict = {}
    for i in input_dict:
        if i < sv_start or i > sv_end:
            output_dict.update({i: input_dict[i]})
    return output_dict


def process_individual_sv(sv, input_bam, reference_genome, output_bam_folder, output_bam, sv_discovery_range, smoothing, test_function):
    hapmap_ref_file = pysam.AlignmentFile(input_bam)
    referece_sequence = pysam.Fastafile(reference_genome)
    modification_dict = get_base_modification_dictionary_new_bam(bam_file=hapmap_ref_file, 
                                                             ref_seq=referece_sequence, 
                                                             chromosome=sv.chr, 
                                                             phase_region=(int(sv.ref_start)-sv_discovery_range, int(sv.ref_end)+sv_discovery_range), 
                                                             sv_supporting_reads=sv.supporting_reads, sv_id=sv.id, output_bam_folder=output_bam_folder, output_bam=output_bam)
    modification_dict = src.smoothing.smoothing_comp(modification_dict, key_diff_threshold=smoothing)
    if test_function == 'ttest':
        ranksum_dict = calculate_ttest(modification_dict)
    elif test_function == 'ranksum':
        ranksum_dict = calculate_ranksum(modification_dict)
    elif test_function == 'fisher':
        ranksum_dict = calculate_fisher(modification_dict)
    else:
        ranksum_dict = calculate_ttest(modification_dict)

    # ranksum_dict = calculate_ttest(modification_dict) # nan if list length is 0
    ranksum_dict = filter_dict_with_sv(ranksum_dict, int(sv.ref_start), int(sv.ref_end))
    cpg_diff_locs = {i: p_value for i, p_value in ranksum_dict.items() if p_value < 0.05}
    cpg_same_locs = {i: p_value for i, p_value in ranksum_dict.items() if p_value > 0.05}
    # Calculate the percentage of different CpG locations
    if len(cpg_diff_locs) + len(cpg_same_locs) > 0:
        cpg_diff_percentage = len(cpg_diff_locs) / (len(cpg_diff_locs) + len(cpg_same_locs))
    else:
        cpg_diff_percentage = None  # or some default value

    return cpg_diff_percentage

def process_sv_wrapper(args):
    return process_individual_sv(*args)

def calculate_methylation_diff_region_bam(sv_vcf, input_bam, reference_genome, 
                                          output_bam_folder, smoothing, threads, output_bam=True, 
                                          min_supporting_read_num=5, 
                                          sv_discovery_range=1000,
                                          test_function='ttest'):
    
    filtered_mosaic_sv = pysam.VariantFile(sv_vcf)
    filtered_mosaic_sv_df = read_vcf_to_df(filtered_mosaic_sv)
    filtered_mosaic_sv_df_filtered = filtered_mosaic_sv_df[filtered_mosaic_sv_df.supporting_reads.apply(len) > min_supporting_read_num]

    args_list = [(test_sv, input_bam, reference_genome, output_bam_folder, output_bam, sv_discovery_range, smoothing, test_function) 
                 for _, test_sv in filtered_mosaic_sv_df_filtered.iterrows()]

    with Pool(threads) as pool:
        cpg_diff_percentages = list(tqdm(pool.imap(process_sv_wrapper, args_list), total=len(args_list), desc="Processing SVs"))

    filtered_mosaic_sv_df_filtered['cpg_diff_percentage'] = cpg_diff_percentages

    return filtered_mosaic_sv_df_filtered


def process_individual_sv_benchmark(sv, input_bam, reference_genome, output_bam_folder, output_bam, sv_discovery_range, benchmark_second_bam, smoothing, test_function):
    hapmap_ref_file = pysam.AlignmentFile(input_bam)
    benchmark_second_bam = pysam.AlignmentFile(benchmark_second_bam)
    referece_sequence = pysam.Fastafile(reference_genome)
    modification_dict = get_base_modification_dictionary_basic_supporting_reads(bam_file=hapmap_ref_file, 
                                                             ref_seq=referece_sequence, 
                                                             chromosome=sv.chr, 
                                                             phase_region=(int(sv.ref_start)-sv_discovery_range, int(sv.ref_end)+sv_discovery_range), 
                                                             sv_supporting_reads=sv.supporting_reads, sv_id=sv.id, output_bam_folder=output_bam_folder, output_bam=output_bam)
    modification_dict = src.smoothing.smoothing_comp(modification_dict, key_diff_threshold=smoothing)
    benchmark_second_bam_sv_methylation_dict = get_base_modification_dictionary_basic_supporting_reads(bam_file=benchmark_second_bam, 
                                                             ref_seq=referece_sequence, 
                                                             chromosome=sv.chr, 
                                                             phase_region=(int(sv.ref_start)-sv_discovery_range, int(sv.ref_end)+sv_discovery_range), 
                                                             sv_supporting_reads=None, sv_id=sv.id, output_bam_folder=output_bam_folder, output_bam=False)
    benchmark_second_bam_sv_methylation_dict = src.smoothing.smoothing_comp(benchmark_second_bam_sv_methylation_dict, key_diff_threshold=smoothing)
    if test_function == 'ttest':
        ranksum_dict = calculate_ttest_basic(modification_dict, benchmark_second_bam_sv_methylation_dict)
    elif test_function == 'ranksum':
        ranksum_dict = calculate_ranksum_basic(modification_dict, benchmark_second_bam_sv_methylation_dict)
    elif test_function == 'fisher':
        ranksum_dict = calculate_fisher_basic(modification_dict, benchmark_second_bam_sv_methylation_dict)
    else:
        ranksum_dict = calculate_ttest_basic(modification_dict, benchmark_second_bam_sv_methylation_dict)

    # ranksum_dict = calculate_ttest(modification_dict) # nan if list length is 0
    ranksum_dict = filter_dict_with_sv(ranksum_dict, int(sv.ref_start), int(sv.ref_end))
    cpg_diff_locs = {i: p_value for i, p_value in ranksum_dict.items() if p_value < 0.05}
    cpg_same_locs = {i: p_value for i, p_value in ranksum_dict.items() if p_value > 0.05}
    # Calculate the percentage of different CpG locations
    if len(cpg_diff_locs) + len(cpg_same_locs) > 0:
        cpg_diff_percentage = len(cpg_diff_locs) / (len(cpg_diff_locs) + len(cpg_same_locs))
    else:
        cpg_diff_percentage = None  
    return cpg_diff_percentage

def process_sv_wrapper_benchmark(args):
    return process_individual_sv_benchmark(*args)


def calculate_methylation_diff_region_benchmark(input_bam, input_vcf, reference_genome, output, sv_discovery_interval, benchmark_second_bam, min_supporting_reads, test_function, smoothing, threads, output_bam):
    input_vcf = pysam.VariantFile(input_vcf)
    benchmarking_vcf = read_vcf_to_df(input_vcf)
    benchmarking_vcf = benchmarking_vcf[benchmarking_vcf.supporting_reads.apply(len) > min_supporting_reads]
    args_list = [(test_sv, input_bam, reference_genome, output, output_bam, sv_discovery_interval, benchmark_second_bam, smoothing, test_function) 
                 for _, test_sv in benchmarking_vcf.iterrows()]

    with Pool(threads) as pool:
        cpg_diff_percentages = list(tqdm(pool.imap(process_sv_wrapper_benchmark, args_list), total=len(args_list), desc="Processing SVs"))
    benchmarking_vcf['cpg_diff_percentage'] = cpg_diff_percentages
    return benchmarking_vcf

