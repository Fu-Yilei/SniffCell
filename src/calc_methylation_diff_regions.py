import pysam, os, re
from src.vcf_to_df import read_vcf_to_df
from src.get_base_modification_dictionary import get_base_modification_dictionary_new_bam
from src.statistic_tests import calculate_ttest, calculate_ranksum, calculate_fisher
from tqdm import tqdm


def filter_dict_with_sv(input_dict, sv_start, sv_end):
    output_dict = {}
    for i in input_dict:
        if i < sv_start or i > sv_end:
            output_dict.update({i: input_dict[i]})
    return output_dict


def calculate_methylation_diff_region_bam(sv_vcf, input_bam, reference_genome, 
                                          output_bam_folder, output_bam=True, 
                                          min_supporting_read_num = 5, 
                                          sv_discovery_range=1000,
                                            test_function = 'ttest'):
    
    hapmap_ref_file = pysam.AlignmentFile(input_bam)
    referece_sequence = pysam.Fastafile(reference_genome)
    filtered_mosaic_sv = pysam.VariantFile(sv_vcf)
    filtered_mosaic_sv_df = read_vcf_to_df(filtered_mosaic_sv)
    cpg_diff_percentages = []
    filtered_mosaic_sv_df_filtered = filtered_mosaic_sv_df[filtered_mosaic_sv_df.supporting_reads.apply(len)>min_supporting_read_num]
    for _, test_sv in tqdm(filtered_mosaic_sv_df_filtered.iterrows(), total=len(filtered_mosaic_sv_df_filtered), desc="Processing SVs"):
        modification_dict = get_base_modification_dictionary_new_bam(bam_file=hapmap_ref_file, 
                                                             ref_seq=referece_sequence, 
                                                             chromosome=test_sv.chr, 
                                                             phase_region=(int(test_sv.ref_start)-sv_discovery_range, int(test_sv.ref_end)+sv_discovery_range), 
                                                             sv_supporting_reads=test_sv.supporting_reads, sv_id=test_sv.id, output_bam_folder=output_bam_folder, output_bam=output_bam)
        if test_function == 'ttest':
            ranksum_dict = calculate_ttest(modification_dict)
        elif test_function == 'ranksum':
            ranksum_dict = calculate_ranksum(modification_dict)
        elif test_function == 'fisher':
            ranksum_dict = calculate_fisher(modification_dict)
        else:
            ranksum_dict = calculate_ttest(modification_dict)

        # ranksum_dict = calculate_ttest(modification_dict) # nan if list length is 0
        ranksum_dict = filter_dict_with_sv(ranksum_dict, int(test_sv.ref_start), int(test_sv.ref_end))
        cpg_diff_locs = {i: p_value for i, p_value in ranksum_dict.items() if p_value < 0.05}
        cpg_same_locs = {i: p_value for i, p_value in ranksum_dict.items() if p_value > 0.05}
        # Calculate the percentage of different CpG locations
        if len(cpg_diff_locs) + len(cpg_same_locs) > 0:
            cpg_diff_percentage = len(cpg_diff_locs) / (len(cpg_diff_locs) + len(cpg_same_locs))
        else:
            cpg_diff_percentage = None  # or some default value
        
        cpg_diff_percentages.append(cpg_diff_percentage)    
    filtered_mosaic_sv_df_filtered['cpg_diff_percentage'] = cpg_diff_percentages

    return filtered_mosaic_sv_df_filtered
