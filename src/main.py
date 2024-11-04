import src.get_base_modification_dictionary
import sys, src.parse_args, src.calc_methylation_diff_regions, pysam, os
import pandas as pd
from tqdm import tqdm
from src import *
def main(argv):
    args = src.parse_args.parse_args(argv)
    input_bam = args.bam
    input_vcf = args.vcf
    reference_genome = args.reference
    output = args.output
    threads = args.threads
    output_bam = args.output_bam
    smoothing = args.smoothing # int
    sv_discovery_interval = args.interval
    region_decide_threshold = args.threshhold
    benchmark_second_bam = args.second_bam
    min_supporting_reads = args.min_supporting
    test_function = args.test_function

    os.makedirs(output, exist_ok=True)
    
    if benchmark_second_bam:
        benchmarking_sv_df = src.calc_methylation_diff_regions.calculate_methylation_diff_region_benchmark(input_bam, 
                                                           input_vcf, reference_genome, 
                                                           output, sv_discovery_interval, 
                                                           benchmark_second_bam, min_supporting_reads, 
                                                           test_function, smoothing=smoothing, 
                                                           threads=threads)
        benchmarking_sv_df.to_csv(os.path.join(output, "benchmarking_sv_df.csv"))

    else:
        sv_methylation_df = src.calc_methylation_diff_regions.calculate_methylation_diff_region_bam(sv_vcf=input_vcf, input_bam=input_bam, 
                                            reference_genome=reference_genome, output_bam_folder=output, smoothing=smoothing, 
                                            output_bam=output_bam, min_supporting_read_num = min_supporting_reads, 
                                            sv_discovery_range=sv_discovery_interval, test_function=test_function, threads=threads)
        sv_methylation_df.to_csv(os.path.join(output, "sv_methylation_df.csv"))


 