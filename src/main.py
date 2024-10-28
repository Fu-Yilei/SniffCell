import sys, src.parse_args, src.calc_methylation_diff_regions, pysam, os
import pandas as pd
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
    benchmark_second_vcf = args.second_vcf
    min_supporting_reads = args.min_supporting

    os.makedirs(output, exist_ok=True)
    # input_bam = pysam.AlignmentFile(input_bam)
    # input_vcf = pysam.VariantFile(input_vcf)
    # reference_genome = pysam.FastaFile(reference_genome)
    
    sv_methylation_df = src.calc_methylation_diff_regions.calculate_methylation_diff_region_bam(sv_vcf=input_vcf, input_bam=input_bam, 
                                                                            reference_genome=reference_genome, output_bam_folder=output, 
                                                                            output_bam=output_bam, min_supporting_read_num = min_supporting_reads, 
                                                                            sv_discovery_range=sv_discovery_interval)

    sv_methylation_df.to_csv(os.path.join(output, "sv_methylation_df.csv"))
    

