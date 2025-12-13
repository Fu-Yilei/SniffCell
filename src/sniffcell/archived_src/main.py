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
from src.deconv.deconv import (
    get_base_modification_dictionary_basic_supporting_reads,
    filter_cell_type_regions,
    assign_read_to_readgroup,
    deconvolution
)
from src.deconv_em import em_haplotype_and_combined
from src.deconv_estimate import estimate_celltype_assignment
from src.figures import assign_variant_with_cell_type_names, plot_cell_type_box_distributions
from src.process_ctdmr_regions import (
    process_individual_region_std,
    process_individual_region_diff,
)
from src.region_deconv import (assign_sv_to_single_celltype_per_region)
from src.__init__ import __version__ as version

os.environ["HTS_LOG_LEVEL"] = "error"


def process_individual_region_std_wrapper(args):
    return process_individual_region_std(*args)
def process_individual_region_diff_wrapper(args):
    return process_individual_region_diff(*args)

def main(argv):
    args = parse_args.parse_args(argv)
    input_bam, input_vcf, reference_genome, output = args.bam, args.vcf, args.reference, args.output
    threads, output_bam, smoothing = args.threads, args.output_bam or False, args.smoothing
    sv_discovery_interval, region_decide_threshold = args.interval, args.threshold
    min_supporting_reads = args.min_supporting
    test_function, primary_only, outlier_thresold = args.test_function, args.primary_only, args.outlier_threshold
    sniffcell_version = version
    sniffcell_command = " ".join(sys.argv)
    os.makedirs(output, exist_ok=True)

    cell_type_dicts, (min_celltype_proportion, max_celltype_proportion) = deconvolution(output, input_bam, reference_genome, args)
    # Variant assignment module
    if args.snp:
        logging.info("SNP methylation analysis is not implemented in this version.")
        return
    else:
        sv_file = os.path.join(output, "sv_methylation_df.csv")
        sv_intermediate_file = os.path.join(output, "sv_methylation_df_intermediate.csv")   
        if os.path.exists(sv_file):
            logging.info("SV methylation data already exists. Skipping calculation.")
            sv_methylation_df = pd.read_csv(sv_file)
        else:
            logging.info("Cell type proportion range for SV assignment: %d to %d", min_celltype_proportion, max_celltype_proportion)
            sv_methylation_df = assign_sv_to_single_celltype_per_region(sv_vcf_path=input_vcf, input_bam_path=input_bam, reference_genome_path=reference_genome, output_bam_folder=output, celltype_proportion_range=(min_celltype_proportion, max_celltype_proportion),
                                                     output_bam=sv_file, min_supporting_read_num=3,
                                                     min_read_supporting_fraction=0.3, cell_type_deconv_range=5000,
                                                     threads=threads, output_file_path=sv_intermediate_file)
            
            sv_methylation_df.to_csv(sv_file)
