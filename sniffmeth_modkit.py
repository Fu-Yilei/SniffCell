#!/usr/bin/env python

import pysam
from src.vcf_to_df import read_vcf_to_df
import logging
from tqdm import tqdm
import os
import argparse

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def sniffmeth_modkit(bamfile, vcf_file, output_bamfile):
    logging.info("Reading VCF file into DataFrame")
    vcf_file = pysam.VariantFile(vcf_file)
    vcf_df = read_vcf_to_df(vcf_file)
    
    logging.info("Opening BAM file")
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    logging.info("Extracting supporting reads from VCF DataFrame")
    all_sv_supporting_reads = set()
    for reads in vcf_df["supporting_reads"]:
        all_sv_supporting_reads.update(reads)
    logging.info(f"Total supporting reads: {len(all_sv_supporting_reads)}")
    
    logging.info("Processing BAM file")
    with pysam.AlignmentFile(output_bamfile, "wb", header=bamfile.header) as out_bamfile:
        for read in tqdm(bamfile, desc="Processing reads"):
            if read.query_name in all_sv_supporting_reads:
                read.set_tag("RG", "SSR")
                out_bamfile.write(read)
    
    logging.info("Closing BAM file")
    bamfile.close()
    logging.info("Processing complete")
def run_modkit(bamfile, output_dir, ref):
    # Construct the command
    command = [
        "modkit", "pileup",
        bamfile,
        output_dir,
        "--ref", ref,
        "--partition-tag", "RG",
        "--prefix", os.path.splitext(os.path.basename(bamfile))[0],
        "--cpg", 
        "--threads", "10"
    ]
    # Run the command
    os.system(" ".join(command))
    


def generate_bedfile(vcf_file, output_bedfile, threshold):
    vcf_file = pysam.VariantFile(vcf_file)
    vcf_df = read_vcf_to_df(vcf_file)
    with open(output_bedfile, "w") as f:
        for _, row in vcf_df.iterrows():
            start = max(0, row['ref_start'] - threshold)
            end = row['ref_end'] + threshold
            f.write(f"{row['chr']}\t{start}\t{end}\n")


def main():
    parser = argparse.ArgumentParser(description="Process BAM and VCF files to identify supporting reads and run modkit.")
    parser.add_argument("bamfile", help="Input BAM file")
    parser.add_argument("vcf_file", help="Input VCF file")
    parser.add_argument("output_bamfile", help="Output BAM file")
    parser.add_argument("output_dir", help="Output directory for modkit")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file")

    args = parser.parse_args()

    sniffmeth_modkit(args.bamfile, args.vcf_file, args.output_bamfile)
    run_modkit(args.output_bamfile, args.output_dir, args.ref)

if __name__ == "__main__":
    main()
