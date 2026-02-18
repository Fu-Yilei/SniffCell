import pandas as pd
import numpy as np
from tqdm import tqdm
import pysam
import multiprocessing
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO)

def filter_snp_df(snp_df, min_confidence):
    """
    Filter SNP DataFrame based on minimum confidence and assigned read fraction.
    """
    filtered_df = snp_df[
        (snp_df["confidence_score"] >= min_confidence) &
        (
            ((snp_df.vaf > 0.0) & (snp_df.vaf < 0.4)) |
            ((snp_df.vaf > 0.6) & (snp_df.vaf < 0.8))
        )
    ].copy()
    logging.info("Filtered SNP number based on minimum confidence: %d total, %d rows remaining", snp_df.shape[0], filtered_df.shape[0])
    return filtered_df.reset_index(drop=True)

def process_chromosome(args):
    import sys
    chrom, chrom_snp_df, bam_file = args
    snp_pos_map = {row["location"]: (row["ref"], row["alt"], row["id"]) for _, row in chrom_snp_df.iterrows()}
    snp_read_ids = {pos: {"ref_reads": set(), "alt_reads": set()} for pos in snp_pos_map}
    # Get number of reads in chrom (optional, for progress bar)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # If .count() is too slow, just set total=None for tqdm (no length bar)
        try:
            total_reads = bam.count(contig=chrom)
        except Exception:
            total_reads = None
        with tqdm(
            bam.fetch(chrom),
            total=total_reads,
            desc=f"Chr {chrom}",
            position=0 if multiprocessing.current_process().name == "MainProcess" else None,
            leave=False,
            file=sys.stdout
        ) as pbar:
            for read in pbar:
                if read.is_unmapped or read.query_sequence is None:
                    continue
                ref_positions = read.get_reference_positions(full_length=True)
                for idx, ref_pos in enumerate(ref_positions):
                    if ref_pos in snp_pos_map and idx < len(read.query_sequence):
                        ref, alt, snp_id = snp_pos_map[ref_pos]
                        base = read.query_sequence[idx]
                        if base == ref:
                            snp_read_ids[ref_pos]["ref_reads"].add(read.query_name)
                        if base == alt:
                            snp_read_ids[ref_pos]["alt_reads"].add(read.query_name)
    # Package results as DataFrame with the same order as chrom_snp_df
    ref_read_ids = []
    alt_read_ids = []
    for _, row in chrom_snp_df.iterrows():
        pos = row["location"]
        ref_read_ids.append(list(snp_read_ids[pos]["ref_reads"]))
        alt_read_ids.append(list(snp_read_ids[pos]["alt_reads"]))
    chrom_snp_df = chrom_snp_df.copy()
    chrom_snp_df["ref_read_ids"] = ref_read_ids
    chrom_snp_df["alt_read_ids"] = alt_read_ids
    return chrom_snp_df


def get_snp_supporting_read_ids_parallel(snp_df, bam_file, nproc=4):
    """
    Main function: batch by chromosome, process each chromosome's SNPs in parallel.
    """
    chrom_groups = [ (chrom, group, bam_file) for chrom, group in snp_df.groupby("chr") ]
    if nproc > 1:
        with multiprocessing.Pool(nproc) as pool:
            results = list(tqdm(pool.imap(process_chromosome, chrom_groups), total=len(chrom_groups)))
    else:
        results = [process_chromosome(args) for args in tqdm(chrom_groups)]
    # Combine results back to a single DataFrame in correct order
    result_df = pd.concat(results, ignore_index=True)
    # Reorder to match original snp_df if needed
    result_df = result_df.reindex(snp_df.index)
    return result_df