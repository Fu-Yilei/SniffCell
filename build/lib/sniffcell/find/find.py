import json
import pandas as pd
import numpy as np

from sniffcell.find import ctdmr
from sniffcell.find.ctdmr import means_from_mapping
import logging
def find_main(args):

    logger = logging.getLogger(__name__)
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    logger.info("Starting find_main")
    npy_file = args.npy
    index_file = args.index
    meta_file = args.meta
    celltypes_file = args.celltypes_file
    celltupes_keys = args.celltypes_keys

    logger.info(f"Loading npy blocks from {npy_file}")
    all_celltype_blocks = np.load(npy_file)
    logger.info(f"Loaded block array with shape {getattr(all_celltype_blocks, 'shape', None)}")

    all_celltypes = []
    logger.info(f"Reading celltype names from {meta_file}")
    with open(meta_file, "r", encoding="utf-8") as meta_f:
        all_celltypes = [line.strip() for line in meta_f if line.strip()]
    logger.info(f"Loaded {len(all_celltypes)} celltype names")

    logger.info(f"Reading CpG index from {index_file}")
    cpg_index = pd.read_csv(index_file, sep="\t", header=None, names=['chr', 'start', 'end','startCpG','endCpG'])
    logger.info(f"CpG index rows: {len(cpg_index)}")

    logger.info("Building methylation DataFrame")
    M_df  = pd.DataFrame(all_celltype_blocks, columns=all_celltypes, index=cpg_index.index)
    logger.info(f"M_df shape: {M_df.shape}")

    logger.info(f"Loading celltype mapping '{celltupes_keys}' from {celltypes_file}")
    with open(celltypes_file, "r", encoding="utf-8") as celltypes_f:
        major_to_subtypes = json.load(celltypes_f)[celltupes_keys]
    logger.info(f"Loaded mapping with {len(major_to_subtypes)} groups")

    logger.info("Computing group means")
    mean_by_group = means_from_mapping(M_df, major_to_subtypes)
    if hasattr(mean_by_group, "shape"):
        logger.info(f"mean_by_group shape: {mean_by_group.shape}")

    logger.info(f"Calling ctdmr.call_ct_specific_dmrs with idx_df={cpg_index.shape}, mean_by_group={len(mean_by_group)} groups, diff_threshold={args.diff_threshold}, min_rows={args.min_rows}")
    dmrs = ctdmr.call_ct_specific_dmrs(
        idx_df=cpg_index,
        mean_by_group=mean_by_group,
        diff_threshold=args.diff_threshold,
        min_rows=args.min_rows,
        min_cpgs=args.min_cpgs,
        min_bp=0,
        direction='both',
        max_gap_bp=args.max_gap_bp,
        bed_out=args.output,
    )
    if hasattr(dmrs, "shape"):
        logger.info(f"DMRs result shape: {dmrs.shape}")
    else:
        try:
            logger.info(f"DMRs result length: {len(dmrs)}")
        except Exception:
            logger.info(f"DMRs result type: {type(dmrs).__name__}")

    logger.info("find_main completed")
    return dmrs