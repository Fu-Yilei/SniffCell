import pysam, re
import pandas as pd
import numpy as np
from src.anno.vcf_to_df import read_vcf_to_df
import logging, os
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial

os.environ["HTS_LOG_LEVEL"] = "error"

logging.basicConfig(level=logging.INFO)

def find_closest_block_for_row(row, deconv_df, global_proportions):
    """
    Find the closest upstream and/or downstream block in `deconv_df` for a given genomic row,
    and compute cell type proportions by averaging with global proportions as needed.
    Parameters
    ----------
    row : pandas.Series
        A row containing at least 'chr' and 'location' fields, representing the genomic position of interest.
    deconv_df : pandas.DataFrame
        DataFrame containing block information with columns including 'chr', 'start', 'end', and 'cell_type_prob_em'.
    global_proportions : array-like
        Global cell type proportions to use as fallback or for averaging.
    Returns
    -------
    dict
        A dictionary containing the closest block's metadata and computed cell type proportions:
            - 'chr': Chromosome of the block or 'global' if no block found.
            - 'start': Block start position or None.
            - 'end': Block end position or None.
            - 'total_reads': Total reads in the block or None.
            - 'assigned_reads': Assigned reads in the block or None.
            - 'cell_type_reads_counts': Cell type read counts (array) or zeros if global.
            - 'cell_type_prob_em': Averaged cell type proportions (array).
            - 'max_distance': Maximum distance in the block or None.
            - 'deviation_from_mean': Deviation from mean or None.
            - 'distance': Distance from the row's location to the closest block boundary or None.
    Notes
    -----
    - If no block is found on the same chromosome, returns a dictionary with 'global' values and global proportions.
    - If only upstream or downstream block exists, averages its proportions with global.
    - If both exist, averages both blocks' and global proportions, and uses metadata from the closer block.
    """


    chr_match = deconv_df[deconv_df["chr"] == row["chr"]].copy()
    if chr_match.empty:
        return {
            'chr': 'global',
            'start': None,
            'end': None,
            'total_reads': None,
            'assigned_reads': None,
            'cell_type_reads_counts': np.zeros(len(global_proportions)),
            'cell_type_prob_em': global_proportions,
            'max_distance': None,
            'deviation_from_mean': None,
            'distance': None
        }

    location = row["location"]
    # Upstream: block with end <= location
    upstream_blocks = chr_match[chr_match["end"] <= location]
    # Downstream: block with start >= location
    downstream_blocks = chr_match[chr_match["start"] >= location]

    closest_upstream = None
    closest_downstream = None

    if not upstream_blocks.empty:
        upstream_idx = (location - upstream_blocks["end"]).abs().idxmin()
        closest_upstream = upstream_blocks.loc[upstream_idx]
    if not downstream_blocks.empty:
        downstream_idx = (downstream_blocks["start"] - location).abs().idxmin()
        closest_downstream = downstream_blocks.loc[downstream_idx]

    # If both are missing, fallback to global
    if closest_upstream is None and closest_downstream is None:
        return {
            'chr': 'global',
            'start': None,
            'end': None,
            'total_reads': None,
            'assigned_reads': None,
            'cell_type_reads_counts': np.zeros(len(global_proportions)),
            'cell_type_prob_em': global_proportions,
            'max_distance': None,
            'deviation_from_mean': None,
            'distance': None
        }

    # If only one exists, use its proportions and average with global
    if closest_upstream is not None and closest_downstream is None:
        props = (np.array(closest_upstream["cell_type_prob_em"]) + np.array(global_proportions)) / 2
        result = closest_upstream.to_dict()
        result["cell_type_prob_em"] = props
        result["distance"] = abs(location - closest_upstream["end"])
        return result
    if closest_downstream is not None and closest_upstream is None:
        props = (np.array(closest_downstream["cell_type_prob_em"]) + np.array(global_proportions)) / 2
        result = closest_downstream.to_dict()
        result["cell_type_prob_em"] = props
        result["distance"] = abs(location - closest_downstream["start"])
        return result

    # If both exist, average their proportions and with global
    props = (np.array(closest_upstream["cell_type_prob_em"]) +
             np.array(closest_downstream["cell_type_prob_em"]) +
             np.array(global_proportions)) / 3
    # Choose the closer block for metadata
    dist_up = abs(location - closest_upstream["end"])
    dist_down = abs(location - closest_downstream["start"])
    if dist_up <= dist_down:
        result = closest_upstream.to_dict()
        result["distance"] = dist_up
    else:
        result = closest_downstream.to_dict()
        result["distance"] = dist_down
    result["cell_type_prob_em"] = props
    return result

def find_closest_blocks(vcf_df, deconv_df, global_proportions, threads=1):
    logging.info("Finding closest blocks in parallel...")
    # Prepare the function so Pool only needs to pass a single row
    worker = partial(find_closest_block_for_row, deconv_df=deconv_df, global_proportions=global_proportions)
    # Convert DataFrame rows to dictionaries so they are picklable
    rows = [row._asdict() if hasattr(row, '_asdict') else row.to_dict() for _, row in vcf_df.iterrows()]
    with Pool(threads) as pool:
        closest_blocks = list(tqdm(
            pool.imap(worker, rows),
            desc="Finding closest blocks",
            total=len(vcf_df)
        ))
    logging.info("Finished finding closest blocks.")
    return closest_blocks
    


def _assign_variant_row(row, cell_type_names, epsilon=1e-6):
    try:
        cell_type_probs = row["closest_cell_type_prob_em"]
        if not isinstance(cell_type_probs, (list, np.ndarray)):
            return ("Unknown", 0.0, "./.")
        if len(cell_type_probs) == 0 or np.all(cell_type_probs == 0):
            return ("Unknown", 0.0, "./.")
        if len(cell_type_probs) != len(cell_type_names):
            raise ValueError(f"Mismatch: Expected {len(cell_type_names)} cell types but found {len(cell_type_probs)} proportions")
        all_proportions = np.concatenate((cell_type_probs, np.array(cell_type_probs) / 2))
        closest_index = np.argmin(np.abs(all_proportions - row["vaf"]))
        closest_value = all_proportions[closest_index]
        abs_diff = abs(row["vaf"] - closest_value)
        confidence_score = 1 - (abs_diff / (closest_value + epsilon))
        cell_type_index = closest_index % len(cell_type_probs)
        is_heterozygous = closest_index >= len(cell_type_probs)
        genotype = "0/1" if is_heterozygous else "1/1"
        return (cell_type_names[cell_type_index], confidence_score, genotype)
    except ValueError as e:
        logging.error(f"ValueError: {e}")
        raise e
    except Exception as e:
        logging.warning(f"Exception: {e}")
        return ("Unknown", None, "./.")

def assign_variant_with_cell_type_names(df, cell_type_names, global_proportions, epsilon=1e-6, threads=1):
    logging.info("Assigning cell type names to variants (multiprocessing)...")
    rows = [row for _, row in df.iterrows()]
    worker = partial(_assign_variant_row, cell_type_names=cell_type_names, epsilon=epsilon)
    with Pool(threads) as pool:
        results = list(tqdm(pool.imap(worker, rows), desc="Assigning cell types", total=len(rows)))

    assigned_cell_types, confidence_scores, genotype_fields = zip(*results)
    df["assigned_cell_type"] = assigned_cell_types
    df["confidence_score"] = confidence_scores
    df["GT"] = genotype_fields
    logging.info("Finished assigning cell type names.")
    return df

def annotate_vcf_by_id_copy(vcf_path, annotated_df, output_vcf_path, cmd_info):
    logging.info(f"Annotating VCF file: {vcf_path}")
    vcf_in = pysam.VariantFile(vcf_path, "r")
    annotated_df = annotated_df.dropna(subset=["id"]) # Solve the issue of empty id
    header = vcf_in.header.copy()
    sniffcell_version = cmd_info[0]
    command_line = cmd_info[1]
    header.add_line(f'##source=SniffCell {sniffcell_version}')
    header.add_line(f'##command={command_line}')

    header.add_line('##INFO=<ID=CELLTYPE,Number=1,Type=String,Description="Predicted cell type for the variant">')
    header.add_line('##INFO=<ID=CONF,Number=1,Type=Float,Description="Confidence score of cell type assignment">')
    header.add_line('##INFO=<ID=HET_STATUS,Number=1,Type=String,Description="Heterozygous (0/1) or Homozygous (1/1) status">')
    header.add_line('##INFO=<ID=CLOSEST_REGION,Number=1,Type=String,Description="Closest related genomic region">')
    
    vcf_out = pysam.VariantFile(output_vcf_path, "w", header=header)
    
    annotation_dict = annotated_df.set_index("id").to_dict(orient="index")
    
    for record in vcf_in.fetch():
        new_record = record
        
        if record.id in annotation_dict:
            match_row = annotation_dict[record.id]
            
            predicted_cell_type = match_row.get("assigned_cell_type", ".")
            conf_value = match_row.get("confidence_score")
            confidence_score = round(conf_value, 3) if conf_value is not None else "."
            het_status = "HET" if match_row.get("GT") == "0/1" else "HOM"
            closest_chr = match_row.get('closest_chr', '.')
            closest_start = match_row.get('closest_start', None)
            closest_end = match_row.get('closest_end', None)
            if closest_start is None or pd.isna(closest_start) or closest_end is None or pd.isna(closest_end):
                closest_region = "."
            else:
                closest_region = f"{closest_chr}_{int(closest_start)}_{int(closest_end)}"
            new_record.header.add_line('##INFO=<ID=CELLTYPE,Number=1,Type=String,Description="Predicted cell type for the variant">')
            new_record.header.add_line('##INFO=<ID=CONF,Number=1,Type=Float,Description="Confidence score of cell type assignment">')
            new_record.header.add_line('##INFO=<ID=HET_STATUS,Number=1,Type=String,Description="Heterozygous (0/1) or Homozygous (1/1) status">')
            new_record.header.add_line('##INFO=<ID=CLOSEST_REGION,Number=1,Type=String,Description="Closest related genomic region">')
            new_record.info["CELLTYPE"] = predicted_cell_type
            new_record.info["CONF"] = confidence_score
            new_record.info["HET_STATUS"] = het_status
            new_record.info["CLOSEST_REGION"] = closest_region
        
        vcf_out.write(new_record)
    vcf_in.close()
    vcf_out.close()
    logging.info(f"Annotated VCF file saved to: {output_vcf_path}")

def calculate_global_celltype_proportion(deconv_df, deconv_region_number, cell_type_list, assignment_method):
    logging.info("Calculating global cell type proportions...")
    if assignment_method == 'std':
        cell_type_proportions = np.mean(deconv_df['cell_type_prob_em'], axis=0)
    elif assignment_method == "diff":
        region_num_for_each_celltype = deconv_region_number // len(cell_type_list) #TODO: need to only access one cell type for the diff method
        cell_type_proportions = np.mean(deconv_df['cell_type_prob_em'], axis=0)
    else:
        raise ValueError(f"Unknown assignment method: {assignment_method}")
    return cell_type_proportions


def estimate_celltype_assignment(vcf_file, sv_methylation_df, deconv_df, cell_type_list, output_vcf, cmd_info, threads, assignment_method="std"):
    logging.info("Starting cell type assignment estimation...")
    global_proportions = calculate_global_celltype_proportion(deconv_df, len(deconv_df), cell_type_list, assignment_method)
    logging.info("Calculating global cell type proportions: %s", global_proportions)
    closest_blocks = find_closest_blocks(sv_methylation_df, deconv_df, global_proportions=global_proportions, threads=threads )
    closest_blocks = [block if block is not None else {} for block in closest_blocks]
    closest_df = pd.DataFrame(closest_blocks)
    closest_df = closest_df.add_prefix("closest_")
    merged_df = pd.concat([sv_methylation_df, closest_df], axis=1)
    updated_df = assign_variant_with_cell_type_names(merged_df, cell_type_list, global_proportions=global_proportions, epsilon=1e-6, threads=threads)
    updated_df.to_csv(os.path.join(os.path.dirname(output_vcf), "sv_methylation_celltype_estimation.csv"))
    annotate_vcf_by_id_copy(vcf_file, updated_df, output_vcf, cmd_info)
    logging.info("Finished cell type assignment estimation.")


def estimate_celltype_assignment_snp(vcf_file, snp_df, deconv_df, cell_type_list, output_vcf, cmd_info, threads, assignment_method="std"):
    logging.info("Starting cell type assignment estimation...")
    global_proportions = calculate_global_celltype_proportion(deconv_df, len(deconv_df), cell_type_list, assignment_method)
    logging.info("Calculating global cell type proportions: %s", global_proportions)
    # For SNPs, always use global proportions (do not use region-specific proportions)
    closest_blocks = [{"chr": "global",
                       "start": None,
                       "end": None,
                       "total_reads": None,
                       "assigned_reads": None,
                       "cell_type_reads_counts": np.zeros(len(global_proportions)),
                       "cell_type_prob_em": global_proportions,
                       "max_distance": None,
                       "deviation_from_mean": None,
                       "distance": None
                      } for _ in range(len(snp_df))]
    closest_blocks = [block if block is not None else {} for block in closest_blocks]
    closest_df = pd.DataFrame(closest_blocks)
    closest_df = closest_df.add_prefix("closest_")
    merged_df = pd.concat([snp_df, closest_df], axis=1)
    updated_df = assign_variant_with_cell_type_names(merged_df, cell_type_list, global_proportions=global_proportions, epsilon=1e-6, threads=threads)
    updated_df.to_csv(os.path.join(os.path.dirname(output_vcf), "snp_methylation_celltype_estimation.csv"))
    return updated_df