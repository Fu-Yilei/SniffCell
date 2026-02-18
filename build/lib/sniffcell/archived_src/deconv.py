from flask import json, logging
import numpy as np
import pandas as pd
import pysam
import re
import warnings
import os
import logging
from tqdm import tqdm
from multiprocessing import Pool, current_process
from src import parse_args, calc_methylation_diff_regions
from src.region_deconv import assign_sv_to_single_celltype_per_region
from src.figures import assign_variant_with_cell_type_names, plot_cell_type_box_distributions
# Suppress FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)


def get_base_modification_dictionary_basic_supporting_reads(
    bam_file, ref_seq, chromosome, phase_region,
):
    """
        basic means for single list data structure
        Return value: a dictionary that contains cpg location and its haplotype related base modification score
    """
    methylation_identifier_0 = ('C', 0, 'm') # We only care about 5mc!!!
    methylation_identifier_1 = ('C', 1, 'm')
    phase_region_start = phase_region[0]
    phase_region_end = phase_region[1]

    phased_block_ref = ref_seq.fetch(
        chromosome, phase_region_start, phase_region_end)
    cg_loc = [
        m.start(0) for m in re.finditer("CG", str(phased_block_ref))
    ]  # Use regular expression to find all CpG locations on the reference
    # record G locatioon of 'CG's
    cg_loc = [x + phase_region_start + 1 for x in cg_loc]
    # Initialize `hp_myth_dict` as a DataFrame
    hp_myth_dict = pd.DataFrame(
        columns=cg_loc  # Set CpG locations as columns
    )    
    """
        Data structure:
            A dataframe that columns are CpG locations and rows are reads.
    """

    phased_block_alignment = bam_file.fetch(
        chromosome, phase_region[0], phase_region[1], multiple_iterators=True
    )
    read_names = []
    haplotypes = []
    read_rows = []
    for reads in phased_block_alignment:
        if not reads.is_secondary and not reads.is_supplementary:
            read_base_ref_loc = reads.get_reference_positions(full_length=True)  
            mm = (reads.modified_bases)
            if (mm != -1) and (mm != {}):
                if methylation_identifier_0 in list(mm.keys()):
                    methylation_identifier = methylation_identifier_0
                elif methylation_identifier_1 in list(mm.keys()):
                    methylation_identifier = methylation_identifier_1
                else:
                    continue
                read_row = {loc: np.nan for loc in hp_myth_dict.columns}
                for i in mm[methylation_identifier]:
                    ref_loc = read_base_ref_loc[i[0]]
                    if ref_loc:
                        if reads.is_forward:
                            mm_ref_loc = ref_loc + 1
                        else:
                            mm_ref_loc = ref_loc
                        if mm_ref_loc in hp_myth_dict.columns:
                            modification_chance = i[1]
                            read_row[mm_ref_loc] = modification_chance/255
                # Get HP tag if present
                try:
                    hp_value = reads.get_tag('HP')
                except (KeyError, AttributeError):
                    hp_value = 0
                read_names.append(reads.query_name)
                haplotypes.append(hp_value)
                read_rows.append(read_row)
    # Build DataFrame with MultiIndex
    if read_rows:
        hp_myth_dict = pd.DataFrame(read_rows, columns=hp_myth_dict.columns)
        hp_myth_dict.index = pd.MultiIndex.from_arrays([read_names, haplotypes], names=["read_name", "haplotype"])
    else:
        hp_myth_dict = pd.DataFrame(columns=hp_myth_dict.columns)

    return hp_myth_dict

def filter_cell_type_regions(cell_types, atlas, top_x=300, by='std'):
    """
        Based on the cell types in the cell_type list, filter the informative DNA methylation region in the atlas
    """
    
    atlas_df = pd.read_csv(atlas, sep='\t')
    atlas_cell_types = atlas_df[['chr', 'start', 'end'] + cell_types].copy()
    if by == 'std':
        atlas_cell_types['std'] = atlas_cell_types[cell_types].std(axis=1)
        # print(atlas_cell_types)
        top_regions = atlas_cell_types.nlargest(top_x, 'std')
        # print(top_regions)
        filtered_atlas_df = top_regions
        return filtered_atlas_df
    elif by =='diff':
        filtered_atlas_df = pd.DataFrame(columns=atlas_cell_types.columns).dropna(how='all')
        region_number_for_each_celltype = top_x // len(cell_types)
        for cell_type in cell_types:
            atlas_cell_types['selected_celltype'] = cell_type
            atlas_cell_types['mean_others'] = atlas_cell_types[cell_types].drop(columns=[cell_type]).mean(axis=1)
            atlas_cell_types['diff'] = abs(atlas_cell_types[cell_type] - atlas_cell_types['mean_others'])
            top_regions = atlas_cell_types.nlargest(region_number_for_each_celltype, 'diff')
            top_regions = top_regions.dropna(how='all')
            filtered_atlas_df = pd.concat([filtered_atlas_df, top_regions])
        return filtered_atlas_df
    else: 
        raise ValueError("Invalid catalog selection parameter. Use 'std' or 'diff'.")
    
def assign_read_to_readgroup(location, cell_types, classification_df, input_bam_file, output_bam_folder, gamma_max_confidence = 0.9):
    """
        Assign reads to read groups based on the classification results
    """
    os.makedirs(output_bam_folder, exist_ok=True)
    # for index, reads_assignemnts in classification_df.iterrows():
    chromosome = location.chr
    phase_region = [location.start, location.end]
    phased_block_alignment = input_bam_file.fetch(
    chromosome, phase_region[0], phase_region[1], multiple_iterators=True)
    cell_types_list = np.zeros(len(cell_types))
    assigned_read_count = 0
    output_bam_file = pysam.AlignmentFile(f"{output_bam_folder}/{chromosome}_{phase_region[0]}_{ phase_region[1]}.bam", "wb", template=input_bam_file)
    for reads in phased_block_alignment:
        if reads.query_name in classification_df.index:
            gamma = classification_df.loc[(reads.query_name, slice(None)), 'gamma'].values[0]
            if isinstance(gamma, (np.ndarray, list)) and len(gamma) == len(cell_types):
                if (max(gamma) > gamma_max_confidence):
                    cell_type = cell_types[np.argmax(gamma)]
                    cell_types_list[np.argmax(gamma)] += 1
                    reads.set_tag("RG", cell_type)
                    output_bam_file.write(reads)
                    assigned_read_count += 1
            else:
                reads.set_tag("RG", "low_confidence")
                output_bam_file.write(reads)
        else:
            reads.set_tag("RG", "unassigned")
            output_bam_file.write(reads)
    # print(assigned_read_count)
    output_bam_file.close()
    pysam.index(f"{output_bam_folder}/{chromosome}_{phase_region[0]}_{phase_region[1]}.bam")
    return assigned_read_count, cell_types_list



def deconvolution(output, input_bam, reference_genome, args):
    '''
        deconvolution module
        This module is used to deconvolute the cell type-specific methylation information from the
        input BAM file and the cell type atlas.
    '''
    deconv_output_location = os.path.join(output, "deconv_bam_output")
    os.makedirs(deconv_output_location, exist_ok=True)
    deconv_atlas, deconv_tissue = os.path.abspath(args.atlas), args.tissue
    deconv_region_number, deconv_method = args.region_number, args.method
    deconv_confidence, deconv_verbose = args.confidence, args.verbose
    min_hp_distance, read_fraction_threshold = args.hp_distance, args.assigned_read_fraction
    with open(os.path.join(os.path.dirname(deconv_atlas), "tissue_celltypes.json"), "r", encoding="utf-8") as f:
        celltype_dict = json.load(f)
    celltypes = celltype_dict[deconv_tissue]["cell_type"]
    celltypes_prior = None
    celltypes_prior = celltype_dict[deconv_tissue].get("prior", None)
    filtered_regions_df = filter_cell_type_regions(celltypes, deconv_atlas, deconv_region_number, by=deconv_method)
    summary_classification_df = pd.DataFrame(
        columns=["chr", "start", "end", "total_reads", "assigned_reads", "cell_type_reads_counts", "cell_type_prob_em"]
    )
    args_list = [
        (
            filtered_regions, celltypes, input_bam, deconv_output_location, reference_genome,
            deconv_confidence, celltypes_prior, deconv_verbose
        )
        for _, filtered_regions in filtered_regions_df.iterrows()
    ]

    with Pool(threads) as p:
        if deconv_method == "std":
            logging.info("Processing cell type differentially methylated regions using standard method.")
            summary_classification_series = list(
                tqdm(p.imap(process_individual_region_std_wrapper, args_list), total=len(args_list),
                        desc="Processing cell type differentially methylated regions", unit="region")
            )
        elif deconv_method == "diff":
            logging.info("Processing cell type differentially methylated regions using k-means 2 cluster mode due to the catelog selection is based on diff")
            summary_classification_series = list(
                tqdm(p.imap(process_individual_region_diff_wrapper, args_list), total=len(args_list),
                        desc="Processing cell type differentially methylated regions", unit="region")
            )

    for summary_classification in summary_classification_series:
        summary_classification_df = pd.concat([summary_classification_df, summary_classification], ignore_index=True)

    if deconv_method == "std":
        summary_classification_df = summary_classification_df[
            (summary_classification_df.max_distance <= min_hp_distance) | (summary_classification_df.max_distance.isna())
        ]
        logging.info("Filtered regions based on max distance between two haplotypes: %d regions remaining", summary_classification_df.shape[0])
        summary_classification_df = summary_classification_df[
            (summary_classification_df.assigned_reads >= read_fraction_threshold * summary_classification_df.total_reads)
        ]
        logging.info("Filtered regions based on assigned reads fraction: %d regions remaining", summary_classification_df.shape[0])
    elif deconv_method == "diff":
        average_celltype_prob_dict = summary_classification_df.groupby('selected_celltype').agg({
            'cell_type_prob': 'mean',
        }).reset_index()
        average_celltype_prob_dict = dict(zip(average_celltype_prob_dict['selected_celltype'], average_celltype_prob_dict['cell_type_prob']))
        for index, row in summary_classification_df.iterrows():
            celltype = row['selected_celltype']
            # Create a list with the current cell type's probability and the other two cell types' probabilities (set to 0)
            celltype_probs = []
            for ct in celltypes:
                if ct == celltype:
                    celltype_probs.append(row['cell_type_prob'])
                else:
                    celltype_probs.append(average_celltype_prob_dict.get(ct, 0.0))
            summary_classification_df.at[index, 'cell_type_prob_em'] = celltype_probs
        logging.info("Filtered regions based on max diff catalog selection method: %d regions remaining", summary_classification_df.shape[0])

    probs_array = np.stack(summary_classification_df["cell_type_prob_em"].values)
    nan_mask = np.isnan(probs_array).any(axis=1)
    summary_classification_df_clean = summary_classification_df[~nan_mask].copy()
    probs_array_clean = probs_array[~nan_mask]
    mean_profile_clean = probs_array_clean.mean(axis=0)
    summary_classification_df_clean["deviation_from_mean"] = np.linalg.norm(
        probs_array_clean - mean_profile_clean, axis=1
    )
    threshold = summary_classification_df_clean["deviation_from_mean"].quantile(outlier_thresold)
    summary_classification_df_filtered = summary_classification_df_clean[
        summary_classification_df_clean["deviation_from_mean"] <= threshold
    ].copy()

    summary_classification_df = summary_classification_df_filtered
    summary_classification_df.to_csv(f"{output}/summary_classification.csv", sep="\t", index=False)

    cell_type_dicts, (min_celltype_proportion, max_celltype_proportion) = assign_variant_with_cell_type_names(summary_classification_df, celltypes)
    plot_cell_type_box_distributions(cell_type_dicts, f"{output}/celltype_prediction_distributions.png")
    logging.info("Cell type proportion estimation figure saved to %s/celltype_prediction_distributions.png", output)
    return  cell_type_dicts, (min_celltype_proportion, max_celltype_proportion)