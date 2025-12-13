import os
import logging
import numpy as np
import pandas as pd
import pysam
from multiprocessing import current_process
from src.deconv_em import em_haplotype_and_combined
from src.deconv.deconv import (
    get_base_modification_dictionary_basic_supporting_reads,
    filter_cell_type_regions,
    assign_read_to_readgroup,
)
from src.region_deconv import kmeans_per_haplotype
def process_individual_region_std(
    filtered_regions, celltypes, bam_file_path, output, ref_seq_path, threshold, celltype_prior, verbose=False
):
    current_process().name = "sniffcell"
    chromosome, phase_region = filtered_regions.chr, [filtered_regions.start, filtered_regions.end]

    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    ref_seq = pysam.FastaFile(ref_seq_path)
    reads_methylation_df = get_base_modification_dictionary_basic_supporting_reads(
        bam_file, ref_seq, chromosome, phase_region
    )

    if reads_methylation_df.empty:
        logging.error(f"No reads found in the region {chromosome}:{phase_region[0]}-{phase_region[1]}")
        return None

    p_init_region = filtered_regions[celltypes]
    alpha_init = celltype_prior if celltype_prior else np.ones(len(celltypes)) / len(celltypes)

    alphas, _, gamma_df = em_haplotype_and_combined(
        df=reads_methylation_df, alpha_init=alpha_init, p_init=p_init_region, max_iter=100, tol=1e-5, random_state=42
    )

    max_distance, alpha_final, gamma_final = None, None, None
    if alphas[0] is None and alphas[1] is None:
        alpha_final = alphas[2]
        gamma_final = gamma_df.copy()
        gamma_final["gamma"] = gamma_df["gamma_list"].apply(
            lambda x: x[-1] if isinstance(x, list) and len(x) > 0 else None
        )
    elif len(alphas) == 3:
        alpha_final = np.mean([i for i in alphas if i is not None], axis=0)
        gamma_final = gamma_df.copy()
        gamma_final["gamma"] = gamma_df["gamma_list"].apply(
            lambda x: np.mean([i for i in x if i is not None], axis=0)
        )
        if alphas[0] is not None and alphas[1] is not None:
            max_distance = max(
                np.linalg.norm(alphas[0] - alphas[1]),
                np.linalg.norm(alphas[0] - alphas[2]),
                np.linalg.norm(alphas[1] - alphas[2]),
            )

    if verbose:
        gamma_final.to_csv(f"{output}/{chromosome}_{phase_region[0]}_{phase_region[1]}.tsv", sep="\t")

    assigned_reads_count, cell_types_list = assign_read_to_readgroup(
        filtered_regions, celltypes, gamma_final, bam_file, output, gamma_max_confidence=threshold
    )

    return pd.DataFrame(
        [
            {
                "chr": chromosome,
                "start": phase_region[0],
                "end": phase_region[1],
                "total_reads": reads_methylation_df.shape[0],
                "assigned_reads": assigned_reads_count,
                "cell_type_reads_counts": cell_types_list,
                "cell_type_prob_em": alpha_final,
                "max_distance": max_distance,
            }
        ]
    )






def process_individual_region_diff(
    filtered_regions, celltypes, bam_file_path, output, ref_seq_path, threshold, celltype_prior, verbose=False
):
    current_process().name = "sniffcell"
    chromosome, phase_region = filtered_regions.chr, [filtered_regions.start, filtered_regions.end]
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    ref_seq = pysam.FastaFile(ref_seq_path)
    reads_methylation_df = get_base_modification_dictionary_basic_supporting_reads(
        bam_file, ref_seq, chromosome, phase_region
    )
    if reads_methylation_df.empty:
        logging.error(f"No reads found in the region {chromosome}:{phase_region[0]}-{phase_region[1]}")
        return None
    else:
        cpg_start, cpg_end = reads_methylation_df.columns[0], reads_methylation_df.columns[-1]
        kmeans_per_haplotype_df = kmeans_per_haplotype(
            reads_methylation_df, cpg_start, cpg_end, n_clusters=2,
        ).set_index(['read_name', 'haplotype'])
        p_init_region = filtered_regions[celltypes]
        reads_methylation_df['cluster'] = kmeans_per_haplotype_df['cluster']
        total_assigned_reads = 0
        reads_assigned_to_selected_celltype_list = []
        reads_assigned_to_others_list = []
        hp_num = 0
        for hap, hap_df in reads_methylation_df.groupby(level='haplotype'):
            if hap == 0:
                continue  # Skip unphased
            features = hap_df.drop(columns='cluster')
            clusters = hap_df['cluster']
            centroids = features.groupby(clusters).mean().dropna(axis=1)
            if centroids.shape[0] < 2:
                logging.warning(f"Not enough clusters for haplotype {hap} in region {chromosome}:{phase_region[0]}-{phase_region[1]}")
                continue  # Only one cluster → can't compute distance
            mean_others = filtered_regions.mean_others
            selected_celltype_name = filtered_regions['selected_celltype']
            selected_celltype_prob = filtered_regions[selected_celltype_name]

            # Compute distances from centroids to selected_celltype_prob and mean_others
            distances_to_selected = centroids.apply(lambda x: np.linalg.norm(x.values - selected_celltype_prob), axis=1)
            distances_to_others = centroids.apply(lambda x: np.linalg.norm(x.values - mean_others), axis=1)

            # Find closest clusters
            closest_to_selected = distances_to_selected.idxmin()
            closest_to_others = distances_to_others.idxmin()
            if distances_to_selected[closest_to_selected] > 1 or distances_to_others[closest_to_others] > 1:
                logging.warning(f"Region {chromosome}:{phase_region[0]}-{phase_region[1]}, haplotype {hap} is discarded due to large centroid distance: {distances_to_selected[closest_to_selected]} (selected) and {distances_to_others[closest_to_others]} (others)")
                continue
            else:
                reads_assigned_to_selected_celltype = hap_df[hap_df['cluster'] == closest_to_selected]
                reads_assigned_to_others = hap_df[hap_df['cluster'] == closest_to_others] 
                total_assigned_reads += hap_df.shape[0]
                reads_assigned_to_selected_celltype_list += list(reads_assigned_to_selected_celltype.index.get_level_values('read_name').unique())
                reads_assigned_to_others_list += list(reads_assigned_to_others.index.get_level_values('read_name').unique())
                hp_num += 1

        # Compute proportion difference among two haplotypes
        if hp_num == 0:
            logging.warning(f"No reads assigned to the selected cell type in region {chromosome}:{phase_region[0]}-{phase_region[1]}")
            return None
        elif hp_num == 2:
            prop_selected = len(reads_assigned_to_selected_celltype_list) / total_assigned_reads if total_assigned_reads > 0 else 0
            prop_others = len(reads_assigned_to_others_list) / total_assigned_reads if total_assigned_reads > 0 else 0
            max_distance = abs(prop_selected - prop_others)
        else:
            prop_selected = len(reads_assigned_to_selected_celltype_list) / total_assigned_reads if total_assigned_reads > 0 else 0
            prop_others = len(reads_assigned_to_others_list) / total_assigned_reads if total_assigned_reads > 0 else 0
            max_distance = None

        return pd.DataFrame(
            [
            {
                "chr": chromosome,
                "start": phase_region[0],
                "end": phase_region[1],
                "total_reads": reads_methylation_df.shape[0],
                "assigned_reads": total_assigned_reads,
                "cell_type_reads_counts": [len(reads_assigned_to_selected_celltype_list), len(reads_assigned_to_others_list)],
                "cell_type_prob": len(reads_assigned_to_selected_celltype_list) / total_assigned_reads if total_assigned_reads > 0 else 0,
                "max_distance": max_distance,
                "selected_celltype": filtered_regions['selected_celltype'],
            }
            ]
        )