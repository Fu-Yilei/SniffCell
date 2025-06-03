import sys, os
import shutil
import subprocess
import pandas as pd
import glob, logging


def run_wgbstools(wgbs_path, bam_file_path, output_dir, atlas_tsv_path, threads):
    """
    Run wgbstools to generate methylation regions.
    """
    # Convert atlas_tsv_path to BED format
    bam_file_prefix = os.path.splitext(os.path.basename(bam_file_path))[0]
    wgbstools_output = os.path.join(output_dir, f"{bam_file_prefix}.pat.gz")
    if os.path.exists(wgbstools_output):
        logging.warning(f"wgbs ouptut directory {output_dir} already exists. It will be removed.")
        return 
    
    atlas_prefix = os.path.splitext(os.path.basename(atlas_tsv_path))[0]
    atlas_bed = os.path.join(output_dir, f"{atlas_prefix}.bed")

    with open(atlas_tsv_path, "r") as tsv_file, open(atlas_bed, "w") as bed_file:
        next(tsv_file)  # Skip the header line
        for line in tsv_file:
            columns = line.strip().split("\t")
            chrom = columns[0]
            start = int(columns[1]) - 1  # Convert to 0-based indexing
            end = columns[2]
            name = columns[5]
            bed_file.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t.\n")
    # output_pat_path = os.path.join(output_dir, "wgbstools.pat.gz")
    subprocess.run(
        [
            wgbs_path, "bam2pat", "-np",
            "-L", atlas_bed,
            "-@", str(threads),
            "--genome", "hg38",
            "--force",
            "--no_beta",
            "-T", output_dir,
            "--out_dir", output_dir,
            bam_file_path
        ],
        check=True
    )
    # return output_pat_path


def run_uxm(uxm_path, file_folder, threads, atlas_file):
    """
    Run uxm to generate methylation regions.
    """
    pat_file_path = glob.glob(os.path.join(file_folder, "*.pat.gz"))[0]
    output_uxm_path = os.path.join(file_folder, "uxm_output.tsv")
    subprocess.run(
        [
            uxm_path, "deconv",
            "--atlas", atlas_file,
            "--output", output_uxm_path,
            "--threads", str(threads),
            f"{pat_file_path}"
        ],
        check=True
    )
    return output_uxm_path

def get_uxm_prior(uxm_output_path, selected_cell_types):
    uxm_df = pd.read_csv(uxm_output_path,  index_col="CellType")
    celltype_proportion_dict = dict(zip(uxm_df.index, uxm_df[uxm_df.columns[0]]))
    proportion_selected_celltypes = {ct: celltype_proportion_dict[ct] for ct in selected_cell_types if ct in celltype_proportion_dict}
    if not proportion_selected_celltypes:
        raise ValueError("No selected cell types found in the uxm output.")
    total_proportion = sum(proportion_selected_celltypes.values())
    if total_proportion == 0:
        raise ValueError("Total proportion of selected cell types is zero.")
    normalized_proportions = {ct: prop / total_proportion for ct, prop in proportion_selected_celltypes.items()} # Normalize proportions to sum to 1
    return normalized_proportions