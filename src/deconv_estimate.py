import pysam, re
import pandas as pd
import numpy as np
from src.vcf_to_df import read_vcf_to_df
import logging, os
os.environ["HTS_LOG_LEVEL"] = "error"

logging.basicConfig(level=logging.INFO)

def find_closest_blocks(vcf_df, deconv_df):
    logging.info("Finding closest blocks...")
    closest_blocks = []
    
    for _, row in vcf_df.iterrows():
        chr_match = deconv_df[deconv_df["chr"] == row["chr"]].copy()  # Ensure it's a new DataFrame
        
        if chr_match.empty:
            closest_blocks.append(None)  # No matching chromosome
            continue

        # Compute distances to the start and end of each block
        chr_match["distance"] = chr_match.apply(lambda x: min(abs(row["location"] - x["start"]), abs(row["location"] - x["end"])), axis=1)
        
        # Find the closest row
        closest_row = chr_match.loc[chr_match["distance"].idxmin()]
        closest_blocks.append(closest_row.to_dict())

    logging.info("Finished finding closest blocks.")
    return closest_blocks

    
def assign_variant_with_cell_type_names(df, cell_type_names, epsilon=1e-6):
    logging.info("Assigning cell type names to variants...")
    assigned_cell_types = []
    confidence_scores = []
    genotype_fields = []

    for _, row in df.iterrows():
        try:
            cell_type_probs = row["closest_cell_type_prob_em"]           
            # Check if cell_type_probs is a list or array
            if not isinstance(cell_type_probs, (list, np.ndarray)):
                # logging.warning(f"Invalid data type for cell_type_probs: {type(cell_type_probs)}")
                assigned_cell_types.append("Unknown")
                confidence_scores.append(0.0)
                genotype_fields.append("./.")  # No valid assignment
                continue

            # Check if cell_type_probs is empty or all zeros
            if len(cell_type_probs) == 0 or np.all(cell_type_probs == 0):
                assigned_cell_types.append("Unknown")
                confidence_scores.append(0.0)
                genotype_fields.append("./.")  # No valid assignment
                continue

            # Ensure the number of proportions matches the number of cell types (non-zero proportions still count)
            if len(cell_type_probs) != len(cell_type_names):
                raise ValueError(f"Mismatch: Expected {len(cell_type_names)} cell types but found {len(cell_type_probs)} proportions")

            # Generate both full and half proportions for heterozygous consideration
            all_proportions = np.concatenate((cell_type_probs, cell_type_probs / 2))

            # Find the closest value (considering full and half proportions)
            closest_index = np.argmin(np.abs(all_proportions - row["vaf"]))
            closest_value = all_proportions[closest_index]

            # Compute confidence score: 1 - (relative error) for better scaling
            abs_diff = abs(row["vaf"] - closest_value)
            confidence_score = 1 - (abs_diff / (closest_value + epsilon))

            # Determine whether the selected proportion was full (homozygous) or half (heterozygous)
            cell_type_index = closest_index % len(cell_type_probs)
            is_heterozygous = closest_index >= len(cell_type_probs)

            # Assign GT field based on match
            genotype = "0/1" if is_heterozygous else "1/1"

            # Use real cell type names
            assigned_cell_types.append(cell_type_names[cell_type_index])
            confidence_scores.append(confidence_score)
            genotype_fields.append(genotype)

        except ValueError as e:
            logging.error(f"ValueError: {e}")
            raise e  # Stop execution if the number of cell types and proportions don't match
        except Exception as e:
            logging.warning(f"Exception: {e}")
            assigned_cell_types.append("Unknown")  # Handle other issues gracefully
            confidence_scores.append(None)
            genotype_fields.append("./.")  # No valid assignment

    df["assigned_cell_type"] = assigned_cell_types
    df["confidence_score"] = confidence_scores
    df["GT"] = genotype_fields
    logging.info("Finished assigning cell type names.")
    return df

def annotate_vcf_by_id_copy(vcf_path, annotated_df, output_vcf_path):
    logging.info(f"Annotating VCF file: {vcf_path}")
    vcf_in = pysam.VariantFile(vcf_path, "r")
    annotated_df = annotated_df.dropna(subset=["id"]) # Solve the issue of empty id
    header = vcf_in.header.copy()
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

def estimate_celltype_assignment(vcf_file, sv_methylation_df, deconv_df, cell_type_list, output_vcf):
    logging.info("Starting cell type assignment estimation...")
    # vcf_file_df = read_vcf_to_df(pysam.VariantFile(vcf_file))
    closest_blocks = find_closest_blocks(sv_methylation_df, deconv_df)
    closest_blocks = [block if block is not None else {} for block in closest_blocks]

    closest_df = pd.DataFrame(closest_blocks)
    closest_df = closest_df.add_prefix("closest_")
    merged_df = pd.concat([sv_methylation_df, closest_df], axis=1)

    updated_df = assign_variant_with_cell_type_names(merged_df, cell_type_list)
    updated_df.to_csv(os.path.join(os.path.dirname(output_vcf), "sv_methylation_celltype_estimation.csv"))
    annotate_vcf_by_id_copy(vcf_file, updated_df, output_vcf)
    logging.info("Finished cell type assignment estimation.")
