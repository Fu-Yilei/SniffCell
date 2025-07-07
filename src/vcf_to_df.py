
import pandas as pd
import numpy as np

def read_vcf_to_df(vcf_file):
    """_summary_

    Args:
        vcf_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    records = []
    for record in vcf_file.fetch():
        if record.info["SVTYPE"] in ["INS", "DEL"]:
            df_record = {
                "chr": record.chrom,
                "location": record.pos,
                "id": record.id,
                "sv_len": record.info.get("SVLEN", "NA"),
                "supporting_reads": record.info.get("RNAMES", "NA"),
                "stdev_len": record.info.get("STDEV_LEN", "NA"),
                "stdev_pos": record.info.get("STDEV_POS", "NA"),
                "vaf": record.info.get("VAF") if "VAF" in record.info else record.info.get("AF", "NA"),
            }

            stdev_pos = df_record["stdev_pos"]
            sv_len = df_record["sv_len"] if record.info["SVTYPE"] == "DEL" else 0
            stdev_len = df_record["stdev_len"] 
            ref_start = df_record["location"] - stdev_pos if stdev_pos != "NA" else np.nan
            if stdev_pos != "NA":
                ref_end = df_record["location"] + stdev_pos + stdev_len - sv_len
            else:
                ref_end = np.nan

            df_record.update({"ref_start": int(ref_start), "ref_end": int(ref_end)})
            records.append(df_record)
    sv_df = pd.DataFrame(records, columns=["chr", "location", "id", "sv_len", "supporting_reads", "stdev_len", "stdev_pos", "ref_start", "ref_end", 'vaf'])
    return sv_df



def read_vcf_to_df_comp(vcf_file):
    """_summary_

    Args:
        vcf_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    records = []
    sample_list = ["HG002", "HG00438", "HG005", "HG02257", "HG02486", "HG02622"]
    for record in vcf_file.fetch():
        if record.info["SVTYPE"] in ["INS", "DEL"]:
            df_record = {
                "chr": record.chrom,
                "location": record.pos,
                "id": record.id,
                "sv_len": record.info.get("SVLEN", "NA"),}
            for sample_name in sample_list:
                if record.samples[sample_name]['GT'] != (0,0):
                    df_record.update({'sample_name' : sample_name})
            records.append(df_record)
    sv_df = pd.DataFrame(records, columns=["chr", "location", "id", "sv_len", "sample_name"])
    return sv_df



def read_snp_vcf_to_df(vcf_in):
    sample_names = list(vcf_in.header.samples)
    records = []

    for record in vcf_in.fetch():
        # Only keep SNVs: single base ref and all alts must also be single base
        if len(record.ref) != 1:
            continue
        if not record.alts:
            continue
        alt_alleles = [alt for alt in record.alts if alt and len(alt) == 1]
        if not alt_alleles:
            continue

        for alt in alt_alleles:
            var_id = f"{record.chrom}_{record.pos}_{record.ref}_{alt}"

            for sample in sample_names:
                af = record.samples[sample].get("AF")
                if isinstance(af, (list, tuple)) and len(af) == len(record.alts):
                    # Match AF to correct ALT index
                    try:
                        alt_index = record.alts.index(alt)
                        af_val = af[alt_index]
                    except ValueError:
                        af_val = None
                else:
                    af_val = af[0] if isinstance(af, (list, tuple)) else af

                # Get phase block if available
                phaseblock = record.samples[sample].get("PS", None)

                records.append({
                    "chr": record.chrom,
                    "location": record.pos,
                    "id": var_id,
                    "ref": record.ref,
                    "alt": alt,
                    "sample": sample,
                    "quality": record.qual,
                    "vaf": af_val,
                    'phaseblock': phaseblock
                })

    return pd.DataFrame(records)
