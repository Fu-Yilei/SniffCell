
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
    sv_df = pd.DataFrame(records, columns=["chr", "location", "id", "sv_len", "supporting_reads", "stdev_len", "stdev_pos", "ref_start", "ref_end"])
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
