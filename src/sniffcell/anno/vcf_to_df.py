
import pandas as pd
import numpy as np
import pysam


def read_vcf_to_df(vcf_file, kanpig_read_names=None):
    """_summary_

    Args:
        vcf_file (_type_): _description_

    Returns:
        _type_: _description_
    """
    records = []
    vcf_file = pysam.VariantFile(vcf_file)
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
    if kanpig_read_names is not None:
        kanpig_df = pd.read_csv(kanpig_read_names, sep="\t", header=None, names=["sv_id", "read_name"])
        sv_to_reads = (
            kanpig_df.groupby("sv_id")["read_name"]
            .apply(list)
            .to_dict()
        )
        for record in records:
            sv_id = record["id"]
            if sv_id in sv_to_reads:
                record["supporting_reads"] = sv_to_reads[sv_id]
            else:
                record["supporting_reads"] = []
    sv_df = pd.DataFrame(records, columns=["chr", "location", "id", "sv_len", "supporting_reads", "stdev_len", "stdev_pos", "ref_start", "ref_end", 'vaf'])
    return sv_df