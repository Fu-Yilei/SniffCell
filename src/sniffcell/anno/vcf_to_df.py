
import pandas as pd
import numpy as np
import pysam


def _first_scalar(value):
    if isinstance(value, (tuple, list, np.ndarray, pd.Series)):
        if len(value) == 0:
            return pd.NA
        return value[0]
    return value


def _safe_info_get(info, key, default=pd.NA):
    try:
        value = info.get(key, default)
    except (KeyError, ValueError):
        return default
    if value is None:
        return default
    return _first_scalar(value)


def _safe_info_get_read_names(info, key="RNAMES") -> list[str]:
    try:
        value = info.get(key, ())
    except (KeyError, ValueError):
        return []
    if value is None:
        return []
    if isinstance(value, (tuple, list, np.ndarray, pd.Series)):
        out: list[str] = []
        for x in value:
            text = str(x).strip()
            if text:
                out.append(text)
        return out
    text = str(value).strip()
    return [text] if text else []


def _safe_int(value):
    try:
        if pd.isna(value):
            return np.nan
    except TypeError:
        pass
    try:
        out = float(value)
    except (TypeError, ValueError):
        return np.nan
    if not np.isfinite(out):
        return np.nan
    return int(out)


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
        svtype = str(_safe_info_get(record.info, "SVTYPE", "")).upper()
        if svtype in ["INS", "DEL"]:
            sv_len = _safe_info_get(record.info, "SVLEN", "NA")
            stdev_len = _safe_info_get(record.info, "STDEV_LEN", "NA")
            stdev_pos = _safe_info_get(record.info, "STDEV_POS", "NA")
            vaf = _safe_info_get(record.info, "VAF", pd.NA)
            try:
                vaf_missing = pd.isna(vaf)
            except TypeError:
                vaf_missing = False
            if vaf_missing:
                vaf = _safe_info_get(record.info, "AF", "NA")

            df_record = {
                "chr": record.chrom,
                "location": record.pos,
                "id": record.id,
                "sv_len": sv_len,
                "supporting_reads": _safe_info_get_read_names(record.info, "RNAMES"),
                "stdev_len": stdev_len,
                "stdev_pos": stdev_pos,
                "vaf": vaf,
            }

            stdev_pos = df_record["stdev_pos"]
            sv_len = df_record["sv_len"] if svtype == "DEL" else 0
            stdev_len = df_record["stdev_len"] 
            ref_start = df_record["location"] - stdev_pos if stdev_pos != "NA" else np.nan
            if stdev_pos != "NA":
                ref_end = df_record["location"] + stdev_pos + stdev_len - sv_len
            else:
                ref_end = np.nan

            df_record.update({"ref_start": _safe_int(ref_start), "ref_end": _safe_int(ref_end)})
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
