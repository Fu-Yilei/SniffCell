
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


def _normalize_format_value(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (tuple, list, np.ndarray, pd.Series)):
        return [_normalize_format_value(x) for x in value]
    return value


def _gt_tuple_to_str(gt_value):
    if gt_value is None:
        return "./."
    alleles = []
    for allele in gt_value:
        if allele is None:
            alleles.append(".")
        else:
            alleles.append(str(allele))
    if not alleles:
        return "./."
    return "/".join(alleles)


def _extract_sample_assignments(record):
    sample_assignments = {}
    called_samples = []
    nonref_samples = []
    alt_support_samples = []

    for sample_name in record.samples:
        sample_data = dict(record.samples[sample_name].items())
        gt_value = sample_data.get("GT")
        gt_str = _gt_tuple_to_str(gt_value)
        ad_value = sample_data.get("AD")
        ad_list = _normalize_format_value(ad_value)
        alt_depth = 0
        if isinstance(ad_list, list) and len(ad_list) > 1:
            alt_depth = int(
                sum(
                    int(x)
                    for x in ad_list[1:]
                    if x is not None and not pd.isna(x)
                )
            )

        is_called = gt_str != "./."
        is_nonref = is_called and any((allele is not None and allele != 0) for allele in (gt_value or ()))
        has_alt_support = alt_depth > 0

        assignment = {
            "GT": gt_str,
            "is_called": is_called,
            "is_nonref": is_nonref,
            "has_alt_support": has_alt_support,
            "alt_depth": alt_depth,
        }
        for key, value in sample_data.items():
            if key == "GT":
                continue
            assignment[key] = _normalize_format_value(value)

        sample_assignments[str(sample_name)] = assignment
        if is_called:
            called_samples.append(str(sample_name))
        if is_nonref:
            nonref_samples.append(str(sample_name))
        if has_alt_support:
            alt_support_samples.append(str(sample_name))

    return {
        "sample_assignments": sample_assignments,
        "called_samples": called_samples,
        "nonref_samples": nonref_samples,
        "alt_support_samples": alt_support_samples,
    }


def read_vcf_to_df(vcf_file, kanpig_read_names=None, include_sample_assignments=False):
    """Parse a structural-variant VCF into a dataframe.

    Args:
        vcf_file: Path to an input VCF/BCF readable by pysam.
        kanpig_read_names: Optional two-column kanpig RNAMES TSV (`sv_id`, `read_name`).
        include_sample_assignments: If True, also parse per-sample FORMAT fields
            into `sample_assignments`, `called_samples`, `nonref_samples`, and
            `alt_support_samples`. This is useful for multi-sample kanpig output
            VCFs and is disabled by default so existing callers are unaffected.

    Returns:
        pandas.DataFrame
    """
    records = []
    vcf_file = pysam.VariantFile(vcf_file)
    for record in vcf_file.fetch():
        svtype = str(_safe_info_get(record.info, "SVTYPE", "")).upper()
        if svtype in ["INS", "DEL"]:
            sv_len = _safe_info_get(record.info, "SVLEN", "NA")
            stdev_len = _safe_info_get(record.info, "STDEV_LEN", "NA")
            stdev_pos = _safe_info_get(record.info, "STDEV_POS", "NA")
            end = getattr(record, "stop", pd.NA)
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
                "sv_type": svtype,
                "sv_len": sv_len,
                "supporting_reads": _safe_info_get_read_names(record.info, "RNAMES"),
                "stdev_len": stdev_len,
                "stdev_pos": stdev_pos,
                "end": end,
                "vaf": vaf,
            }

            if include_sample_assignments:
                df_record.update(_extract_sample_assignments(record))

            stdev_pos = df_record["stdev_pos"]
            sv_len = df_record["sv_len"] if svtype == "DEL" else 0
            stdev_len = df_record["stdev_len"] 
            end = _safe_int(df_record["end"])
            if stdev_pos != "NA":
                ref_start = df_record["location"] - stdev_pos
                ref_end = df_record["location"] + stdev_pos + stdev_len - sv_len
            elif not np.isnan(end):
                ref_start = _safe_int(df_record["location"])
                ref_end = max(_safe_int(df_record["location"]), end)
            else:
                ref_start = _safe_int(df_record["location"])
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
    columns = [
        "chr",
        "location",
        "id",
        "sv_type",
        "sv_len",
        "supporting_reads",
        "stdev_len",
        "stdev_pos",
        "ref_start",
        "ref_end",
        "vaf",
    ]
    if include_sample_assignments:
        columns.extend(
            [
                "sample_assignments",
                "called_samples",
                "nonref_samples",
                "alt_support_samples",
            ]
        )
    sv_df = pd.DataFrame(records, columns=columns)
    return sv_df


def read_kanpig_vcf_to_df(vcf_file, kanpig_read_names=None):
    """Parse a kanpig multi-sample VCF into a dataframe with sample assignments.

    This is a convenience wrapper around `read_vcf_to_df(..., include_sample_assignments=True)`.
    The returned dataframe includes:
      - `sample_assignments`: per-sample FORMAT payloads keyed by sample name
      - `called_samples`: samples with a non-missing GT
      - `nonref_samples`: samples with a non-reference GT
      - `alt_support_samples`: samples with non-zero alternate depth
    """
    return read_vcf_to_df(
        vcf_file,
        kanpig_read_names=kanpig_read_names,
        include_sample_assignments=True,
    )
