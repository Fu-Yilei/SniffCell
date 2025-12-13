from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam
import pandas as pd
from typing import Tuple, Iterable
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm


def statistical_test_on_methylation(methyl_matrix1, methyl_matrix2, test_type="t-test"):
    """
    Perform statistical test on two methylation matrices.

    Args:
        methyl_matrix1 (pd.DataFrame): Methylation matrix for group 1.
        methyl_matrix2 (pd.DataFrame): Methylation matrix for group 2.
        test_type (str): Type of statistical test to perform. Options are "t-test", "mannwhitneyu", "fisher".

    Returns:
        dict: A dictionary with CpG positions as keys and p-values as values.
    """
    from scipy.stats import ttest_ind, mannwhitneyu, fisher_exact
    import numpy as np

    p_values = {}
    common_cpgs = set(methyl_matrix1.columns).intersection(set(methyl_matrix2.columns))

    for cpg in common_cpgs:
        group1_values = methyl_matrix1[cpg].dropna(how='all').values
        group2_values = methyl_matrix2[cpg].dropna(how='all').values
        # print(group1_values, group2_values)
        if len(group1_values) < 2 or len(group2_values) < 2:
            p_values[cpg] = np.nan
            continue

        if test_type == "t-test":
            stat, p = ttest_ind(group1_values, group2_values, equal_var=False)
        elif test_type == "mannwhitneyu":
            stat, p = mannwhitneyu(group1_values, group2_values, alternative='two-sided')
        elif test_type == "fisher":
            # Create contingency table
            table = [
                [np.sum(group1_values >= 0.5), np.sum(group1_values < 0.5)],
                [np.sum(group2_values >= 0.5), np.sum(group2_values < 0.5)]
            ]
            stat, p = fisher_exact(table)
        else:
            raise ValueError(f"Unsupported test type: {test_type}")

        p_values[cpg] = p

    return p_values


def _split_supporting_vs_other(
    methyl_df: pd.DataFrame,
    supporting_reads: Iterable[str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Given a methylation matrix with MultiIndex rows (read_name, haplotype),
    return (supporting_reads_df, other_reads_df).
    """
    read_names = methyl_df.index.get_level_values("read_name")
    sup_set = set(str(r) for r in supporting_reads)
    mask = read_names.isin(sup_set)
    return methyl_df[mask], methyl_df[~mask]

def _drop_cpgs_overlapping_sv(
    methyl_matrix: pd.DataFrame,
    sv_start: int,
    sv_end: int,
    nan_threshold: float = 0.4,
    read_nan_threshold: float = 0.2,
) -> pd.DataFrame:
    """
    Remove CpG columns that:
      1. Fall inside the SV span [sv_start, sv_end]
      2. Have > nan_threshold fraction of NaN values (per column)
    Also remove reads (rows) with > read_nan_threshold fraction of NaN values.

    Assumes DataFrame columns are integer CpG positions (as emitted by methyl_matrix_from_bam).
    """
    cols = pd.Index(methyl_matrix.columns)

    # Coerce to integer CpG positions if needed
    if cols.dtype.kind not in {"i", "u"}:
        try:
            cols_int = pd.to_numeric(cols, errors="raise")
        except Exception:
            return methyl_matrix  # skip filtering if columns aren’t numeric
    else:
        cols_int = cols

    # (1) Drop CpGs overlapping the SV region
    keep_mask = (cols_int < sv_start) | (cols_int > sv_end)
    filtered = methyl_matrix.loc[:, keep_mask]

    # (3) Drop reads with too many NaNs (per read)
    nan_frac_rows = filtered.isna().mean(axis=1)
    keep_nan_mask_rows = nan_frac_rows <= read_nan_threshold
    filtered = filtered.loc[keep_nan_mask_rows, :]

    # (2) Drop CpGs with too many NaNs (per CpG)
    nan_frac_cols = filtered.isna().mean()
    keep_nan_mask_cols = nan_frac_cols <= nan_threshold
    filtered = filtered.loc[:, keep_nan_mask_cols]
    return filtered



def get_statistical_tests_around_sv(
    sv_methylation_data,
    haplotype_majority_threshold: float = 0.7,
    test_type: str = "t-test"
):
    """
    For each SV's methylation data, perform statistical tests between supporting reads and other reads.

    Haplotype handling:
      - Use all supporting reads to determine the major haplotype (1 or 2).
      - Ignore unphased reads (-1) when counting for haplotype identity, 
        but they still contribute to total read count for n_phased (denominator).
      - If the major haplotype constitutes >= threshold of all supporting reads,
        restrict both supporting and other reads to that haplotype.
      - Otherwise, use all reads regardless of haplotype.

    Args:
        sv_methylation_data (dict): Output from get_methylation_around_sv.
        haplotype_majority_threshold (float): Fraction threshold for haplotype dominance.
        test_type (str): Statistical test type ("t-test", "mannwhitneyu", "fisher").

    Returns:
        dict: {sv_id: {cpg_position: p_value}}
    """
    import logging
    import pandas as pd

    results = {}

    for sv_id, methyl_data in sv_methylation_data.items():
        mm_sup = methyl_data["methyl_matrix_supporting"]
        mm_oth = methyl_data["methyl_matrix_other"]

        sup_use = mm_sup
        oth_use = mm_oth

        # Check haplotype structure
        can_hap_filter = (
            isinstance(mm_sup.index, pd.MultiIndex)
            and "haplotype" in (mm_sup.index.names or [])
            and isinstance(mm_oth.index, pd.MultiIndex)
            and "haplotype" in (mm_oth.index.names or [])
        )

        if can_hap_filter and len(mm_sup) > 0:
            hap_sup = mm_sup.index.get_level_values("haplotype")
            # Normalize to int
            try:
                hap_sup = hap_sup.astype(int)
            except Exception:
                hap_sup = hap_sup.map(lambda x: {"1": 1, "2": 2, "-1": -1}.get(str(x), -1))

            # Count all reads
            n_total = len(hap_sup)
            n_hp1 = int((hap_sup == 1).sum())
            n_hp2 = int((hap_sup == 2).sum())

            # Identify major haplotype based on all reads
            if n_hp1 >= n_hp2:
                major_hap = 1
                frac_major = n_hp1 / n_total
            else:
                major_hap = 2
                frac_major = n_hp2 / n_total

            if frac_major >= haplotype_majority_threshold:
                sup_use = mm_sup.xs(major_hap, level="haplotype", drop_level=True)
                try:
                    oth_use = mm_oth.xs(major_hap, level="haplotype", drop_level=True)
                except KeyError:
                    oth_use = pd.DataFrame(columns=mm_oth.columns)
                logging.info(
                    f"SV {sv_id}: using major haplotype {major_hap} "
                    f"({frac_major:.2f} of all supporting reads)."
                )
            else:
                logging.warning(
                    f"SV {sv_id}: no clear haplotype majority ({frac_major:.2f}); "
                    f"using all reads."
                )
        else:
            if not can_hap_filter:
                logging.info(f"SV {sv_id}: no haplotype information; using all reads.")
        # Perform the statistical test per CpG
        p_values = statistical_test_on_methylation(
            methyl_matrix1=sup_use,
            methyl_matrix2=oth_use,
            test_type=test_type
        )
        results[sv_id] = p_values

    return results

def _worker_sv(rec, input_bam, reference_fasta, flank_size, min_supporting_reads):
    chrom = rec['chr']
    start = max(0, rec['ref_start'] - flank_size)
    end = rec['ref_end'] + flank_size
    sv_id = rec['id']
    region_id = f"{chrom}:{start}-{end}"
    if len(rec['supporting_reads']) < min_supporting_reads:
        return None
    methyl_matrix = methyl_matrix_from_bam(
        bam_path=input_bam,
        fasta_path=reference_fasta,
        chrom=chrom,
        start=start,
        end=end
    )
    methyl_matrix = _drop_cpgs_overlapping_sv(
        methyl_matrix,
        sv_start=rec['ref_start'],
        sv_end=rec['ref_end'],
    )
    sv_supporting_reads = rec['supporting_reads']
    # print(methyl_matrix)
    methyl_matrix_supporting_reads, methyl_matrix_other_reads = _split_supporting_vs_other(
        methyl_df=methyl_matrix,
        supporting_reads=sv_supporting_reads
    )
    # print(methyl_matrix_supporting_reads, methyl_matrix_other_reads)

    return sv_id, {
        'region_id': region_id,
        'methyl_matrix_supporting': methyl_matrix_supporting_reads,
        'methyl_matrix_other': methyl_matrix_other_reads
    }

def get_methylation_around_sv(
    input_bam,
    reference_fasta,
    sv_df,
    min_supporting_reads=3,
    flank_size=5000,
    n_threads=4   # <---- specify your thread (process) count here
):
    """
    For each SV in sv_df, extract methylation matrix around the SV region
    (flanked by flank_size on both sides) from the input BAM file.

    Args:
        input_bam (str): Path to the input BAM file.
        reference_fasta (str): Path to the reference FASTA file.
        sv_df (pd.DataFrame): DataFrame containing SV information with columns
                              ['chr', 'ref_start', 'ref_end'].
        flank_size (int): Number of base pairs to flank on both sides of the SV.
    Returns:
        dict: A dictionary where keys are SV identifiers and values are methylation matrices.
    """
    sv_methylation_data = {}
    records = sv_df.to_dict('records')

    with ProcessPoolExecutor(max_workers=n_threads) as ex:
        futures = [
            ex.submit(
                _worker_sv, rec, input_bam, reference_fasta, flank_size, min_supporting_reads
            )
            for rec in records
        ]
        # print(futures)
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Processing SVs"):
            res = fut.result()
            if res is None:
                continue
            sv_id, payload = res
            sv_methylation_data[sv_id] = payload

    return sv_methylation_data