import numpy as np
import pandas as pd

def filter_bed_based_on_variants(bed_df: pd.DataFrame, sv_df: pd.DataFrame, window: int = 5000) -> pd.DataFrame:
    # 1) Normalize column names
    bed = bed_df.copy() 
    sv  = sv_df.rename(columns={'chr': 'chr'}).copy()
    if 'supporting_reads' not in sv.columns:
        exit("VCF DataFrame missing 'supporting_reads' column required for SV assignment.")
    # 2) Normalize chromosome naming (strip or add 'chr' to match)
    def _norm_chr(x):
        x = str(x)
        return x[3:] if x.startswith('chr') else x
    bed['chr'] = bed['chr'].map(_norm_chr)
    sv['chr']  = sv['chr'].map(_norm_chr)

    # 3) Enforce integer dtype
    bed[['start','end']] = bed[['start','end']].astype(np.int64)
    sv[['ref_start','ref_end']] = sv[['ref_start','ref_end']].astype(np.int64)

    # 4) Convert VCF (1-based inclusive) -> BED-style half-open
    #    [ref_start-1, ref_end) in 0-based half-open. Equivalently:
    #    start0 = ref_start - 1; end0_exclusive = ref_end
    sv_start0 = (sv['ref_start'].to_numpy(np.int64) - 1)
    sv_end0   = sv['ref_end'].to_numpy(np.int64)      # already exclusive after conversion

    # Handle insertions robustly: if END < POS (some callers), clamp end to POS.
    bad = sv_end0 < (sv_start0 + 1)
    if np.any(bad):
        sv_end0[bad] = sv_start0[bad] + 1

    sv = sv.assign(_start=sv_start0, _end=sv_end0)

    # 5) Optional: make chr categorical for faster groupbys
    bed['chr'] = bed['chr'].astype('category')
    sv['chr']  = sv['chr'].astype('category')

    out_mask = np.zeros(len(bed), dtype=bool)

    # Build per-chrom sorted arrays of SV starts/ends (already sorted? we still ensure monotonic)
    sv_idx = {}
    for chrom, sdf in sv.groupby('chr', sort=False, observed=False):
        starts = sdf['_start'].to_numpy(np.int64)
        ends   = sdf['_end'].to_numpy(np.int64)
        # Ensure sorted (cheap if already sorted)
        order = np.argsort(starts, kind='mergesort')
        starts = starts[order]
        ends   = ends[order]  # maintain pairing order
        # We also want an array of ends sorted to use searchsorted independently.
        ends_sorted = np.sort(ends, kind='mergesort')
        sv_idx[chrom] = (starts, ends, ends_sorted)


    for chrom, bdf in bed.groupby('chr', sort=False, observed=False):
        if chrom not in sv_idx:
            continue
        starts, _, ends_sorted = sv_idx[chrom]
        if starts.size == 0:
            continue

        # Core (original) BED interval
        core_start = bdf['start'].to_numpy(np.int64)
        core_end   = bdf['end'].to_numpy(np.int64)

        # Padded interval
        bed_starts = core_start - window
        bed_ends   = core_end   + window

        # ---------- 1) Overlap with *padded* interval (half-open, as before) ----------
        n_start_lt_end_pad = np.searchsorted(starts,      bed_ends,   side='left')
        n_end_le_start_pad = np.searchsorted(ends_sorted, bed_starts, side='right')
        overlap_padded = (n_start_lt_end_pad - n_end_le_start_pad) > 0

        # ---------- 2) Overlap with *core* interval (treat breakpoint as overlap) ------
        # Closed-interval style: sv_start <= core_end AND sv_end >= core_start
        n_start_le_end_core = np.searchsorted(starts,      core_end,   side='right')
        n_end_lt_start_core = np.searchsorted(ends_sorted, core_start, side='left')
        overlap_core = (n_start_le_end_core - n_end_lt_start_core) > 0

        # ---------- 3) We want SV in padding, but NOT in core --------------------------
        keep = overlap_padded & (~overlap_core)
        out_mask[bdf.index] = keep


    return bed_df.loc[out_mask]
