from __future__ import annotations

import numpy as np
import pandas as pd


def validate_breakpoint_exclusion_frac(value: float) -> float:
    frac = float(value)
    if frac < 0:
        raise ValueError("breakpoint_exclusion_frac must be >= 0")
    return frac


def compute_breakpoint_exclusion_bp(
    sv_df: pd.DataFrame,
    *,
    breakpoint_exclusion_frac: float,
    sv_len_col: str = "sv_len",
    start_col: str = "_start",
    end_col: str = "_end",
) -> np.ndarray:
    frac = validate_breakpoint_exclusion_frac(breakpoint_exclusion_frac)
    if sv_df.empty or frac == 0.0:
        return np.zeros(len(sv_df), dtype=np.int64)

    if sv_len_col in sv_df.columns:
        sv_len = pd.to_numeric(sv_df[sv_len_col], errors="coerce").abs()
    else:
        sv_len = pd.Series(np.nan, index=sv_df.index, dtype="float64")

    start = pd.to_numeric(sv_df[start_col], errors="coerce") if start_col in sv_df.columns else pd.Series(np.nan, index=sv_df.index, dtype="float64")
    end = pd.to_numeric(sv_df[end_col], errors="coerce") if end_col in sv_df.columns else pd.Series(np.nan, index=sv_df.index, dtype="float64")
    span = (end - start).abs()

    effective_len = sv_len.where(sv_len.notna(), span)
    effective_len = effective_len.fillna(0.0).clip(lower=0.0)
    return np.ceil(effective_len.to_numpy(dtype=float) * frac).astype(np.int64)
