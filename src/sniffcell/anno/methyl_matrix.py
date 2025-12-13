import pandas as pd
import pysam, re
import numpy as np
from typing import Optional, List, Tuple, Union
from scipy import sparse
from math import log, exp

def _cpg_c_sites(fa: pysam.FastaFile, chrom: str, start: int, end: int):
    return [m.start() + start for m in re.finditer(r"CG", fa.fetch(chrom, start, end))]

def _combine_m_h(pm: float, ph: float, mode: str = "union") -> float:
    """
    Combine 5mC (pm) and 5hmC (ph) probabilities into one methylated probability.
    NaNs are treated as absence for that channel.
    """
    if (pm is None or np.isnan(pm)) and (ph is None or np.isnan(ph)):
        return np.nan
    pm = 0.0 if (pm is None or np.isnan(pm)) else pm
    ph = 0.0 if (ph is None or np.isnan(ph)) else ph

    if mode == "union":              # p(m or h)
        return 1.0 - (1.0 - pm) * (1.0 - ph)
    elif mode == "max":              # max(pm, ph)
        return pm if pm >= ph else ph
    elif mode == "mean":             # (pm + ph)/2
        return 0.5 * (pm + ph)
    elif mode == "logit_sum":        # invlogit(logit(pm) + logit(ph))
        eps = 1e-6
        def logit(p):
            p = min(max(p, eps), 1 - eps)
            return log(p / (1 - p))
        l = logit(pm) + logit(ph)
        return 1.0 / (1.0 + exp(-l))
    else:
        # fallback: max
        return pm if pm >= ph else ph

def methyl_matrix_from_bam(
    bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
    min_read_length: int = 0, include_secondary: bool = False,
    include_supplementary: bool = False, include_unmapped: bool = False,
    as_sparse: bool = False, return_positions: bool = False,
    wanted_keys: Optional[set] = None,
    combine_mode: str = "union",   # <-- NEW: how to combine m & h
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, List[int]]]:

    # default: include both 5mC and 5hmC on both strands
    wanted_keys = wanted_keys or {
        ('C', 0, 'm'), ('C', 1, 'm'),
        ('C', 0, 'h'), ('C', 1, 'h'),
    }

    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(fasta_path) as fa:
        cpgs = _cpg_c_sites(fa, chrom, start, end)
        if not cpgs:
            idx = pd.MultiIndex.from_arrays([[], []], names=["read_name", "haplotype"])
            out = pd.DataFrame(index=idx)
            return (out, cpgs) if return_positions else out
        col_index = {p: j for j, p in enumerate(cpgs)}

        # We’ll accumulate per-cell (rid, j) the best pm and ph seen
        # cell_probs[(rid, j)] = [pm, ph]
        cell_probs = {}

        row_ids, haps, key2row = [], [], {}

        for r in bam.fetch(chrom, start, end, multiple_iterators=True):
            if (r.is_unmapped and not include_unmapped) or \
               (r.is_secondary and not include_secondary) or \
               (r.is_supplementary and not include_supplementary) or \
               (min_read_length and (r.query_length or 0) < min_read_length):
                continue

            try:
                hp = r.get_tag("HP")
            except (KeyError, AttributeError):
                hp = -1

            mb = getattr(r, "modified_bases", None)
            if not mb:
                continue
            refpos = r.get_reference_positions(full_length=True)
            if refpos is None:
                continue

            rk = (r.query_name, hp)
            rid = key2row.setdefault(rk, len(row_ids))
            if rid == len(row_ids):
                row_ids.append(r.query_name)
                haps.append(hp)

            for k, mods in mb.items():
                if k not in wanted_keys:
                    continue
                is_m = (k[2] == 'm')
                is_h = (k[2] == 'h')
                for qpos, score in mods:
                    if 0 <= qpos < len(refpos):
                        p = refpos[qpos]
                        if p is None:
                            continue
                        p_c = p - 1 if r.is_reverse else p
                        j = col_index.get(p_c)
                        if j is None:
                            continue
                        val = score / 255.0
                        pm, ph = cell_probs.get((rid, j), [np.nan, np.nan])
                        if is_m:
                            pm = val if (np.isnan(pm) or val > pm) else pm
                        elif is_h:
                            ph = val if (np.isnan(ph) or val > ph) else ph
                        cell_probs[(rid, j)] = [pm, ph]

        idx = pd.MultiIndex.from_arrays([row_ids, haps], names=["read_name", "haplotype"])
        if not row_ids:
            out = pd.DataFrame(columns=cpgs, index=idx)
            return (out, cpgs) if return_positions else out

        # Build COO from combined values
        if cell_probs:
            data_i, data_j, data_v = [], [], []
            for (i, j), (pm, ph) in cell_probs.items():
                v = _combine_m_h(pm, ph, combine_mode)
                if not np.isnan(v):
                    data_i.append(i); data_j.append(j); data_v.append(v)
            coo = sparse.coo_matrix((data_v, (data_i, data_j)), shape=(len(row_ids), len(cpgs)))
        else:
            coo = sparse.coo_matrix((len(row_ids), len(cpgs)))

        if as_sparse:
            df = pd.DataFrame.sparse.from_spmatrix(coo, index=idx, columns=cpgs)
        else:
            arr = np.full((len(row_ids), len(cpgs)), np.nan, float)
            if coo.nnz:
                arr[coo.row, coo.col] = coo.data
            df = pd.DataFrame(arr, index=idx, columns=cpgs)

        return (df, cpgs) if return_positions else df
