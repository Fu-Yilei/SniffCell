import pandas as pd
import pysam, re
import numpy as np
from typing import Optional, List, Tuple, Union
from scipy import sparse

def _cpg_c_sites(fa: pysam.FastaFile, chrom: str, start: int, end: int):
    return [m.start() + start for m in re.finditer(r"CG", fa.fetch(chrom, start, end))]

def methyl_matrix_from_bam(
    bam_path: str, fasta_path: str, chrom: str, start: int, end: int,
    min_read_length: int = 0, include_secondary: bool = False,
    include_supplementary: bool = False, include_unmapped: bool = False,
    as_sparse: bool = False, return_positions: bool = False,
    wanted_keys: Optional[set] = None,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, List[int]]]:

    wanted_keys = wanted_keys or {('C', 0, 'm'), ('C', 1, 'm')}

    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(fasta_path) as fa:
        cpgs = _cpg_c_sites(fa, chrom, start, end)
        if not cpgs:
            idx = pd.MultiIndex.from_arrays([[], []], names=["read_name", "haplotype"])
            out = pd.DataFrame(index=idx)
            return (out, cpgs) if return_positions else out
        col_index = {p: j for j, p in enumerate(cpgs)}

        data_i, data_j, data_v = [], [], []
        row_ids, haps, key2row = [], [], {}

        for r in bam.fetch(chrom, start, end, multiple_iterators=True):
            if (r.is_unmapped and not include_unmapped) or \
               (r.is_secondary and not include_secondary) or \
               (r.is_supplementary and not include_supplementary) or \
               (min_read_length and (r.query_length or 0) < min_read_length):
                continue

            # --- Safe HP tag handling ---
            try:
                hp = r.get_tag("HP")
            except (KeyError, AttributeError):
                hp = -1
            # ----------------------------

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
                for qpos, score in mods:
                    if 0 <= qpos < len(refpos):
                        p = refpos[qpos]
                        if p is None:
                            continue
                        p_c = p - 1 if r.is_reverse else p
                        j = col_index.get(p_c)
                        if j is not None:
                            data_i.append(rid)
                            data_j.append(j)
                            data_v.append(score / 255.0)

        idx = pd.MultiIndex.from_arrays([row_ids, haps], names=["read_name", "haplotype"])
        if not row_ids:
            out = pd.DataFrame(columns=cpgs, index=idx)
            return (out, cpgs) if return_positions else out

        coo = sparse.coo_matrix((data_v, (data_i, data_j)), shape=(len(row_ids), len(cpgs)))
        if as_sparse:
            df = pd.DataFrame.sparse.from_spmatrix(coo, index=idx, columns=cpgs)
        else:
            arr = np.full((len(row_ids), len(cpgs)), np.nan, float)
            if coo.nnz:
                arr[coo.row, coo.col] = coo.data
            df = pd.DataFrame(arr, index=idx, columns=cpgs)

        return (df, cpgs) if return_positions else df