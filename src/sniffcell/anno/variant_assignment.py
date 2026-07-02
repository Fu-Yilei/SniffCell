import numpy as np
import pandas as pd
from sniffcell.anno.breakpoint_exclusion import (
    compute_breakpoint_exclusion_bp,
    validate_breakpoint_exclusion_frac,
)

def _split_pipe_values(value) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [x.strip() for x in text.split("|") if x.strip()]


def _parse_support_read_names(value) -> list[str]:
    """
    Parse SV supporting read payloads from VCF-derived values.
    Accepts list-like objects or delimiter-separated strings.
    """
    null_tokens = {"", "NA", "NAN", "NONE", "."}

    if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
        out: list[str] = []
        for item in value:
            text = str(item).strip()
            if text and text.upper() not in null_tokens:
                out.append(text)
        return out

    if isinstance(value, str):
        text = value.strip()
        if not text or text.upper() in null_tokens:
            return []
        for delim in (",", "|", ";"):
            if delim in text:
                return [t.strip() for t in text.split(delim) if t.strip() and t.strip().upper() not in null_tokens]
        return [text]

    try:
        if pd.isna(value):
            return []
    except TypeError:
        pass

    text = str(value).strip()
    if not text or text.upper() in null_tokens:
        return []
    return [text]


def _normalize_binary_code(code_value, schema_value) -> object:
    """
    Normalize binary code strings and preserve leading zeros based on schema length.
    Example: schema A|B|C|D with code 111 -> 0111.
    """
    if pd.isna(code_value):
        return pd.NA

    if isinstance(code_value, (int, float)) and not isinstance(code_value, bool):
        code_num = float(code_value)
        if code_num.is_integer():
            code = str(int(code_num))
        else:
            code = str(code_value).strip()
    else:
        code = str(code_value).strip()

    if not code:
        return ""

    if code.endswith(".0") and code[:-2].isdigit():
        code = code[:-2]

    if pd.notna(schema_value):
        labels = [ct.strip() for ct in str(schema_value).split("|") if ct.strip()]
        if labels and set(code).issubset({"0", "1"}) and len(code) < len(labels):
            code = code.zfill(len(labels))

    return code


def _decode_linked_celltypes_from_row(row: pd.Series) -> tuple[str, list[str]]:
    """
    Decode linked cell types from row-level code payload.
    Priority:
      1) explicit code_order + code columns
      2) code_token payload (schema::bits or bits only)
      3) fallback to best_group if is_best_group=True
    """
    code_order = row.get("code_order", pd.NA)
    code = row.get("code", pd.NA)
    code_token = row.get("code_token", pd.NA)
    best_group = row.get("best_group", pd.NA)
    best_group_leaves = row.get("best_group_leaves", pd.NA)
    is_best_group = bool(row.get("is_best_group", False))

    schema = ""
    bits = ""

    if pd.notna(code_order):
        schema = str(code_order).strip()
    if pd.notna(code):
        bits = str(code).strip()

    if pd.notna(code_token):
        token = str(code_token).strip()
        if "::" in token:
            token_schema, token_bits = token.split("::", 1)
            if not schema:
                schema = token_schema.strip()
            if not bits:
                bits = token_bits.strip()
        elif not bits:
            bits = token

    linked: list[str] = []
    if schema and bits and set(bits).issubset({"0", "1"}):
        labels = [ct.strip() for ct in schema.split("|") if ct.strip()]
        if len(bits) < len(labels):
            bits = bits.zfill(len(labels))
        if len(labels) == len(bits):
            linked = [ct for ct, bit in zip(labels, bits) if bit == "1"]

    if (not linked) and is_best_group and pd.notna(best_group):
        bg_leaves = _split_pipe_values(best_group_leaves)
        if bg_leaves:
            linked = bg_leaves
        else:
            bg = str(best_group).strip()
            if bg:
                linked = [bg]
                if not schema:
                    schema = bg

    return schema, linked


def _build_group_leaf_sets(evidence: pd.DataFrame) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = {}

    for group_col, leaves_col in [("best_group", "best_group_leaves"), ("other_group", "other_group_leaves")]:
        if group_col not in evidence.columns or leaves_col not in evidence.columns:
            continue
        # Drop duplicates first — for 11M-row DataFrames this reduces iterations from millions to dozens.
        unique_pairs = (
            evidence[[group_col, leaves_col]]
            .drop_duplicates()
            .dropna(subset=[group_col])
        )
        for _, row in unique_pairs.iterrows():
            group = str(row[group_col]).strip()
            if not group:
                continue
            leaves = _split_pipe_values(row[leaves_col])
            if not leaves:
                continue
            leaves_set = set(leaves)
            # Keep the smallest declared set as the most specific mapping for this label.
            if group not in mapping or len(leaves_set) < len(mapping[group]):
                mapping[group] = leaves_set

    if "best_group" in evidence.columns:
        for group in evidence["best_group"].dropna().astype(str).unique():
            group = group.strip()
            if group and group not in mapping:
                mapping[group] = {group}

    return mapping


def _resolve_hierarchy_labels(linked_leaves: list[str], group_leaf_sets: dict[str, set[str]]) -> list[str]:
    leaf_set = set(linked_leaves)
    if not leaf_set:
        return []

    candidates: list[str] = []
    for group, leaves in group_leaf_sets.items():
        if leaf_set == leaves:
            candidates.append(group)

    if not candidates:
        out: list[str] = []
        seen = set()
        for ct in linked_leaves:
            if ct not in seen:
                seen.add(ct)
                out.append(ct)
        return out

    return sorted(candidates)


def _norm_chr(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    return text[3:] if text.lower().startswith("chr") else text


def _split_code_token_schema_bits(code_token: object) -> tuple[str, str]:
    text = "" if pd.isna(code_token) else str(code_token).strip()
    if not text:
        return "", ""
    if "::" in text:
        schema, bits = text.split("::", 1)
        schema = schema.strip()
        bits = bits.strip()
    else:
        schema = ""
        bits = text
    if bits.endswith(".0") and bits[:-2].isdigit():
        bits = bits[:-2]
    return schema, bits


def _bitwise_intersection(bits_list: list[str]) -> str:
    if not bits_list:
        return ""
    width = max(len(b) for b in bits_list)
    normalized = [b.zfill(width) for b in bits_list]
    out = list(normalized[0])
    for bits in normalized[1:]:
        out = ["1" if (a == "1" and b == "1") else "0" for a, b in zip(out, bits)]
    return "".join(out)


def _summarize_code_conflicts(read_votes: pd.DataFrame, sv_id_col: str) -> pd.DataFrame:
    cols = ["has_hard_conflict", "intersection_code"]
    if read_votes.empty or "code_token" not in read_votes.columns:
        return pd.DataFrame(columns=cols).set_index(pd.Index([], name=sv_id_col))

    rows = []
    for sv_id, g in read_votes.groupby(sv_id_col, sort=False):
        schema_to_bits: dict[str, list[str]] = {}
        schema_vote_count: dict[str, int] = {}
        for token in g["code_token"].dropna():
            schema, bits = _split_code_token_schema_bits(token)
            if not bits or not set(bits).issubset({"0", "1"}):
                continue
            schema_to_bits.setdefault(schema, []).append(bits)
            schema_vote_count[schema] = schema_vote_count.get(schema, 0) + 1

        has_hard_conflict = False
        intersection_code = pd.NA

        if schema_vote_count:
            # Use the most represented schema as the primary intersection view.
            top_schema = sorted(schema_vote_count.items(), key=lambda x: (-x[1], x[0]))[0][0]
            top_bits = schema_to_bits.get(top_schema, [])
            if top_bits:
                inter = _bitwise_intersection(top_bits)
                intersection_code = f"{top_schema}::{inter}" if top_schema else inter

            # Hard conflict: within any schema, constraints intersect to empty set.
            for bits_list in schema_to_bits.values():
                if not bits_list:
                    continue
                inter = _bitwise_intersection(bits_list)
                unique_codes = len(set(bits_list))
                if inter and ("1" not in inter) and unique_codes > 1:
                    has_hard_conflict = True
                    break

        rows.append(
            {
                sv_id_col: sv_id,
                "has_hard_conflict": has_hard_conflict,
                "intersection_code": intersection_code,
            }
        )

    return pd.DataFrame(rows).set_index(sv_id_col)


def _build_sv_region_map(
    region_df: pd.DataFrame,
    sv_regions: pd.DataFrame,
    *,
    sv_id_col: str,
    window: int,
) -> pd.DataFrame:
    out_cols = [sv_id_col, "chr_norm", "_start", "_end"]
    if region_df.empty or sv_regions.empty:
        return pd.DataFrame(columns=out_cols)

    links: list[pd.DataFrame] = []
    for chrom, cdf in region_df.groupby("chr_norm", sort=False, observed=False):
        sv_chrom = sv_regions[sv_regions["chr_norm"] == chrom]
        if sv_chrom.empty:
            continue

        starts = cdf["_start"].to_numpy(np.int64)
        ends = cdf["_end"].to_numpy(np.int64)

        for _, sv_row in sv_chrom.iterrows():
            sv_start = int(sv_row["_start"])
            sv_end = int(sv_row["_end"])

            # Keep DMRs within padded distance from SV, but not directly overlapping the core DMR.
            overlap_padded = (sv_start < (ends + window)) & (sv_end > (starts - window))
            overlap_core = (sv_start <= ends) & (sv_end >= starts)
            keep = overlap_padded & (~overlap_core)
            if not np.any(keep):
                continue

            matched = cdf.iloc[np.flatnonzero(keep)][["chr_norm", "_start", "_end"]].copy()
            matched[sv_id_col] = sv_row[sv_id_col]
            links.append(matched[[sv_id_col, "chr_norm", "_start", "_end"]])

    if not links:
        return pd.DataFrame(columns=out_cols)

    return pd.concat(links, ignore_index=True).drop_duplicates(ignore_index=True)


def _build_sv_regions(
    sv_df: pd.DataFrame,
    *,
    sv_id_col: str,
    breakpoint_exclusion_frac: float,
) -> pd.DataFrame:
    sv_region_cols = [sv_id_col, "chr", "location"]
    for optional_col in ("ref_start", "ref_end", "sv_len"):
        if optional_col in sv_df.columns:
            sv_region_cols.append(optional_col)

    sv_regions = sv_df[sv_region_cols].drop_duplicates(subset=[sv_id_col]).copy()
    sv_regions["chr_norm"] = sv_regions["chr"].map(_norm_chr)
    sv_loc = pd.to_numeric(sv_regions["location"], errors="coerce")
    ref_start = (
        pd.to_numeric(sv_regions["ref_start"], errors="coerce")
        if "ref_start" in sv_regions.columns else pd.Series(np.nan, index=sv_regions.index, dtype="float64")
    )
    ref_end = (
        pd.to_numeric(sv_regions["ref_end"], errors="coerce")
        if "ref_end" in sv_regions.columns else pd.Series(np.nan, index=sv_regions.index, dtype="float64")
    )
    sv_regions["_start"] = np.where(ref_start.notna(), ref_start - 1, sv_loc - 1)
    sv_regions["_end"] = np.where(ref_end.notna(), ref_end, sv_loc)
    valid_sv = (
        sv_regions["chr_norm"].ne("")
        & pd.notna(sv_regions["_start"])
        & pd.notna(sv_regions["_end"])
    )
    sv_regions = sv_regions.loc[valid_sv].copy()
    if sv_regions.empty:
        return pd.DataFrame(columns=[sv_id_col, "chr_norm", "_start", "_end", "_exclude_bp"])

    sv_regions[["_start", "_end"]] = sv_regions[["_start", "_end"]].astype("int64")
    bad = sv_regions["_end"] < (sv_regions["_start"] + 1)
    if bad.any():
        sv_regions.loc[bad, "_end"] = sv_regions.loc[bad, "_start"] + 1

    sv_regions["_exclude_bp"] = compute_breakpoint_exclusion_bp(
        sv_regions,
        breakpoint_exclusion_frac=breakpoint_exclusion_frac,
    )
    return sv_regions[[sv_id_col, "chr_norm", "_start", "_end", "_exclude_bp"]].copy()


def _filter_mapped_breakpoint_exclusion(
    mapped: pd.DataFrame,
    sv_regions: pd.DataFrame,
    *,
    sv_id_col: str,
) -> pd.DataFrame:
    if mapped.empty or sv_regions.empty:
        return mapped

    sv_lookup = (
        sv_regions[[sv_id_col, "_start", "_end", "_exclude_bp"]]
        .drop_duplicates(subset=[sv_id_col])
        .rename(columns={"_start": "_sv_start", "_end": "_sv_end"})
    )
    filtered = mapped.merge(sv_lookup, on=sv_id_col, how="left")

    dmr_start = pd.to_numeric(filtered["start"], errors="coerce")
    dmr_end = pd.to_numeric(filtered["end"], errors="coerce")
    sv_start = pd.to_numeric(filtered["_sv_start"], errors="coerce")
    sv_end = pd.to_numeric(filtered["_sv_end"], errors="coerce")
    exclude_bp = pd.to_numeric(filtered["_exclude_bp"], errors="coerce").fillna(0.0)

    overlap_excluded = (
        dmr_start.notna()
        & dmr_end.notna()
        & sv_start.notna()
        & sv_end.notna()
        & ((sv_start - exclude_bp) <= dmr_end)
        & ((sv_end + exclude_bp) >= dmr_start)
    )
    return filtered.loc[~overlap_excluded, mapped.columns].copy()


def _summarize_celltype_links(
    evidence: pd.DataFrame,
    sv_id_col: str,
    unique_reads_for_overlap: bool = True,
) -> pd.DataFrame:
    """
    Build per-SV readable cell-type summaries from row-level codes.
    Each read linked to an SV contributes one winning vote (unless unique_reads_for_overlap=False).
    """
    cols = [
        "primary_celltype",
        "majority_linked_celltypes",
        "majority_excluded_celltypes",
        "majority_schema",
        "linked_celltypes",
        "linked_celltype_counts",
        "linked_celltype_fractions",
        "is_multi_celltype_link",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=cols).set_index(pd.Index([], name=sv_id_col))

    group_leaf_sets = _build_group_leaf_sets(evidence)

    decoded_rows = []
    for _, row in evidence.iterrows():
        schema, linked = _decode_linked_celltypes_from_row(row)
        if not linked:
            continue
        resolved = _resolve_hierarchy_labels(linked, group_leaf_sets)
        if not resolved:
            continue
        decoded_rows.append(
            {
                sv_id_col: row[sv_id_col],
                "read": row["read"],
                "schema": schema,
                "linked": tuple(resolved),
                "linked_key": "|".join(resolved),
            }
        )
    if not decoded_rows:
        return pd.DataFrame(columns=cols).set_index(pd.Index([], name=sv_id_col))

    decoded = pd.DataFrame(decoded_rows)
    decoded["read"] = decoded["read"].astype(str)

    if unique_reads_for_overlap:
        per_read = (
            decoded.groupby([sv_id_col, "read", "linked_key"], sort=False)
               .size().rename("count").reset_index()
        )
        per_read = per_read.sort_values(["count", "linked_key"], ascending=[False, True], kind="stable")
        read_votes = (
            per_read.drop_duplicates(subset=[sv_id_col, "read"], keep="first")
            .merge(
                decoded[[sv_id_col, "read", "linked_key", "schema", "linked"]].drop_duplicates(),
                on=[sv_id_col, "read", "linked_key"],
                how="left",
            )
        )
    else:
        read_votes = decoded.copy()
        read_votes["count"] = 1

    if read_votes.empty:
        return pd.DataFrame(columns=cols).set_index(pd.Index([], name=sv_id_col))

    rows = []
    for sv_id, g in read_votes.groupby(sv_id_col, sort=False):
        cell_counts: dict[str, int] = {}
        n_votes = int(len(g))

        for _, row in g.iterrows():
            for celltype in row["linked"]:
                cell_counts[celltype] = cell_counts.get(celltype, 0) + 1

        if not cell_counts:
            continue

        ranked_cells = sorted(cell_counts.items(), key=lambda x: (-x[1], x[0]))
        linked = [ct for ct, _ in ranked_cells]
        counts = [cnt for _, cnt in ranked_cells]

        primary = linked[0] if linked else ""
        counts_str = ";".join(f"{ct}:{cnt}" for ct, cnt in zip(linked, counts))
        fracs_str = (
            ";".join(f"{ct}:{(cnt / n_votes):.3f}" for ct, cnt in zip(linked, counts))
            if n_votes > 0 else ""
        )

        rows.append(
            {
                sv_id_col: sv_id,
                "primary_celltype": primary,
                "majority_linked_celltypes": primary,
                "majority_excluded_celltypes": "",
                "majority_schema": "",
                "linked_celltypes": "|".join(linked),
                "linked_celltype_counts": counts_str,
                "linked_celltype_fractions": fracs_str,
                "is_multi_celltype_link": len(linked) > 1,
            }
        )

    return pd.DataFrame(rows).set_index(sv_id_col)

def _compute_per_read_consensus_df(
    assignment_df: pd.DataFrame,
    per_read_min_agreement: float = 0.66,
) -> pd.DataFrame:
    """
    Compute a single consensus cell-type code for each read from all its ctDMR votes.

    Uses bitwise intersection: codes "1101" and "0001" are NOT a conflict — their
    intersection "0001" is the consensus. A genuine conflict (intersection collapses
    to all zeros with >1 unique code) falls back to plurality majority; if plurality
    fraction >= per_read_min_agreement the plurality code wins, otherwise the read is
    marked mixed and excluded from SV-level voting.

    The threshold (per_read_min_agreement) is only active in the conflict fallback path.
    When intersection succeeds, all codes are mathematically compatible (agreement = 100%).

    Returns a DataFrame with one row per read and columns:
        readname, read_consensus_code, read_agreement_pct,
        read_is_mixed, read_n_ctdmrs, read_used_intersection
    """
    if "code_token" not in assignment_df.columns:
        raise ValueError("assignment_df must have a 'code_token' column")

    flat = assignment_df.reset_index()
    idx_name = assignment_df.index.name or "index"
    if "readname" not in flat.columns:
        flat = flat.rename(columns={idx_name: "readname"})
    flat["readname"] = flat["readname"].astype(str)

    n_ctdmrs = flat.groupby("readname", sort=False).size().rename("read_n_ctdmrs")

    split_result = flat["code_token"].map(_split_code_token_schema_bits)
    flat["_schema"] = [r[0] for r in split_result]
    flat["_bits"]   = [r[1] for r in split_result]

    valid = flat["_bits"].str.match(r"^[01]+$", na=False)
    parsed = flat.loc[valid, ["readname", "_schema", "_bits"]].copy()

    if parsed.empty:
        result = pd.DataFrame({
            "readname":            n_ctdmrs.index.astype(str),
            "read_consensus_code": None,
            "read_agreement_pct":  0.0,
            "read_is_mixed":       True,
            "read_used_intersection": False,
        })
        result["read_n_ctdmrs"] = result["readname"].map(n_ctdmrs).fillna(0).astype(int)
        return result

    # Group by (readname, schema) → list of bits, then pick dominant schema per read.
    schema_agg = (
        parsed.groupby(["readname", "_schema"], sort=False)["_bits"]
              .agg(list)
              .reset_index()
              .rename(columns={"_bits": "_bits_list"})
    )
    schema_agg["_vote_count"] = schema_agg["_bits_list"].map(len)
    schema_agg = schema_agg.sort_values(
        ["_vote_count", "_schema"], ascending=[False, True], kind="stable"
    )
    dominant = schema_agg.drop_duplicates(subset=["readname"], keep="first").copy()

    dominant["_intersection"]   = dominant["_bits_list"].map(_bitwise_intersection)
    dominant["_unique_codes"]   = dominant["_bits_list"].map(lambda bl: len(set(bl)))
    dominant["_has_clean"]      = dominant["_intersection"].str.contains("1", regex=False, na=False)
    dominant["_genuine_conflict"] = ~dominant["_has_clean"] & (dominant["_unique_codes"] > 1)

    # --- Vectorized non-conflict path (typically >99% of reads) -----------------
    nc = dominant[~dominant["_genuine_conflict"]]
    has_schema_nc = nc["_schema"].str.len() > 0
    consensus_nc = nc["_schema"].where(~has_schema_nc, nc["_schema"] + "::" + nc["_intersection"])
    consensus_nc = consensus_nc.where(has_schema_nc, nc["_intersection"])
    nc_result = pd.DataFrame({
        "readname":               nc["readname"].values,
        "read_consensus_code":    consensus_nc.values,
        "read_agreement_pct":     1.0,
        "read_is_mixed":          False,
        "read_used_intersection": True,
    })

    # --- Python loop only for the rare genuine-conflict reads -------------------
    conflict_records = []
    for _, row in dominant[dominant["_genuine_conflict"]].iterrows():
        schema    = row["_schema"]
        bits_list = row["_bits_list"]
        code_counts: dict[str, int] = {}
        for b in bits_list:
            code_counts[b] = code_counts.get(b, 0) + 1
        plurality_bits = min(code_counts, key=lambda k: (-code_counts[k], k))
        plurality_pct  = code_counts[plurality_bits] / len(bits_list)
        is_mixed = plurality_pct < per_read_min_agreement
        if not is_mixed:
            consensus_token = f"{schema}::{plurality_bits}" if schema else plurality_bits
        else:
            consensus_token = None
        conflict_records.append({
            "readname":               row["readname"],
            "read_consensus_code":    consensus_token,
            "read_agreement_pct":     plurality_pct,
            "read_is_mixed":          is_mixed,
            "read_used_intersection": False,
        })

    parts = [nc_result]
    if conflict_records:
        parts.append(pd.DataFrame(conflict_records))
    result = pd.concat(parts, ignore_index=True) if len(parts) > 1 else nc_result

    # Reads with only NaN/non-binary codes are not in dominant → mark mixed.
    covered = set(result["readname"])
    missing = [r for r in n_ctdmrs.index if r not in covered]
    if missing:
        result = pd.concat([result, pd.DataFrame({
            "readname":            missing,
            "read_consensus_code": None,
            "read_agreement_pct":  0.0,
            "read_is_mixed":       True,
            "read_used_intersection": False,
        })], ignore_index=True)

    result = result.merge(n_ctdmrs.reset_index(), on="readname", how="left")
    result["read_n_ctdmrs"]        = result["read_n_ctdmrs"].fillna(0).astype(int)
    result["read_is_mixed"]        = result["read_is_mixed"].astype(bool)
    result["read_used_intersection"] = result["read_used_intersection"].astype(bool)
    result["read_agreement_pct"]   = result["read_agreement_pct"].astype(float)
    return result


def assign_sv_celltypes(
    sv_df: pd.DataFrame,
    read_assignment_df: pd.DataFrame,
    *,
    min_overlap_pct: float = 0,
    min_agreement_pct: float = 0,
    window: int = 5000,
    breakpoint_exclusion_frac: float = 0.0,
    reads_col: str = "supporting_reads",
    sv_id_col: str = "id",
    unique_reads_for_overlap: bool = True,
    per_read_min_agreement: float = 0.66,
    use_read_names: bool = False,
) -> pd.DataFrame:
    """
    Link SVs to cell-type codes derived from per-read assignments and attach coordinates.

    Output includes:
      - SV coordinates and metadata (`sv_chr`, `sv_pos`, `vaf`, `sv_len`)
      - CpG coordinates (`cpg_chr`, `cpg_start`, `cpg_end`)
      - vote-based code summaries (`majority_code`, `assigned_code`, `code_counts`)
      - human-readable decoded summaries (`linked_celltypes`, counts/fractions, multi-link flag)
      - `n_mixed_reads`: reads excluded from voting due to conflicting ctDMR codes

    SV linking is coordinate-based (default) or read-name-based:
      - When use_read_names=False (default): SV and DMR chromosomes must match
        (chr-normalized), and a DMR is linked to an SV if it is within `window` bp
        padding around the SV interval.
      - When use_read_names=True: the spatial window is bypassed entirely. Each
        SV-supporting read is directly joined to its per-read classification using
        all ctDMR evidence genome-wide. This is correct when sv_df already carries
        explicit read name lists, because the reads were classified against
        targeted ctDMRs and the classifying ctDMR may lie
        anywhere along the read, not necessarily within `window` bp of the SV
        breakpoint.
      - ctDMRs overlapping the SV core expanded by
        `breakpoint_exclusion_frac * abs(SVLEN)` on each side are excluded.

    Per-read consensus (controlled by `per_read_min_agreement`):
      Before SV-level voting, each read's ctDMR codes are reduced to one consensus via
      bitwise intersection. If intersection collapses to all-zeros (genuine conflict),
      the plurality code is accepted only when its fraction >= per_read_min_agreement;
      otherwise the read is marked mixed and excluded from voting but counted in
      `n_mixed_reads`.
    """
    if window < 0:
        raise ValueError("window must be >= 0")
    breakpoint_exclusion_frac = validate_breakpoint_exclusion_frac(breakpoint_exclusion_frac)

    # ---- required SV columns ----
    needed = {"chr", "location", sv_id_col, reads_col}
    if not needed.issubset(sv_df.columns):
        missing = needed - set(sv_df.columns)
        raise ValueError(f"sv_df missing required columns: {sorted(missing)}")

    # ---- SV coords (chr/location) ----
    sv_coords = (
        sv_df[[sv_id_col, "chr", "location"]]
        .drop_duplicates(subset=[sv_id_col])
        .set_index(sv_id_col)
        .rename(columns={"chr": "sv_chr", "location": "sv_pos"})
    )
    sv_coords["sv_pos"] = pd.to_numeric(sv_coords["sv_pos"], errors="coerce").astype("Int64")

    # ---- optional SV extras (vaf/sv_len) ----
    cols_have = [c for c in ("vaf", "sv_len") if c in sv_df.columns]
    if cols_have:
        sv_extra = (
            sv_df[[sv_id_col] + cols_have]
            .drop_duplicates(subset=[sv_id_col])
            .set_index(sv_id_col)
        )
    else:
        sv_extra = pd.DataFrame(index=sv_coords.index)
    # ensure both columns exist
    for c in ("vaf", "sv_len"):
        if c not in sv_extra.columns:
            sv_extra[c] = pd.NA
    # types
    sv_extra["vaf"] = pd.to_numeric(sv_extra["vaf"], errors="coerce")
    sv_extra["sv_len"] = pd.to_numeric(sv_extra["sv_len"], errors="coerce").astype("Int64")

    # ---- normalize read assignment ----
    assignment = read_assignment_df.copy()
    assignment.index = assignment.index.astype(str)
    if "code" not in assignment.columns:
        raise ValueError("read_assignment_df must contain a 'code' column")
    assignment["code"] = assignment["code"].astype("string")
    if "best_group" not in assignment.columns:
        assignment["best_group"] = pd.Series(pd.NA, index=assignment.index, dtype="string")
    else:
        assignment["best_group"] = assignment["best_group"].astype("string")
    if "is_best_group" not in assignment.columns:
        assignment["is_best_group"] = False
    else:
        assignment["is_best_group"] = assignment["is_best_group"].astype("boolean").fillna(False).astype(bool)
    if "code_order" not in assignment.columns:
        assignment["code_order"] = pd.Series(pd.NA, index=assignment.index, dtype="string")
    else:
        assignment["code_order"] = assignment["code_order"].astype("string")
    if "best_group_leaves" not in assignment.columns:
        assignment["best_group_leaves"] = pd.Series(pd.NA, index=assignment.index, dtype="string")
    else:
        assignment["best_group_leaves"] = assignment["best_group_leaves"].astype("string")
    assignment["code"] = assignment.apply(
        lambda r: _normalize_binary_code(r.get("code", pd.NA), r.get("code_order", pd.NA)),
        axis=1,
    ).astype("string")
    for col in ("chr", "start", "end"):
        if col not in assignment.columns:
            assignment[col] = pd.NA

    # If multiple code schemas are present (e.g., multiple BEDs with different cell-type sets),
    # qualify code by schema so identical raw codes are not conflated.
    non_empty_orders = assignment["code_order"].dropna()
    non_empty_orders = non_empty_orders[non_empty_orders.str.len() > 0]
    use_qualified_code = non_empty_orders.nunique(dropna=True) > 1
    if use_qualified_code:
        has_order = assignment["code_order"].str.len().fillna(0) > 0
        assignment["code_token"] = assignment["code"]
        assignment.loc[has_order, "code_token"] = (
            assignment.loc[has_order, "code_order"] + "::" + assignment.loc[has_order, "code"]
        )
    else:
        assignment["code_token"] = assignment["code"]

    # ---- per-read consensus (globally across all ctDMRs, before SV-level voting) ----
    per_read_consensus = _compute_per_read_consensus_df(
        assignment,
        per_read_min_agreement=per_read_min_agreement,
    )

    # ---- SV supporting-read set ----
    sv = sv_df[[sv_id_col, reads_col]].copy()
    all_sv_ids = pd.Index(pd.unique(sv[sv_id_col]), name=sv_id_col)
    support_pairs = (
        sv.assign(_support_reads=sv[reads_col].apply(_parse_support_read_names))
          [[sv_id_col, "_support_reads"]]
          .explode("_support_reads")
          .dropna(subset=["_support_reads"])
    )
    if support_pairs.empty:
        support_pairs = pd.DataFrame(columns=[sv_id_col, "read"])
    else:
        support_pairs = (
            support_pairs.assign(read=support_pairs["_support_reads"].astype(str).str.strip())
            .loc[lambda d: d["read"].ne("")]
            [[sv_id_col, "read"]]
            .drop_duplicates(ignore_index=True)
        )

    n_supporting = (
        support_pairs.groupby(sv_id_col, sort=False)["read"]
        .nunique()
        .reindex(all_sv_ids, fill_value=0)
        .astype(int)
        .rename("n_supporting")
    )

    assignment["read"] = assignment.index.astype(str)
    sv_regions = _build_sv_regions(
        sv_df,
        sv_id_col=sv_id_col,
        breakpoint_exclusion_frac=breakpoint_exclusion_frac,
    )
    explicit_sv_col = sv_id_col if sv_id_col in assignment.columns else ("sv_id" if "sv_id" in assignment.columns else None)

    if use_read_names:
        # Bypass spatial window: directly join support reads with their full
        # per-read classification. Every ctDMR a supporting read overlaps
        # (anywhere in the genome) contributes evidence for the SV. This is
        # correct when reads were classified against targeted ctDMRs and the
        # classifying ctDMR may be
        # far from the SV breakpoint on a long read.
        if support_pairs.empty:
            mapped = assignment.iloc[0:0].copy()
            mapped[sv_id_col] = pd.Series(dtype="object")
        else:
            mapped = support_pairs.merge(assignment, on="read", how="inner")
    elif explicit_sv_col is not None:
        mapped = assignment[assignment[explicit_sv_col].notna()].copy()
        if explicit_sv_col != sv_id_col:
            mapped = mapped.rename(columns={explicit_sv_col: sv_id_col})
    else:
        assignment["chr_norm"] = assignment["chr"].map(_norm_chr)
        assignment["_start"] = pd.to_numeric(assignment["start"], errors="coerce")
        assignment["_end"] = pd.to_numeric(assignment["end"], errors="coerce")
        valid_region = (
            assignment["chr_norm"].ne("")
            & assignment["_start"].notna()
            & assignment["_end"].notna()
        )
        region_df = (
            assignment.loc[valid_region, ["chr_norm", "_start", "_end"]]
            .astype({"_start": "int64", "_end": "int64"})
            .drop_duplicates(ignore_index=True)
        )

        region_map = _build_sv_region_map(
            region_df,
            sv_regions,
            sv_id_col=sv_id_col,
            window=window,
        )
        if region_map.empty:
            mapped = assignment.iloc[0:0].copy()
            mapped[sv_id_col] = pd.Series(dtype="object")
        else:
            assignment_for_merge = assignment.copy()
            assignment_for_merge = assignment_for_merge.dropna(subset=["_start", "_end"])
            assignment_for_merge[["_start", "_end"]] = assignment_for_merge[["_start", "_end"]].astype("int64")
            mapped = assignment_for_merge.merge(
                region_map,
                on=["chr_norm", "_start", "_end"],
                how="inner",
            )

    if not mapped.empty:
        if breakpoint_exclusion_frac > 0.0:
            mapped = _filter_mapped_breakpoint_exclusion(
                mapped,
                sv_regions,
                sv_id_col=sv_id_col,
            )
        mapped["read"] = mapped["read"].astype(str)
        if not use_read_names:
            # spatial path: support_pairs filter not yet applied
            mapped = mapped.merge(support_pairs, on=[sv_id_col, "read"], how="inner")

    evidence_cols = [
        sv_id_col,
        "read",
        "code_token",
        "code",
        "code_order",
        "best_group",
        "best_group_leaves",
        "is_best_group",
        "chr",
        "start",
        "end",
    ]
    for optional_col in ("other_group", "other_group_leaves", "hyper_group_leaves", "hypo_group_leaves"):
        if optional_col in mapped.columns:
            evidence_cols.append(optional_col)

    evidence = mapped[mapped["code_token"].notna()][evidence_cols].drop_duplicates()

    # Resolve potentially multiple region-level votes for the same (SV, read).
    # n_overlapped always counts unique supporting reads with at least one linked evidence row.
    n_overlapped = evidence.groupby(sv_id_col, sort=False)["read"].nunique().rename("n_overlapped")

    # Apply per-read consensus: replace per-ctDMR codes with each read's consensus code,
    # exclude mixed reads from voting, and count them separately per SV.
    evidence_with_consensus = evidence.merge(
        per_read_consensus[["readname", "read_consensus_code", "read_is_mixed"]],
        left_on="read",
        right_on="readname",
        how="left",
    )
    _mixed_mask = evidence_with_consensus["read_is_mixed"].fillna(False)
    n_mixed_reads = (
        evidence_with_consensus[_mixed_mask]
        .groupby(sv_id_col, sort=False)["read"]
        .nunique()
        .reindex(all_sv_ids, fill_value=0)
        .astype(int)
        .rename("n_mixed_reads")
    )
    non_mixed = evidence_with_consensus[~_mixed_mask].copy()
    # Use consensus code; fall back to original code_token only when consensus is NA.
    non_mixed["code_token"] = non_mixed["read_consensus_code"].where(
        non_mixed["read_consensus_code"].notna(),
        non_mixed["code_token"],
    )
    # One vote per (SV, read): consensus already resolved multi-ctDMR duplicates.
    read_votes = (
        non_mixed[[sv_id_col, "read", "code_token"]]
        .drop_duplicates(subset=[sv_id_col, "read"])
    )
    n_vote_units = (
        read_votes.groupby(sv_id_col, sort=False).size()
        .reindex(all_sv_ids, fill_value=0)
        .rename("_n_vote_units")
    )

    # ---- code counts & majority ----
    if read_votes.empty:
        code_counts_str = pd.Series(dtype="string", name="code_counts")
        majority = pd.DataFrame(columns=["majority_code","majority_count"]).astype(
            {"majority_code":"string","majority_count":"int64"}
        )
        conflict_summary = _summarize_code_conflicts(read_votes, sv_id_col=sv_id_col)
        readable_summary = _summarize_celltype_links(
            evidence,
            sv_id_col=sv_id_col,
            unique_reads_for_overlap=unique_reads_for_overlap,
        )
    else:
        cc = (
            read_votes.assign(code_token=read_votes["code_token"].astype("string"))
                      .groupby([sv_id_col, "code_token"], sort=False)
                      .size().rename("count").reset_index()
        )
        cc = cc.sort_values(["count","code_token"], ascending=[False, True], kind="stable")
        code_counts_str = (
            cc.assign(pair=cc["code_token"] + ":" + cc["count"].astype(str))
              .groupby(sv_id_col, sort=False)["pair"].agg(";".join)
              .rename("code_counts").astype("string")
        )
        majority = (
            cc.drop_duplicates(subset=[sv_id_col], keep="first")
              .set_index(sv_id_col)
              .rename(columns={"code_token":"majority_code","count":"majority_count"})
        )
        majority["majority_code"] = majority["majority_code"].astype("string")
        conflict_summary = _summarize_code_conflicts(read_votes, sv_id_col=sv_id_col)
        readable_summary = _summarize_celltype_links(
            evidence,
            sv_id_col=sv_id_col,
            unique_reads_for_overlap=unique_reads_for_overlap,
        )
    # ---- CpG coords (mode) ----
    coord_df = (
        evidence[[sv_id_col, "chr", "start", "end"]]
        .assign(
            start=pd.to_numeric(evidence["start"], errors="coerce"),
            end=pd.to_numeric(evidence["end"], errors="coerce"),
        )
        .dropna(subset=["chr", "start", "end"])
    )

    if coord_df.empty:
        cpg_coords = pd.DataFrame(columns=["cpg_chr", "cpg_start", "cpg_end"]).set_index(
            pd.Index([], name=sv_id_col)
        )
    else:
        coord_counts = (
            coord_df
            .groupby([sv_id_col, "chr", "start", "end"], sort=False)
            .size()
            .rename("cnt")
            .reset_index()
            .sort_values(["cnt", "chr", "start", "end"], ascending=[False, True, True, True], kind="stable")
        )
        cpg_coords = (
            coord_counts
            .drop_duplicates(subset=[sv_id_col], keep="first")
            .set_index(sv_id_col)[["chr", "start", "end"]]
            .rename(columns={"chr": "cpg_chr", "start": "cpg_start", "end": "cpg_end"})
        )
    for col in ["cpg_start","cpg_end"]:
        if col in cpg_coords:
            cpg_coords[col] = pd.to_numeric(cpg_coords[col], errors="coerce").astype("Int64")

    # ---- assemble output ----
    out = (
        pd.DataFrame(index=n_supporting.index)
          .join(n_supporting)
          .join(n_overlapped)
          .join(n_mixed_reads)
          .join(n_vote_units)
          .join(majority[["majority_code","majority_count"]])
          .join(code_counts_str)
          .join(conflict_summary)
          .join(readable_summary)
          .reset_index()
          .rename(columns={sv_id_col: "id"})
    )

    out["n_overlapped"] = out["n_overlapped"].fillna(0).astype(int)
    out["n_mixed_reads"] = out["n_mixed_reads"].fillna(0).astype(int)
    out["_n_vote_units"] = out["_n_vote_units"].fillna(0).astype(int)
    out["majority_count"] = out["majority_count"].fillna(0).astype(int)
    out["majority_code"] = out.get("majority_code", pd.Series([pd.NA]*len(out))).astype("object")
    out.loc[out["majority_count"] == 0, "majority_code"] = pd.NA
    out["code_counts"] = out["code_counts"].fillna("").astype("string")
    if "has_hard_conflict" not in out.columns:
        out["has_hard_conflict"] = pd.Series(pd.array([pd.NA] * len(out), dtype="boolean"))
    else:
        out["has_hard_conflict"] = out["has_hard_conflict"].astype("boolean")
    if "intersection_code" not in out.columns:
        out["intersection_code"] = pd.Series([pd.NA] * len(out), dtype="object")
    out["intersection_code"] = out["intersection_code"].astype("object")
    for c in [
        "majority_schema",
        "primary_celltype",
        "majority_linked_celltypes",
        "majority_excluded_celltypes",
        "linked_celltypes",
        "linked_celltype_counts",
        "linked_celltype_fractions",
    ]:
        if c not in out.columns:
            out[c] = ""
        out[c] = out[c].fillna("").astype("string")
    if "is_multi_celltype_link" not in out.columns:
        out["is_multi_celltype_link"] = pd.Series(pd.array([pd.NA] * len(out), dtype="boolean"))
    else:
        out["is_multi_celltype_link"] = out["is_multi_celltype_link"].astype("boolean")

    out["overlap_pct"] = out["n_overlapped"] / out["n_supporting"].replace(0, 1)
    out["majority_pct"] = out["majority_count"] / out["_n_vote_units"].replace(0, 1)
    cond = (
        (out["overlap_pct"] >= min_overlap_pct)
        & (out["majority_pct"] >= min_agreement_pct)
        & (~out["has_hard_conflict"])
    )
    out["assigned_code"] = out["majority_code"].where(cond).astype("object")
    out = out.drop(columns=["_n_vote_units"])

    out = (
        out.set_index("id")
           .join(sv_coords, how="left")
           .join(cpg_coords, how="left")
           .join(sv_extra, how="left")
           .reset_index()
    )

    # final dtypes
    for col in ["sv_pos","cpg_start","cpg_end","sv_len"]:
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")

    return out[
        [
            "id","n_supporting","n_overlapped","n_mixed_reads","overlap_pct",
            "majority_code","majority_pct","assigned_code","code_counts",
            "majority_schema","primary_celltype","majority_linked_celltypes","majority_excluded_celltypes",
            "linked_celltypes","linked_celltype_counts","linked_celltype_fractions",
            "is_multi_celltype_link","has_hard_conflict","intersection_code",
            "sv_chr","sv_pos","cpg_chr","cpg_start","cpg_end",
            "vaf","sv_len",
        ]
    ]
