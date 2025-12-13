import pandas as pd

def assign_sv_celltypes(
    sv_df: pd.DataFrame,
    read_assignment_df: pd.DataFrame,
    *,
    min_overlap_pct: float = 0,
    min_agreement_pct: float = 0,
    reads_col: str = "supporting_reads",
    sv_id_col: str = "id",
    unique_reads_for_overlap: bool = True,
) -> pd.DataFrame:
    """
    Link SVs to cell-type codes derived from per-read assignments and attach coordinates.

    Adds seven extra columns to the output:
      - sv_chr, sv_pos                 (from sv_df: 'chr' and 'location')
      - cpg_chr, cpg_start, cpg_end    (mode of (#chr, start, end) among reads linked to each SV)
      - vaf, sv_len                    (from sv_df if present; else <NA>)
    """

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
    for col in ("chr", "start", "end"):
        if col not in assignment.columns:
            assignment[col] = pd.NA

    # ---- helper ----
    def _len_support(x):
        if isinstance(x, (list, tuple)):
            return len(x)
        if isinstance(x, str):
            return len([t for t in x.split(",") if t]) if "," in x else (1 if x else 0)
        return 0

    sv = sv_df[[sv_id_col, reads_col]].copy()
    n_supporting = (
        sv.groupby(sv_id_col, sort=False)[reads_col]
          .apply(lambda s: int(s.apply(_len_support).sum()))
          .rename("n_supporting")
    )
    # print(sv)
    # explode to (sv_id, read)
    exploded = sv.assign(
        read=sv[reads_col].apply(
            lambda x: list(x) if isinstance(x, (list, tuple))
            else ([t for t in x.split(",") if t] if isinstance(x, str) else [])
        )
    )[[sv_id_col, "read"]].explode("read", ignore_index=True)

    if exploded.empty:
        out = sv[[sv_id_col]].drop_duplicates().rename(columns={sv_id_col: "id"})
        out["n_supporting"] = 0
        out["n_overlapped"] = 0
        out["overlap_pct"] = 0.0
        out["majority_code"] = pd.Series(pd.array([pd.NA] * len(out), dtype="string"))
        out["majority_pct"] = 0.0
        out["assigned_code"] = pd.Series(pd.array([pd.NA] * len(out), dtype="string"))
        out["code_counts"] = pd.Series(pd.array([""] * len(out), dtype="string"))
        out = (
            out.set_index("id")
               .join(sv_coords, how="left")
               .join(sv_extra, how="left")
               .reset_index()
        )
        out["cpg_chr"] = pd.NA
        out["cpg_start"] = pd.NA
        out["cpg_end"] = pd.NA
        return out[
            [
                "id","n_supporting","n_overlapped","overlap_pct",
                "majority_code","majority_pct","assigned_code","code_counts",
                "sv_chr","sv_pos","cpg_chr","cpg_start","cpg_end",
                "vaf","sv_len",
            ]
        ]

    exploded["read"] = exploded["read"].astype(str)
    merged = exploded.join(assignment[["code", "chr", "start", "end"]], on="read", how="left")

    if unique_reads_for_overlap:
        merged = merged.drop_duplicates([sv_id_col, "read"])

    overlapped = merged[merged["code"].notna()].copy()
    n_overlapped = overlapped.groupby(sv_id_col, sort=False).size().rename("n_overlapped")
    # ---- code counts & majority ----
    if overlapped.empty:
        code_counts_str = pd.Series(dtype="string", name="code_counts")
        majority = pd.DataFrame(columns=["majority_code","majority_count"]).astype(
            {"majority_code":"string","majority_count":"int64"}
        )
    else:
        cc = (
            overlapped.assign(code=overlapped["code"].astype("string"))
                      .groupby([sv_id_col, "code"], sort=False)
                      .size().rename("count").reset_index()
        )
        cc = cc.sort_values(["count","code"], ascending=[False, True], kind="stable")
        code_counts_str = (
            cc.assign(pair=cc["code"] + ":" + cc["count"].astype(str))
              .groupby(sv_id_col, sort=False)["pair"].agg(";".join)
              .rename("code_counts").astype("string")
        )
        # print(cc)
        majority = (
            cc.drop_duplicates(subset=[sv_id_col], keep="first")
              .set_index(sv_id_col)
              .rename(columns={"code":"majority_code","count":"majority_count"})
        )
        majority["majority_code"] = majority["majority_code"].astype("string")
    # ---- CpG coords (mode) ----
    coord_df = (
        merged[[sv_id_col, "chr", "start", "end"]]
        .assign(
            start=pd.to_numeric(merged["start"], errors="coerce"),
            end=pd.to_numeric(merged["end"], errors="coerce"),
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
          .join(majority[["majority_code","majority_count"]])
          .join(code_counts_str)
          .reset_index()
          .rename(columns={sv_id_col: "id"})
    )

    out["n_overlapped"] = out["n_overlapped"].fillna(0).astype(int)
    out["majority_count"] = out["majority_count"].fillna(0).astype(int)
    out["majority_code"] = out.get("majority_code", pd.Series([pd.NA]*len(out))).astype("object")
    out.loc[out["majority_count"] == 0, "majority_code"] = pd.NA
    out["code_counts"] = out["code_counts"].fillna("").astype("string")

    out["overlap_pct"] = out["n_overlapped"] / out["n_supporting"].replace(0, 1)
    out["majority_pct"] = out["majority_count"] / out["n_overlapped"].replace(0, 1)
    cond = (out["overlap_pct"] >= min_overlap_pct) & (out["majority_pct"] >= min_agreement_pct)
    out["assigned_code"] = out["majority_code"].where(cond).astype("object")

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
            "id","n_supporting","n_overlapped","overlap_pct",
            "majority_code","majority_pct","assigned_code","code_counts",
            "sv_chr","sv_pos","cpg_chr","cpg_start","cpg_end",
            "vaf","sv_len",
        ]
    ]
