"""Assay-aware ctDMR summaries and native-layout rendering.

The legacy modifiedC renderer/extractor remain the default for legacy atlases.
Paired m/h widths encode relative assay means, not genomic subintervals or
canonical-C fractions. Missing channels never become zero-valued evidence.
"""

import numpy as np
import pandas as pd

from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam

COLORS = {"5mC": "#D84A45", "5hmC": "#287BC1", "modifiedC": "#88739B"}
READ_BOX_HEIGHT = 0.48


def has_separate_assays(dmrs):
    return "modification" in dmrs and dmrs.modification.isin(["5mC", "5hmC"]).any()


def marker_columns(frame, chromosome="chr"):
    return [chromosome, "start", "end"] + (
        ["modification"] if "modification" in frame else []
    )


def genomic_window(dmrs, chrom, start, end, count):
    """Include nearest N distinct intervals per side of the variant, plus base window."""
    if count < 0:
        raise ValueError("flanking_ctdmrs must be >= 0")
    if dmrs.empty:
        return start, end
    rows = dmrs[dmrs.chr_norm.eq(str(chrom).removeprefix("chr"))]
    intervals = rows[["start", "end"]].drop_duplicates()
    left = (
        intervals[intervals.end <= start].sort_values("end").tail(count)
        if count
        else intervals.iloc[:0]
    )
    right = (
        intervals[intervals.start >= end].sort_values("start").head(count)
        if count
        else intervals.iloc[:0]
    )
    selected = pd.concat([left, right])
    if selected.empty:
        return start, end
    return min(start, int(selected.start.min())), max(end, int(selected.end.max()))


def paired_fractions(m, h):
    if not np.isfinite(m) or not np.isfinite(h) or m < 0 or h < 0 or m + h <= 0:
        return None
    return m / (m + h), h / (m + h)


def interval_lanes(dmrs):
    """Greedy interval packing; exact m/h pairs share one interval/one lane."""
    intervals = sorted(set(zip(dmrs.start.astype(int), dmrs.end.astype(int))))
    ends = []
    lanes = {}
    for key in intervals:
        lane = next((i for i, end in enumerate(ends) if end <= key[0]), len(ends))
        if lane == len(ends):
            ends.append(key[1])
        else:
            ends[lane] = key[1]
        lanes[key] = lane
    return lanes, max(1, len(ends))


def extract_assay_means(
    *,
    sv_id,
    bam_path,
    reference_path,
    dmrs,
    support_assignment_df,
    decoded_assignment_df,
    logger,
    bam_handle=None,
    fasta_handle=None,
):
    """Use native extraction, matching observed CpGs before pairing assay means."""
    from . import viz as native

    outputs = []
    for chrom, chromosome_rows in dmrs.groupby("chr", sort=False):
        matrices = {}
        lo, hi = int(chromosome_rows.start.min()), int(chromosome_rows.end.max())
        for mod in chromosome_rows.modification.unique():
            matrix, _ = methyl_matrix_from_bam(
                bam_path,
                reference_path,
                chrom=str(chrom),
                start=lo,
                end=hi,
                return_positions=True,
                modification=mod,
                read_name_whitelist=set(support_assignment_df.read_name),
                bam_handle=bam_handle,
                fasta_handle=fasta_handle,
            )
            if isinstance(matrix.index, pd.MultiIndex):
                matrix.index = matrix.index.get_level_values("read_name").astype(str)
            # Match legacy treatment of repeated alignments of one read.
            matrices[mod] = matrix.groupby(level=0, sort=False).mean()
        pairs = {
            key
            for key, rows in chromosome_rows.groupby(["start", "end"])
            if set(rows.modification) == {"5mC", "5hmC"}
        }
        for mod, rows in chromosome_rows.groupby("modification", sort=False):
            # Reuse native metadata and missing-value contract without mixing assays.
            for row in rows.itertuples(index=False):
                columns = [
                    p for p in matrices[mod].columns if row.start <= int(p) < row.end
                ]
                for meta in support_assignment_df.itertuples(index=False):
                    values = (
                        matrices[mod].loc[meta.read_name, columns].dropna()
                        if meta.read_name in matrices[mod].index
                        else pd.Series(dtype=float)
                    )
                    paired = (row.start, row.end) in pairs
                    if paired:
                        other = matrices["5hmC" if mod == "5mC" else "5mC"]
                        available = (
                            other.loc[meta.read_name].dropna().index
                            if meta.read_name in other.index
                            else []
                        )
                        values = values.loc[values.index.intersection(available)]
                    matching = decoded_assignment_df
                    if not matching.empty:
                        matching = matching[
                            (matching.read_name == meta.read_name)
                            & (matching.chr_norm == native._norm_chr(chrom))
                            & (matching.start == row.start)
                            & (matching.end == row.end)
                        ]
                        if "modification" in matching:
                            matching = matching[matching.modification.eq(mod)]
                    assigned = (
                        "|".join(sorted(set(matching.assigned_celltypes)))
                        if not matching.empty
                        else ""
                    )
                    outputs.append(
                        dict(
                            sv_id=sv_id,
                            read_name=meta.read_name,
                            assignment_status=meta.assignment_status,
                            is_assigned=meta.is_assigned,
                            assigned_celltypes=meta.assigned_celltypes,
                            dmr_assigned_celltypes=assigned,
                            dmr_was_assigned=bool(assigned),
                            chr=chrom,
                            start=int(row.start),
                            end=int(row.end),
                            label=row.label,
                            best_group=row.best_group,
                            best_group_leaves=row.best_group_leaves,
                            best_dir=row.best_dir,
                            modification=mod,
                            mean_methylation=float(values.mean())
                            if len(values)
                            else np.nan,
                            n_cpg_observed=len(values),
                            n_cpg_in_dmr=len(columns),
                            paired_interval=paired,
                        )
                    )
    return pd.DataFrame(outputs)


def read_consensus(assignment_df, summary_path=None, min_agreement=0.66):
    """Saved consensus takes priority; fallback shares the scorer's exact code."""
    if summary_path:
        return pd.read_csv(summary_path, sep="\t", dtype={"readname": str})
    from sniffcell.deconv.deconv import _build_read_summary

    return _build_read_summary(assignment_df, per_read_min_agreement=min_agreement)


def apply_consensus(support, consensus):
    out = support.copy()
    if consensus.empty:
        return out
    calls = consensus.set_index("readname")
    for idx, row in out.iterrows():
        if row.read_name not in calls.index:
            continue
        call = calls.loc[row.read_name]
        mixed = str(call.is_mixed).lower() in {"true", "1"}
        label = call.get("linked_leaf_celltypes", "")
        label = "" if pd.isna(label) else str(label)
        assigned = bool(label) and not mixed
        out.loc[idx, "assigned_celltypes"] = label if assigned else ""
        out.loc[idx, "is_assigned"] = assigned
        out.loc[idx, "assignment_status"] = (
            "assigned"
            if assigned
            else "unassigned_mixed"
            if mixed
            else "unassigned_no_evidence"
        )
    return out


def _color(mod, value):
    from matplotlib.colors import to_rgb

    return tuple((1 - value) * np.ones(3) + value * np.array(to_rgb(COLORS[mod])))


def _draw_boxes(ax, a, b, y, height, values, *, paired, text=False, zorder=7):
    from matplotlib.patches import Rectangle

    if paired:
        m, h = values.get("5mC", np.nan), values.get("5hmC", np.nan)
        if not np.isfinite(m) or not np.isfinite(h):
            return
        fractions = paired_fractions(m, h)
        split = a + (b - a) * fractions[0] if fractions else a
        for left, right, mod, value in [(a, split, "5mC", m), (split, b, "5hmC", h)]:
            if fractions:
                ax.add_patch(
                    Rectangle(
                        (left, y - height / 2),
                        right - left,
                        height,
                        facecolor=COLORS[mod],
                        edgecolor="none",
                        zorder=zorder,
                    )
                )
            if text and right - left > 0.025 * (ax.get_xlim()[1] - ax.get_xlim()[0]):
                ax.text(
                    (left + right) / 2,
                    y,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white",
                    zorder=zorder + 2,
                )
        ax.add_patch(
            Rectangle(
                (a, y - height / 2),
                b - a,
                height,
                facecolor="none" if fractions else "#eeeeee",
                edgecolor="#333333",
                lw=0.7,
                zorder=zorder + 1,
            )
        )
    else:
        for mod, value in values.items():
            if not np.isfinite(value):
                continue
            ax.add_patch(
                Rectangle(
                    (a, y - height / 2),
                    b - a,
                    height,
                    facecolor=_color(mod, value),
                    edgecolor=COLORS[mod],
                    lw=0.7,
                    zorder=zorder,
                )
            )


def plot_assay_panel(**kwargs):
    """Use native read/CIGAR geometry with a genomic assay reference panel."""
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter
    from . import viz as native

    dmrs, methyl = kwargs["dmrs"], kwargs["methyl_df"]
    callouts = kwargs["linked_ctdmr_callouts"]
    # Even explicitly requested dual callouts retain genomic spacing in this view.
    if not callouts.empty:
        dmrs = pd.concat([dmrs, callouts]).drop_duplicates(marker_columns(dmrs))
        kwargs["region_start"] = min(kwargs["region_start"], int(dmrs.start.min()))
        kwargs["region_end"] = max(kwargs["region_end"], int(dmrs.end.max()))
    refs = native._reference_celltype_mean_columns(dmrs)
    lanes, n_lanes = interval_lanes(dmrs)
    spacing = max(1.0, n_lanes * 0.54 + 0.2)
    base = dict(
        kwargs,
        dmrs=dmrs,
        methyl_df=pd.DataFrame(),
        linked_ctdmr_callouts=pd.DataFrame(),
    )
    fig, axr, axd, tracks = native._plot_legacy_sv_panel(
        **base,
        _return_axes=True,
        _read_spacing=spacing,
        _assigned_support_color="#303840",
    )
    lo, hi = kwargs["region_start"], kwargs["region_end"]
    fig.set_size_inches(
        18, max(9, 5 + 0.16 * len(tracks) * spacing + 0.30 * len(refs) * n_lanes)
    )
    fig.subplots_adjust(left=0.12, right=0.79, top=0.79, bottom=0.16, hspace=0.24)
    axr.set_xlim(lo, hi)
    lookup = (
        {
            (
                str(r.read_name),
                int(r.start),
                int(r.end),
                str(r.modification),
            ): r.mean_methylation
            for r in methyl.itertuples(index=False)
        }
        if not methyl.empty
        else {}
    )
    groups = list(dmrs.groupby(["start", "end"], sort=True))
    for name, (y, read_start, read_end, supporting) in tracks.items():
        if not supporting:
            continue
        for (a, b), rows in groups:
            left, right = max(a, read_start), min(b, read_end)
            if right <= left:
                continue
            offset = (lanes[(a, b)] - (n_lanes - 1) / 2) * 0.54
            values = {
                mod: lookup.get((name, a, b, mod), np.nan) for mod in rows.modification
            }
            _draw_boxes(
                axr,
                left,
                right,
                y + offset,
                READ_BOX_HEIGHT,
                values,
                paired=set(values) == {"5mC", "5hmC"},
                text=True,
            )
    support = kwargs["support_assignment_df"].set_index("read_name")
    for name, (y, a, b, is_support) in tracks.items():
        if name not in support.index:
            continue
        row = support.loc[name]
        label = (
            str(row.assigned_celltypes)
            if row.is_assigned
            else "Ambiguous"
            if "mixed" in str(row.assignment_status)
            else "Unassigned"
        )
        axr.text(
            1.01,
            y,
            label,
            transform=axr.get_yaxis_transform(),
            va="center",
            fontsize=10,
            fontweight="bold" if is_support else "normal",
        )
        if is_support:
            axr.text(
                1.19,
                y,
                "Expanded"
                if str(kwargs["sv"].get("svtype", "")).startswith("expansion")
                else "Support",
                transform=axr.get_yaxis_transform(),
                va="center",
                fontsize=9,
                color="#5D3476",
            )
    axr.text(
        1.01,
        1.03,
        "Cell-type call",
        transform=axr.transAxes,
        fontsize=11,
        fontweight="bold",
    )
    axr.text(
        1.19, 1.03, "Variant", transform=axr.transAxes, fontsize=11, fontweight="bold"
    )
    axr.set_ylabel(
        "" if kwargs.get("applied_support_haplotype") else "Reads", fontsize=15
    )
    axr.tick_params(axis="x", labelbottom=False)
    axd.clear()
    axd.set_xlim(lo, hi)
    axd.axvspan(
        kwargs["sv"]["start"],
        kwargs["sv"]["end"],
        color="#3182bd",
        alpha=0.16,
        zorder=0,
    )
    row_height = max(1.0, n_lanes * 0.65)
    for (a, b), rows in groups:
        paired = set(rows.modification) == {"5mC", "5hmC"}
        lane_offset = (lanes[(a, b)] - (n_lanes - 1) / 2) * 0.65
        for i, ref in enumerate(refs):
            y = i * row_height + 0.5 + lane_offset
            values = {
                r.modification: float(r[ref]) if pd.notna(r[ref]) else np.nan
                for _, r in rows.iterrows()
            }
            _draw_boxes(axd, a, b, y, 0.48, values, paired=paired)
            fractions = (
                paired_fractions(values.get("5mC", np.nan), values.get("5hmC", np.nan))
                if paired
                else None
            )
            split = a + (b - a) * fractions[0] if fractions else a
            for mod, value in values.items():
                if not np.isfinite(value):
                    continue
                x0, x1 = (
                    ((a, split) if mod == "5mC" else (split, b))
                    if paired and fractions
                    else (a, b)
                )
                inside = x1 - x0 >= 0.020 * (hi - lo)
                text_y = y if inside else y + (-0.34 if mod != "5hmC" else 0.34)
                axd.text(
                    (x0 + x1) / 2,
                    text_y,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if paired and inside else "#333333",
                    zorder=10,
                )
        label = (
            "m | h"
            if paired
            else {"5mC": "m", "5hmC": "h", "modifiedC": "modC"}[
                rows.iloc[0].modification
            ]
        )
        axd.text((a + b) / 2, -0.38, label, ha="center", va="bottom", fontsize=8.5)
    axd.set_ylim(max(1, len(refs)) * row_height + 0.2, -0.7)
    axd.set_yticks(
        [i * row_height + 0.5 for i in range(len(refs))],
        [c[5:] for c in refs],
        fontsize=12,
    )
    axd.set_ylabel("Atlas reference", fontsize=14)
    axd.set_title(
        f"{len(groups)} genomic intervals · {len(dmrs)} assay rows · paired widths = relative 5mC : 5hmC",
        loc="left",
        fontsize=11,
        pad=12,
    )
    axd.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x):,}"))
    axd.tick_params(axis="x", labelsize=9)
    axd.set_xlabel(f"{kwargs['sv']['chrom']} coordinate (bp)", fontsize=13)
    axd.grid(axis="x", alpha=0.12)
    for legend in list(fig.legends):
        legend.remove()
    for item in list(fig.texts):
        if item is not fig._suptitle:
            item.remove()
    fig._suptitle.set_fontsize(17)
    if kwargs["sv"].get("sample_label"):
        fig._suptitle.set_text(
            kwargs["sv"]["sample_label"] + " · " + fig._suptitle.get_text()
        )
    excluded = kwargs["shown_reads"].attrs.get("excluded_read_count", 0)
    shown_support = sum(row[3] for row in tracks.values())
    hp = kwargs.get("applied_support_haplotype")
    fig.text(
        0.12,
        0.943,
        f"Haplotype: {'HP' + str(hp) if hp else 'all'} · {len(tracks)} reads shown · {shown_support}/{len(kwargs['sv']['supporting_reads'])} variant-support reads shown",
        fontsize=11,
    )
    fig.text(
        0.12,
        0.917,
        f"Locus: {kwargs['sv']['chrom']}:{kwargs['sv']['start'] + 1:,}–{kwargs['sv']['end']:,}"
        + (
            f" · Display-only exclusion: {excluded} read(s); calls unchanged"
            if excluded
            else ""
        ),
        fontsize=10,
    )
    handles = [
        Line2D([0], [0], color="#303840", marker=">", label="Assigned support"),
        Line2D([0], [0], color="#ff7f0e", ls="--", label="Unassigned support"),
        Line2D([0], [0], color="#bdbdbd", label="Other read"),
    ]
    handles += [Patch(facecolor=color, label=mod) for mod, color in COLORS.items()]
    handles += [
        Patch(facecolor="white", edgecolor="#555555", label="Low single-assay mean"),
        Patch(facecolor="#7b3294", label="CIGAR insertion"),
        Patch(facecolor="#fff7bc", edgecolor="#e66101", label="CIGAR deletion"),
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.115, 0.89),
        ncol=5,
        frameon=False,
        fontsize=9,
    )
    fig.text(
        0.12,
        0.075,
        "Paired boxes: red = m/(m+h), blue = h/(m+h); canonical C is excluded from this ratio. Numbers are absolute means.\nSingle-assay boxes: white-to-color mean. Equal read-box heights; overlapping intervals occupy separate lanes.",
        fontsize=9,
    )
    fig.text(
        0.12,
        0.028,
        "Read boxes summarize observed CpGs in ctDMRs; paired means use identical CpGs, without imputation.\nCell-type labels are predictions. CIGAR insertion width is not a direct repeat-length scale.",
        fontsize=9,
    )
    kwargs["output_path"].parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(kwargs["output_path"], dpi=kwargs["dpi"], bbox_inches="tight")
    plt.close(fig)
