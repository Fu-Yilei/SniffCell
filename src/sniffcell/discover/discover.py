from __future__ import annotations

import csv
import fcntl
import gzip
import hashlib
import json
import logging
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

from sniffcell.discover.harmonize_variants import write_harmonized_variants


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MOSAIC_FILTER_EXPR = "INFO/MOSAIC=1"
DEFAULT_SV_POST_MIN_DP = 5
DEFAULT_SV_POST_MIN_TARGET_ALT_AD = 2
DEFAULT_SV_POST_OTHER_MAX_ALT_AD = 0
DEFAULT_SNV_POST_MIN_DP = 5
DEFAULT_SNV_POST_MAX_DP = 50
DEFAULT_SNV_POST_MIN_DP_ABSENCE = 15
DEFAULT_SNV_POST_MIN_GQ = 21
DEFAULT_SNV_POST_MIN_OTHER_AF = 0.9
DEFAULT_STAGE_ORDER = (
    "sniffles",
    "sniffles_filter",
    "kanpig",
    "collapse",
    "medaka",
    "tdb_create",
    "tdb_merge",
    "clair3",
    "modkit",
)
# `trgt` is the PacBio HiFi-friendly TR genotyper. It sits parallel to `medaka`
# and is auto-substituted in when --clair3-platform=hifi; kept out of
# DEFAULT_STAGE_ORDER so the "all" alias stays a single TR-genotyping pass.
VALID_STAGES = DEFAULT_STAGE_ORDER + ("trgt",)
TR_GENOTYPER_STAGES = {"medaka", "trgt"}
STAGE_ALIASES = {
    "all": set(DEFAULT_STAGE_ORDER),
    "sv": {"sniffles", "sniffles_filter", "kanpig", "collapse"},
    "tdb": {"tdb_create", "tdb_merge"},
    "snv": {"clair3"},
    "mods": {"modkit"},
}
GROUP_SCOPED_STAGES = {
    "sniffles",
    "sniffles_filter",
    "kanpig",
    "medaka",
    "trgt",
    "tdb_create",
    "clair3",
    "modkit",
}
STAGE_RUNTIME_TOOLS = {
    "sniffles": {"sniffles", "bcftools"},
    "sniffles_filter": {"bcftools"},
    "kanpig": {"kanpig", "bcftools"},
    "collapse": {"bcftools", "truvari", "kanpig", "bgzip", "tabix"},
    "medaka": {"medaka"},
    "trgt": {"trgt", "samtools"},
    "tdb_create": {"tdb"},
    "tdb_merge": {"tdb"},
    "clair3": {"clair3"},
    "modkit": {"modkit", "tabix"},
}


@dataclass(frozen=True)
class SplitGroup:
    name: str
    bam_path: str
    bai_path: str
    read_summary_path: str | None = None


@dataclass
class RunContext:
    sample_id: str
    deconv_dir: Path
    split_dir: Path
    run_id: str
    run_root: Path
    manifest_dir: Path
    status_dir: Path
    commands_dir: Path
    logs_dir: Path
    slurm_dir: Path
    groups: list[SplitGroup]
    selected_groups: list[str]
    stages: tuple[str, ...]
    scheduler: str
    dry_run: bool
    force: bool
    rerun_failed: bool
    reference: Path
    tr_bed: Path
    sex: str
    tool_paths: dict[str, str]
    params: dict[str, object]


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _expand_path(value: str | os.PathLike[str]) -> Path:
    return Path(value).expanduser().resolve()


def _ensure_tool(tool: str, explicit_path: str | None = None) -> str:
    if explicit_path:
        path = _expand_path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(f"Executable not found: {path}")
        return str(path)
    found = shutil.which(tool)
    if not found:
        raise FileNotFoundError(f"Required executable not found in PATH: {tool}")
    return found


@lru_cache(maxsize=None)
def _tool_help_text(tool_path: str) -> str:
    try:
        completed = subprocess.run(
            [tool_path, "--help"],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return ""
    return "\n".join(part for part in (completed.stdout, completed.stderr) if part)


def _sniffles_supports_sample_name(tool_path: str) -> bool:
    return "--sample-name" in _tool_help_text(tool_path)


def _resolve_tool_optional(tool: str, explicit_path: str | None, required: bool) -> str:
    """Like _ensure_tool but only raises when required=True."""
    if explicit_path:
        path = _expand_path(explicit_path)
        if not path.exists():
            raise FileNotFoundError(f"Executable not found: {path}")
        return str(path)
    found = shutil.which(tool)
    if not found and required:
        raise FileNotFoundError(f"Required executable not found in PATH: {tool}")
    return found or tool


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _read_json(path: Path) -> dict:
    text = path.read_text()
    if not text.strip():
        return {}
    decoder = json.JSONDecoder()
    payload, _ = decoder.raw_decode(text)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return payload


def _command_hash(cmd: Sequence[str] | str) -> str:
    if isinstance(cmd, str):
        text = cmd
    else:
        text = shlex.join(cmd)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _write_command(path: Path, cmd: Sequence[str] | str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(cmd, str):
        text = cmd
    else:
        text = shlex.join(cmd)
    path.write_text(text + "\n")


def _run_command(
    *,
    cmd: Sequence[str] | str,
    stdout_path: Path,
    stderr_path: Path,
    dry_run: bool,
    shell: bool = False,
    env: dict[str, str] | None = None,
) -> None:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return
    with stdout_path.open("a") as stdout_handle, stderr_path.open("a") as stderr_handle:
        subprocess.run(cmd, shell=shell, check=True, stdout=stdout_handle, stderr=stderr_handle, env=env)


def _run_command_capture(cmd: Sequence[str]) -> str:
    return subprocess.check_output(cmd, text=True).strip()


def _normalize_chrom(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    lowered = text.lower()
    return lowered[3:] if lowered.startswith("chr") else lowered


def _vcf_has_records(vcf_path: Path) -> bool:
    openers = [gzip.open, open] if str(vcf_path).endswith(".gz") else [open]
    for opener in openers:
        try:
            with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    if line.strip():
                        return True
            return False
        except (gzip.BadGzipFile, OSError, EOFError):
            continue
    return False


def _text_file_has_records(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            return True
    return False


def _write_empty_vcf(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(header)
        return
    path.write_text(header, encoding="utf-8")


def _deconv_run_manifest_path(ctx: RunContext) -> Path:
    return ctx.deconv_dir / "deconv_run_manifest.json"


def _deconv_region_bed(ctx: RunContext) -> Path | None:
    manifest_path = _deconv_run_manifest_path(ctx)
    if not manifest_path.exists():
        return None
    payload = _read_json(manifest_path)
    regional_inputs = payload.get("regional_inputs")
    if not isinstance(regional_inputs, dict):
        return None
    expanded_bed = str(regional_inputs.get("expanded_bed", "")).strip()
    if not expanded_bed:
        return None
    candidate = _expand_path(expanded_bed)
    return candidate if candidate.exists() else None


def _resolve_targeted_tr_bed(ctx: RunContext) -> Path:
    region_bed = _deconv_region_bed(ctx)
    if region_bed is None:
        return ctx.tr_bed
    output_bed = ctx.manifest_dir / "targeted_tr.bed"
    if output_bed.exists():
        return output_bed

    regions_by_chrom: dict[str, list[tuple[int, int]]] = {}
    with region_bed.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = _normalize_chrom(fields[0])
            start = int(fields[1])
            end = int(fields[2])
            regions_by_chrom.setdefault(chrom, []).append((start, end))

    output_bed.parent.mkdir(parents=True, exist_ok=True)
    with ctx.tr_bed.open(encoding="utf-8") as source, output_bed.open("w", encoding="utf-8") as sink:
        for raw_line in source:
            line = raw_line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = _normalize_chrom(fields[0])
            if chrom not in regions_by_chrom:
                continue
            start = int(fields[1])
            end = int(fields[2])
            if any(start < region_end and end > region_start for region_start, region_end in regions_by_chrom[chrom]):
                sink.write(line + "\n")
    return output_bed


def _clair3_env(ctx: RunContext) -> dict[str, str]:
    env = os.environ.copy()
    clair3_bin = Path(ctx.tool_paths["clair3"]).expanduser().resolve()
    bin_dir = clair3_bin.parent
    existing_path = env.get("PATH", "")
    env["PATH"] = f"{bin_dir}{os.pathsep}{existing_path}" if existing_path else str(bin_dir)
    if bin_dir.name == "bin":
        env.setdefault("CONDA_PREFIX", str(bin_dir.parent))
    return env


def _require_paths(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required path(s): " + ", ".join(missing))


def _sanitize_token(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in value.strip())
    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")
    return cleaned.strip("_") or "group"


def _parse_stages(stage_text: str | None) -> tuple[str, ...]:
    if not stage_text:
        return DEFAULT_STAGE_ORDER
    requested: list[str] = []
    for raw in stage_text.split(","):
        token = raw.strip()
        if not token:
            continue
        if token in STAGE_ALIASES:
            for stage in DEFAULT_STAGE_ORDER:
                if stage in STAGE_ALIASES[token] and stage not in requested:
                    requested.append(stage)
            continue
        if token not in VALID_STAGES:
            raise ValueError(f"Unsupported stage: {token}")
        if token not in requested:
            requested.append(token)
    if not requested:
        raise ValueError("No valid stages were selected")
    return tuple(requested)


def _apply_platform_tr_substitution(
    stages: tuple[str, ...], platform: str | None
) -> tuple[str, ...]:
    """For HiFi platforms, swap medaka→trgt so a single TR genotyper runs.

    Only applies when medaka is present and trgt is not already requested.
    Returns the stages tuple unchanged when no substitution is needed.
    """
    if not platform or str(platform).lower() != "hifi":
        return stages
    if "trgt" in stages or "medaka" not in stages:
        return stages
    return tuple("trgt" if stage == "medaka" else stage for stage in stages)


def _required_tools_for_stages(stages: Sequence[str]) -> set[str]:
    required = {"python"}
    for stage in stages:
        required.update(STAGE_RUNTIME_TOOLS.get(stage, set()))
    return required


def _infer_sample_id(deconv_dir: Path) -> str:
    parent = deconv_dir.parent
    if parent.name:
        return parent.name
    raise ValueError(f"Could not infer sample ID from deconv dir: {deconv_dir}")


def _discover_groups(split_dir: Path) -> list[SplitGroup]:
    manifest_path = split_dir / "requested_group_splits.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing requested split manifest: {manifest_path}")
    with manifest_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No split BAM rows found in {manifest_path}")
    groups: list[SplitGroup] = []
    for row in rows:
        group_name = row.get("requested_group")
        bam_path = row.get("bam_path")
        if not group_name or not bam_path:
            raise ValueError(f"Manifest row is missing requested_group or bam_path: {row}")
        bam = _expand_path(bam_path)
        bai = _expand_path(str(bam) + ".bai")
        groups.append(
            SplitGroup(
                name=group_name,
                bam_path=str(bam),
                bai_path=str(bai),
                read_summary_path=row.get("read_summary_path") or None,
            )
        )
    return groups


def _select_groups(discovered: list[SplitGroup], group_text: str | None) -> list[str]:
    available = {group.name: group for group in discovered}
    if not group_text:
        if len(discovered) != 2:
            raise ValueError(
                "discover v1 requires exactly two split groups unless --groups is used"
            )
        return [group.name for group in discovered]
    requested = [token.strip() for token in group_text.split(",") if token.strip()]
    if not requested:
        raise ValueError("No valid group names were provided via --groups")
    missing = [group for group in requested if group not in available]
    if missing:
        raise ValueError(f"Unknown requested group(s): {', '.join(missing)}")
    return requested


def _build_context(args) -> RunContext:
    deconv_dir = _expand_path(args.deconv_dir)
    split_dir = _expand_path(args.split_dir) if args.split_dir else deconv_dir / "deconv_requested_group_splits"
    groups = _discover_groups(split_dir)
    selected_groups = _select_groups(groups, args.groups)
    if args.scheduler == "slurm" and len(selected_groups) != 2:
        raise ValueError("Slurm mode currently expects exactly two selected groups")
    sample_id = args.sample_id or _infer_sample_id(deconv_dir)
    run_id = args.run_id or f"discover_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    run_root = split_dir / "discover" / run_id
    stages = _parse_stages(args.stages)
    stages = _apply_platform_tr_substitution(stages, getattr(args, "clair3_platform", None))
    required_tools = _required_tools_for_stages(stages)
    tool_paths = {
        "bcftools": _resolve_tool_optional("bcftools", args.bcftools_bin, required="bcftools" in required_tools),
        "bgzip": _resolve_tool_optional("bgzip", args.bgzip_bin, required="bgzip" in required_tools),
        "kanpig": _resolve_tool_optional("kanpig", args.kanpig_bin, required="kanpig" in required_tools),
        "medaka": _resolve_tool_optional("medaka", args.medaka_bin, required="medaka" in required_tools),
        "modkit": _resolve_tool_optional("modkit", args.modkit_bin, required="modkit" in required_tools),
        "python": _ensure_tool("python", sys.executable),
        "samtools": _resolve_tool_optional("samtools", getattr(args, "samtools_bin", None), required="samtools" in required_tools),
        "sniffles": _resolve_tool_optional("sniffles", args.sniffles_bin, required="sniffles" in required_tools),
        "tabix": _resolve_tool_optional("tabix", args.tabix_bin, required="tabix" in required_tools),
        "tdb": _resolve_tool_optional("tdb", args.tdb_bin, required="tdb" in required_tools),
        "trgt": _resolve_tool_optional("trgt", getattr(args, "trgt_bin", None), required="trgt" in required_tools),
        "truvari": _resolve_tool_optional("truvari", args.truvari_bin, required="truvari" in required_tools),
        "clair3": _resolve_tool_optional("run_clair3.sh", args.clair3_bin, required="clair3" in required_tools),
    }
    reference = _expand_path(args.reference)
    tr_bed = _expand_path(args.tr_bed)
    _require_paths([deconv_dir, split_dir, reference, tr_bed])
    for group in groups:
        _require_paths([Path(group.bam_path), Path(group.bai_path)])
    params = {
        "collapse_use": args.collapse_use,
        "kanpig_passonly": bool(args.kanpig_passonly),
        "kanpig_sample_name_template": args.kanpig_sample_name_template,
        "kanpig_seqsim": args.kanpig_seqsim,
        "kanpig_sizesim": args.kanpig_sizesim,
        "medaka_model": args.medaka_model,
        "medaka_padding": args.medaka_padding,
        "medaka_phasing": args.medaka_phasing,
        "medaka_sample_name_template": args.medaka_sample_name_template,
        "mods_mode": args.mods_mode,
        "mosaic_filter_expression": args.sniffles_mosaic_filter_expression,
        "sniffles_cluster_merge_len": args.sniffles_cluster_merge_len,
        "slurm_account": args.slurm_account,
        "sniffles_include_germline": True,
        "sniffles_mosaic": True,
        "sniffles_output_rnames": True,
        "threads": args.threads,
        "tdb_create_force": bool(args.tdb_create_force),
        "tdb_create_mem": args.tdb_create_mem,
        "truvari_passonly": bool(args.truvari_passonly),
        "truvari_pctseq": args.truvari_pctseq,
        "truvari_pctsize": args.truvari_pctsize,
        "truvari_refdist": args.truvari_refdist,
        "clair3_platform": args.clair3_platform,
        "clair3_model_path": args.clair3_model_path,
        "trgt_sample_name_template": getattr(args, "trgt_sample_name_template", "{sample_id}.{group}"),
        "trgt_karyotype": getattr(args, "trgt_karyotype", None),
        "tr_margin_bp": getattr(args, "tr_margin_bp", 100),
        "tr_min_supporting_reads": getattr(args, "tr_min_supporting_reads", 2),
    }
    return RunContext(
        sample_id=sample_id,
        deconv_dir=deconv_dir,
        split_dir=split_dir,
        run_id=run_id,
        run_root=run_root,
        manifest_dir=run_root / "manifest",
        status_dir=run_root / "status",
        commands_dir=run_root / "commands",
        logs_dir=run_root / "logs",
        slurm_dir=run_root / "slurm",
        groups=groups,
        selected_groups=selected_groups,
        stages=stages,
        scheduler=args.scheduler,
        dry_run=bool(args.dry_run),
        force=bool(args.force),
        rerun_failed=bool(args.rerun_failed),
        reference=reference,
        tr_bed=tr_bed,
        sex=args.sex,
        tool_paths=tool_paths,
        params=params,
    )


def _group_lookup(ctx: RunContext, name: str) -> SplitGroup:
    for group in ctx.groups:
        if group.name == name:
            return group
    raise KeyError(name)


def _task_id(stage: str, group_name: str | None = None) -> str:
    if group_name is None:
        return stage
    return f"{stage}.{_sanitize_token(group_name)}"


def _task_paths(ctx: RunContext, stage: str, group_name: str | None = None) -> tuple[Path, Path, Path]:
    task = _task_id(stage, group_name)
    return (
        ctx.commands_dir / f"{task}.command.txt",
        ctx.logs_dir / f"{task}.out",
        ctx.logs_dir / f"{task}.err",
    )


def _done_path(ctx: RunContext, stage: str, group_name: str | None = None) -> Path:
    return ctx.status_dir / f"{_task_id(stage, group_name)}.done.json"


def _failed_path(ctx: RunContext, stage: str, group_name: str | None = None) -> Path:
    return ctx.status_dir / f"{_task_id(stage, group_name)}.failed.json"


def _all_exist(paths: Iterable[Path]) -> bool:
    return all(path.exists() for path in paths)


def _should_skip(
    *,
    ctx: RunContext,
    stage: str,
    outputs: Iterable[Path],
    group_name: str | None = None,
) -> bool:
    done = _done_path(ctx, stage, group_name)
    failed = _failed_path(ctx, stage, group_name)
    if ctx.force:
        return False
    if failed.exists() and not ctx.rerun_failed:
        raise RuntimeError(
            f"Task previously failed and rerun was not requested: {failed.name}"
        )
    return done.exists() and _all_exist(outputs)


def _record_done(
    *,
    ctx: RunContext,
    stage: str,
    outputs: Iterable[Path],
    command: Sequence[str] | str,
    group_name: str | None = None,
) -> None:
    payload = {
        "stage": stage,
        "group_name": group_name,
        "finished_at": _now_utc(),
        "outputs": [str(path) for path in outputs],
        "command_hash": _command_hash(command),
    }
    _write_json(_done_path(ctx, stage, group_name), payload)
    failed = _failed_path(ctx, stage, group_name)
    if failed.exists():
        failed.unlink()


def _record_failure(
    *,
    ctx: RunContext,
    stage: str,
    command: Sequence[str] | str,
    error: BaseException,
    group_name: str | None = None,
) -> None:
    payload = {
        "stage": stage,
        "group_name": group_name,
        "failed_at": _now_utc(),
        "error_type": type(error).__name__,
        "error": str(error),
        "command_hash": _command_hash(command),
    }
    _write_json(_failed_path(ctx, stage, group_name), payload)


def _symlink_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        dst.symlink_to(src)
    except OSError:
        shutil.copy2(src, dst)


def _sample_name(template: str, sample_id: str, group_name: str) -> str:
    return template.format(sample_id=sample_id, group=group_name)


def _sniffles_stage_dir(ctx: RunContext, group_name: str) -> Path:
    return ctx.run_root / "sv" / "sniffles" / _sanitize_token(group_name)


def _kanpig_stage_dir(ctx: RunContext, group_name: str) -> Path:
    return ctx.run_root / "sv" / "kanpig" / _sanitize_token(group_name)


def _collapse_stage_dir(ctx: RunContext) -> Path:
    return ctx.run_root / "sv" / "truvari_collapse" / (
        f"{_sanitize_token(ctx.selected_groups[0])}_vs_{_sanitize_token(ctx.selected_groups[1])}"
    )


def _sv_post_stage_dir(ctx: RunContext) -> Path:
    return ctx.run_root / "sv" / "sv_post_processing" / (
        f"{_sanitize_token(ctx.selected_groups[0])}_vs_{_sanitize_token(ctx.selected_groups[1])}"
    )


def _medaka_stage_dir(ctx: RunContext, group_name: str) -> Path:
    return ctx.run_root / "medaka_tandem" / f"{_sanitize_token(group_name)}.medaka"


def _trgt_stage_dir(ctx: RunContext, group_name: str) -> Path:
    return ctx.run_root / "trgt" / f"{_sanitize_token(group_name)}.trgt"


def _trgt_output_prefix(ctx: RunContext, group_name: str, sample_name: str) -> Path:
    return _trgt_stage_dir(ctx, group_name) / f"{sample_name}.trgt"


def _trgt_output_vcf(ctx: RunContext, group_name: str, sample_name: str) -> Path:
    return _trgt_output_prefix(ctx, group_name, sample_name).with_suffix(".vcf.gz")


def _trgt_output_spanning_sorted_bam(ctx: RunContext, group_name: str, sample_name: str) -> Path:
    prefix = _trgt_output_prefix(ctx, group_name, sample_name)
    return prefix.parent / f"{prefix.name}.spanning.sorted.bam"


def _karyotype_for_trgt(ctx: RunContext) -> str:
    override = ctx.params.get("trgt_karyotype")
    if override:
        return str(override)
    sex = (ctx.sex or "").lower()
    if sex == "male":
        return "XY"
    if sex == "female":
        return "XX"
    return "XX"


def _tdb_stage_dir(ctx: RunContext) -> Path:
    return ctx.run_root / "medaka_tandem" / "tdb"


def _tr_post_stage_dir(ctx: RunContext) -> Path:
    return ctx.run_root / "medaka_tandem" / "tr_post_processing" / (
        f"{_sanitize_token(ctx.selected_groups[0])}_vs_{_sanitize_token(ctx.selected_groups[1])}"
    )


def _snv_post_stage_dir(ctx: RunContext) -> Path:
    return ctx.run_root / "snv" / "snv_post_processing" / (
        f"{_sanitize_token(ctx.selected_groups[0])}_vs_{_sanitize_token(ctx.selected_groups[1])}"
    )


def _modkit_stage_dir(ctx: RunContext, group_name: str) -> Path:
    return ctx.run_root / "modkit" / _sanitize_token(group_name)


def _clair3_stage_dir(ctx: RunContext, group_name: str) -> Path:
    return ctx.run_root / "snv" / "clair3" / _sanitize_token(group_name)


def _existing_split_sniffles_vcf(ctx: RunContext, group_name: str) -> Path:
    return ctx.split_dir / f"{group_name}.sniffles.vcf.gz"


def _existing_split_medaka_vcf(ctx: RunContext, group_name: str) -> Path:
    return ctx.split_dir / "medaka_tandem" / f"{group_name}.medaka" / "medaka_to_ref.TR.vcf"


def _clair3_merge_output_path(ctx: RunContext, group_name: str) -> Path:
    return _clair3_stage_dir(ctx, group_name) / "merge_output.vcf.gz"


def _clair3_pileup_output_path(ctx: RunContext, group_name: str) -> Path:
    return _clair3_stage_dir(ctx, group_name) / "pileup.vcf.gz"


def _sniffles_output_paths(ctx: RunContext, group_name: str) -> tuple[Path, Path]:
    stage_dir = _sniffles_stage_dir(ctx, group_name)
    return stage_dir / "sniffles.raw.vcf.gz", stage_dir / "sniffles.raw.snf"


def _sniffles_filter_output_path(ctx: RunContext, group_name: str) -> Path:
    return _sniffles_stage_dir(ctx, group_name) / "sniffles.mosaic_only.vcf.gz"


def _kanpig_output_path(ctx: RunContext, group_name: str) -> Path:
    return _kanpig_stage_dir(ctx, group_name) / "kanpig.mosaic.vcf.gz"


def _medaka_output_vcf(ctx: RunContext, group_name: str) -> Path:
    return _medaka_stage_dir(ctx, group_name) / "medaka_to_ref.TR.vcf"


def _tdb_output_path(ctx: RunContext, group_name: str) -> Path:
    return _tdb_stage_dir(ctx) / f"{_sanitize_token(group_name)}.tdb"


def _merged_tdb_output_path(ctx: RunContext) -> Path:
    return _tdb_stage_dir(ctx) / f"{_sanitize_token(ctx.sample_id)}.merged.tdb"


def _existing_split_merged_tdb(ctx: RunContext) -> Path:
    return ctx.split_dir / "medaka_tandem" / f"{_sanitize_token(ctx.sample_id)}.medaka.tdb"


def _resolve_trimmed_reads_input(ctx: RunContext, group_name: str) -> Path:
    run_fasta = _medaka_stage_dir(ctx, group_name) / "trimmed_reads.fasta"
    if run_fasta.exists():
        return run_fasta
    existing_fasta = ctx.split_dir / "medaka_tandem" / f"{group_name}.medaka" / "trimmed_reads.fasta"
    if existing_fasta.exists():
        return existing_fasta
    return run_fasta


def _write_run_manifest(ctx: RunContext) -> None:
    payload = {
        "created_at": _now_utc(),
        "sample_id": ctx.sample_id,
        "deconv_dir": str(ctx.deconv_dir),
        "split_dir": str(ctx.split_dir),
        "run_id": ctx.run_id,
        "run_root": str(ctx.run_root),
        "scheduler": ctx.scheduler,
        "dry_run": ctx.dry_run,
        "force": ctx.force,
        "rerun_failed": ctx.rerun_failed,
        "stages": list(ctx.stages),
        "selected_groups": ctx.selected_groups,
        "reference": str(ctx.reference),
        "tr_bed": str(ctx.tr_bed),
        "sex": ctx.sex,
        "tool_paths": ctx.tool_paths,
        "params": ctx.params,
        "groups": [asdict(group) for group in ctx.groups],
    }
    _write_json(ctx.manifest_dir / "discover_run_manifest.json", payload)


def _write_task_manifest(ctx: RunContext, rows: list[dict[str, str]]) -> None:
    path = ctx.manifest_dir / "discover_task_manifest.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "task_id",
        "stage",
        "group_name",
        "scope",
        "inputs",
        "outputs",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _update_status(ctx: RunContext, task_id: str, payload: dict[str, object]) -> None:
    status_path = ctx.status_dir / "discover_status.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    with status_path.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        text = handle.read()
        if text.strip():
            decoder = json.JSONDecoder()
            state, _ = decoder.raw_decode(text)
            if not isinstance(state, dict):
                raise ValueError(f"Expected JSON object in {status_path}")
        else:
            state = {}
        state[task_id] = payload
        handle.seek(0)
        handle.truncate()
        handle.write(json.dumps(state, indent=2, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _selected_task_specs(ctx: RunContext) -> list[tuple[str, str | None]]:
    specs: list[tuple[str, str | None]] = []
    for stage in ctx.stages:
        if stage in GROUP_SCOPED_STAGES:
            specs.extend((stage, group_name) for group_name in ctx.selected_groups)
            if stage == "clair3" and len(ctx.selected_groups) == 2:
                specs.append(("snv_post_processing", None))
            continue
        if stage == "collapse":
            specs.extend(
                [
                    ("collapse_inputs", None),
                    ("collapse", None),
                    ("sv_post_processing", None),
                ]
            )
            continue
        if stage == "tdb_merge":
            specs.append(("tdb_merge", None))
            continue
        raise ValueError(f"Unsupported stage: {stage}")
    return specs


def _clear_force_rerun_state(ctx: RunContext) -> None:
    if not ctx.force:
        return
    status_path = ctx.status_dir / "discover_status.json"
    status_payload = _read_json(status_path) if status_path.exists() else {}
    changed = False
    seen: set[tuple[str, str | None]] = set()
    for stage, group_name in _selected_task_specs(ctx):
        spec = (stage, group_name)
        if spec in seen:
            continue
        seen.add(spec)
        task = _task_id(stage, group_name)
        if task in status_payload:
            status_payload.pop(task, None)
            changed = True
        for marker in (_done_path(ctx, stage, group_name), _failed_path(ctx, stage, group_name)):
            if marker.exists():
                marker.unlink()
                changed = True
    if changed:
        _write_json(status_path, status_payload)


def _run_task(
    *,
    ctx: RunContext,
    stage: str,
    command: Sequence[str] | str,
    outputs: Iterable[Path],
    group_name: str | None = None,
    shell: bool = False,
    env: dict[str, str] | None = None,
) -> None:
    outputs = list(outputs)
    task = _task_id(stage, group_name)
    command_path, stdout_path, stderr_path = _task_paths(ctx, stage, group_name)
    _write_command(command_path, command)
    if _should_skip(ctx=ctx, stage=stage, outputs=outputs, group_name=group_name):
        _update_status(
            ctx,
            task,
            {"state": "skipped", "updated_at": _now_utc(), "outputs": [str(x) for x in outputs]},
        )
        return
    _update_status(
        ctx,
        task,
        {"state": "running", "started_at": _now_utc(), "outputs": [str(x) for x in outputs]},
    )
    try:
        _run_command(
            cmd=command,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dry_run=ctx.dry_run,
            shell=shell,
            env=env,
        )
        if not ctx.dry_run:
            missing = [str(path) for path in outputs if not path.exists()]
            if missing:
                raise RuntimeError(f"Task completed without expected outputs: {', '.join(missing)}")
            _record_done(ctx=ctx, stage=stage, outputs=outputs, command=command, group_name=group_name)
            _update_status(
                ctx,
                task,
                {"state": "completed", "finished_at": _now_utc(), "outputs": [str(x) for x in outputs]},
            )
        else:
            _update_status(
                ctx,
                task,
                {"state": "dry_run", "finished_at": _now_utc(), "outputs": [str(x) for x in outputs]},
            )
    except Exception as exc:  # pragma: no cover - exercised via integration
        _record_failure(ctx=ctx, stage=stage, command=command, error=exc, group_name=group_name)
        _update_status(
            ctx,
            task,
            {"state": "failed", "finished_at": _now_utc(), "error": str(exc)},
        )
        raise


def _index_vcf(ctx: RunContext, vcf_path: Path) -> None:
    _run_task(
        ctx=ctx,
        stage="index_vcf",
        group_name=vcf_path.stem,
        command=[ctx.tool_paths["bcftools"], "index", "-t", "-f", str(vcf_path)],
        outputs=[Path(str(vcf_path) + ".tbi")],
    )


def _run_sniffles(ctx: RunContext, group_name: str) -> None:
    group = _group_lookup(ctx, group_name)
    stage_dir = _sniffles_stage_dir(ctx, group_name)
    stage_dir.mkdir(parents=True, exist_ok=True)
    output_vcf, output_snf = _sniffles_output_paths(ctx, group_name)
    sample_name = _sample_name(str(ctx.params["kanpig_sample_name_template"]), ctx.sample_id, group_name)
    cmd = [
        ctx.tool_paths["sniffles"],
        "--input",
        group.bam_path,
        "--reference",
        str(ctx.reference),
        "--vcf",
        str(output_vcf),
        "--snf",
        str(output_snf),
        "--threads",
        str(ctx.params["threads"]),
        "--mosaic",
        "--mosaic-include-germline",
        "--output-rnames",
        "--allow-overwrite",
        "--cluster-merge-len",
        str(ctx.params["sniffles_cluster_merge_len"]),
    ]
    if _sniffles_supports_sample_name(ctx.tool_paths["sniffles"]):
        cmd[9:9] = ["--sample-name", sample_name]
    _run_task(ctx=ctx, stage="sniffles", group_name=group_name, command=cmd, outputs=[output_vcf, output_snf])
    if not ctx.dry_run:
        subprocess.run(
            [ctx.tool_paths["bcftools"], "index", "-t", "-f", str(output_vcf)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def _resolve_filter_input(ctx: RunContext, group_name: str) -> Path:
    run_vcf, _ = _sniffles_output_paths(ctx, group_name)
    if run_vcf.exists():
        return run_vcf
    if ctx.dry_run:
        return run_vcf
    split_vcf = _existing_split_sniffles_vcf(ctx, group_name)
    if split_vcf.exists():
        return split_vcf
    raise FileNotFoundError(f"No Sniffles VCF found for group {group_name}")


def _run_sniffles_filter(ctx: RunContext, group_name: str) -> None:
    input_vcf = _resolve_filter_input(ctx, group_name)
    output_vcf = _sniffles_filter_output_path(ctx, group_name)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    expr = str(ctx.params["mosaic_filter_expression"])
    cmd = [
        ctx.tool_paths["bcftools"],
        "view",
        "-f",
        "PASS",
        "-i",
        expr,
        "-Oz",
        "-o",
        str(output_vcf),
        str(input_vcf),
    ]
    _run_task(ctx=ctx, stage="sniffles_filter", group_name=group_name, command=cmd, outputs=[output_vcf])
    if not ctx.dry_run:
        subprocess.run(
            [ctx.tool_paths["bcftools"], "index", "-t", "-f", str(output_vcf)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def _run_kanpig(ctx: RunContext, group_name: str) -> None:
    group = _group_lookup(ctx, group_name)
    input_vcf = _sniffles_filter_output_path(ctx, group_name)
    if not input_vcf.exists() and not ctx.dry_run:
        raise FileNotFoundError(f"Kanpig input is missing for {group_name}: {input_vcf}")
    stage_dir = _kanpig_stage_dir(ctx, group_name)
    stage_dir.mkdir(parents=True, exist_ok=True)
    raw_vcf = stage_dir / "kanpig.mosaic.vcf"
    output_vcf = _kanpig_output_path(ctx, group_name)
    rnames_tsv = stage_dir / "kanpig.rnames.tsv"
    sample_name = _sample_name(str(ctx.params["kanpig_sample_name_template"]), ctx.sample_id, group_name)
    cmd = [
        ctx.tool_paths["kanpig"],
        "mosaic",
        "--input",
        str(input_vcf),
        "--reference",
        str(ctx.reference),
        "--reads",
        group.bam_path,
        "--sample",
        sample_name,
        "--seqsim",
        str(ctx.params["kanpig_seqsim"]),
        "--sizesim",
        str(ctx.params["kanpig_sizesim"]),
        "-t",
        str(ctx.params["threads"]),
        "--rnames",
        str(rnames_tsv),
        "-o",
        str(raw_vcf),
    ]
    if ctx.params["kanpig_passonly"]:
        cmd.append("--passonly")
    outputs = [output_vcf, rnames_tsv]
    if not ctx.dry_run and not _vcf_has_records(input_vcf):
        task = _task_id("kanpig", group_name)
        command_path, _, _ = _task_paths(ctx, "kanpig", group_name)
        empty_cmd = ["skip-empty-kanpig-input", str(input_vcf), str(output_vcf)]
        _write_command(command_path, empty_cmd)
        if _should_skip(ctx=ctx, stage="kanpig", outputs=outputs, group_name=group_name):
            _update_status(
                ctx,
                task,
                {"state": "skipped", "updated_at": _now_utc(), "outputs": [str(x) for x in outputs]},
            )
            return
        _update_status(
            ctx,
            task,
            {
                "state": "running",
                "started_at": _now_utc(),
                "outputs": [str(x) for x in outputs],
                "empty_input_vcf": str(input_vcf),
            },
        )
        _symlink_or_copy(input_vcf, output_vcf)
        input_index = Path(str(input_vcf) + ".tbi")
        if input_index.exists():
            _symlink_or_copy(input_index, Path(str(output_vcf) + ".tbi"))
        rnames_tsv.write_text("", encoding="utf-8")
        _record_done(ctx=ctx, stage="kanpig", outputs=outputs, command=empty_cmd, group_name=group_name)
        _update_status(
            ctx,
            task,
            {
                "state": "completed_empty",
                "finished_at": _now_utc(),
                "outputs": [str(x) for x in outputs],
                "empty_input_vcf": str(input_vcf),
            },
        )
        return
    _run_task(ctx=ctx, stage="kanpig", group_name=group_name, command=cmd, outputs=[raw_vcf, rnames_tsv])
    if not ctx.dry_run:
        subprocess.run(
            [ctx.tool_paths["bcftools"], "sort", "-Oz", "-o", str(output_vcf), str(raw_vcf)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            [ctx.tool_paths["bcftools"], "index", "-t", "-f", str(output_vcf)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        raw_vcf.unlink(missing_ok=True)


def _resolve_collapse_input(ctx: RunContext, group_name: str) -> Path:
    if ctx.params["collapse_use"] == "sniffles":
        candidate = _sniffles_filter_output_path(ctx, group_name)
        if candidate.exists() or ctx.dry_run:
            return candidate
        raise FileNotFoundError(f"Collapse input is missing for {group_name}: {candidate}")
    candidate = _kanpig_output_path(ctx, group_name)
    if candidate.exists() or ctx.dry_run:
        return candidate
    raise FileNotFoundError(f"Collapse input is missing for {group_name}: {candidate}")


def _run_collapse(ctx: RunContext) -> None:
    stage_dir = _collapse_stage_dir(ctx)
    stage_dir.mkdir(parents=True, exist_ok=True)
    input_a = _resolve_collapse_input(ctx, ctx.selected_groups[0])
    input_b = _resolve_collapse_input(ctx, ctx.selected_groups[1])
    sampleless_a = stage_dir / f"{_sanitize_token(ctx.selected_groups[0])}.sites.vcf.gz"
    sampleless_b = stage_dir / f"{_sanitize_token(ctx.selected_groups[1])}.sites.vcf.gz"
    merged_input = stage_dir / "collapse.inputs.vcf.gz"
    raw_output = stage_dir / "collapsed.vcf"
    output_vcf = stage_dir / "collapsed.sorted.vcf.gz"
    removed_vcf = stage_dir / "removed.vcf"
    merge_cmd = "\n".join(
        [
            "set -euo pipefail",
            shlex.join(
                [
                    ctx.tool_paths["bcftools"],
                    "view",
                    "-G",
                    "-Oz",
                    "-o",
                    str(sampleless_a),
                    str(input_a),
                ]
            ),
            shlex.join([ctx.tool_paths["bcftools"], "index", "-t", "-f", str(sampleless_a)]),
            shlex.join(
                [
                    ctx.tool_paths["bcftools"],
                    "view",
                    "-G",
                    "-Oz",
                    "-o",
                    str(sampleless_b),
                    str(input_b),
                ]
            ),
            shlex.join([ctx.tool_paths["bcftools"], "index", "-t", "-f", str(sampleless_b)]),
            shlex.join(
                [
                    ctx.tool_paths["bcftools"],
                    "concat",
                    "-a",
                    "-Oz",
                    "-o",
                    str(merged_input),
                    str(sampleless_a),
                    str(sampleless_b),
                ]
            ),
        ]
    )
    _run_task(ctx=ctx, stage="collapse_inputs", command=merge_cmd, outputs=[merged_input], shell=True)
    if not ctx.dry_run:
        subprocess.run(
            [ctx.tool_paths["bcftools"], "index", "-t", "-f", str(merged_input)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    cmd = [
        ctx.tool_paths["truvari"],
        "collapse",
        "-i",
        str(merged_input),
        "-o",
        str(raw_output),
        "-c",
        str(removed_vcf),
        "-f",
        str(ctx.reference),
        "-r",
        str(ctx.params["truvari_refdist"]),
        "-p",
        str(ctx.params["truvari_pctseq"]),
        "-P",
        str(ctx.params["truvari_pctsize"]),
    ]
    if ctx.params["truvari_passonly"]:
        cmd.append("--passonly")
    _run_task(ctx=ctx, stage="collapse", command=cmd, outputs=[raw_output, removed_vcf])
    if not ctx.dry_run:
        subprocess.run(
            [ctx.tool_paths["bcftools"], "sort", "-Oz", "-o", str(output_vcf), str(raw_output)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            [ctx.tool_paths["bcftools"], "index", "-t", "-f", str(output_vcf)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        raw_output.unlink(missing_ok=True)
        _write_json(
            stage_dir / "collapse.report.json",
            {
                "created_at": _now_utc(),
                "input_a": str(input_a),
                "input_b": str(input_b),
                "collapsed_vcf": str(output_vcf),
                "removed_vcf": str(removed_vcf),
            },
        )
        _run_sv_post_processing(ctx)


def _run_sv_post_processing(ctx: RunContext) -> None:
    if len(ctx.selected_groups) != 2:
        raise ValueError("sv_post_processing requires exactly two selected groups")
    for group_name in ctx.selected_groups:
        split_vcf = _existing_split_sniffles_vcf(ctx, group_name)
        if split_vcf.exists():
            continue
        run_vcf, _ = _sniffles_output_paths(ctx, group_name)
        if not run_vcf.exists():
            raise FileNotFoundError(
                f"sv_post_processing could not find a Sniffles VCF for group {group_name}"
            )
        _symlink_or_copy(run_vcf, split_vcf)
    output_dir = _sv_post_stage_dir(ctx)
    summary_path = output_dir / "summary.json"
    cmd = [
        sys.executable,
        "-m",
        "sniffcell.discover.sv_post_processing",
        "--split-dir",
        str(ctx.split_dir),
        "--reference",
        str(ctx.reference),
        "--groups",
        ",".join(ctx.selected_groups),
        "--output-dir",
        str(output_dir),
        "--bcftools-bin",
        ctx.tool_paths["bcftools"],
        "--bgzip-bin",
        ctx.tool_paths["bgzip"],
        "--tabix-bin",
        ctx.tool_paths["tabix"],
        "--truvari-bin",
        ctx.tool_paths["truvari"],
        "--kanpig-bin",
        ctx.tool_paths["kanpig"],
        "--threads",
        str(ctx.params["threads"]),
        "--kanpig-seqsim",
        str(ctx.params["kanpig_seqsim"]),
        "--kanpig-sizesim",
        str(ctx.params["kanpig_sizesim"]),
        "--mosaic-filter",
        "--mosaic-filter-expression",
        str(ctx.params["mosaic_filter_expression"]),
        "--min-dp",
        str(DEFAULT_SV_POST_MIN_DP),
        "--min-target-alt-ad",
        str(DEFAULT_SV_POST_MIN_TARGET_ALT_AD),
        "--other-max-alt-ad",
        str(DEFAULT_SV_POST_OTHER_MAX_ALT_AD),
        "--sample-id",
        ctx.sample_id,
    ]
    _run_task(
        ctx=ctx,
        stage="sv_post_processing",
        command=cmd,
        outputs=[summary_path],
    )


def _write_empty_medaka_outputs(stage_dir: Path, output_vcf: Path) -> None:
    stage_dir.mkdir(parents=True, exist_ok=True)
    _write_empty_vcf(output_vcf)
    trimmed_reads = stage_dir / "trimmed_reads.fasta"
    if not trimmed_reads.exists():
        trimmed_reads.write_text("", encoding="utf-8")


def _run_medaka(ctx: RunContext, group_name: str) -> None:
    group = _group_lookup(ctx, group_name)
    stage_dir = _medaka_stage_dir(ctx, group_name)
    stage_dir.parent.mkdir(parents=True, exist_ok=True)
    output_vcf = stage_dir / "medaka_to_ref.TR.vcf"
    existing_vcf = _existing_split_medaka_vcf(ctx, group_name)
    if existing_vcf.exists() and not ctx.force:
        command = ["reuse-existing-medaka-vcf", str(existing_vcf), str(output_vcf)]
        command_path, _, _ = _task_paths(ctx, "medaka", group_name)
        _write_command(command_path, command)
        if _should_skip(ctx=ctx, stage="medaka", outputs=[output_vcf], group_name=group_name):
            _update_status(
                ctx,
                _task_id("medaka", group_name),
                {"state": "skipped", "updated_at": _now_utc(), "outputs": [str(output_vcf)]},
            )
            return
        _update_status(
            ctx,
            _task_id("medaka", group_name),
            {
                "state": "running",
                "started_at": _now_utc(),
                "outputs": [str(output_vcf)],
                "reused_from": str(existing_vcf),
            },
        )
        try:
            if not ctx.dry_run:
                _symlink_or_copy(existing_vcf, output_vcf)
                _record_done(ctx=ctx, stage="medaka", outputs=[output_vcf], command=command, group_name=group_name)
                _update_status(
                    ctx,
                    _task_id("medaka", group_name),
                    {
                        "state": "completed",
                        "finished_at": _now_utc(),
                        "outputs": [str(output_vcf)],
                        "reused_from": str(existing_vcf),
                    },
                )
            else:
                _update_status(
                    ctx,
                    _task_id("medaka", group_name),
                    {"state": "dry_run", "finished_at": _now_utc(), "outputs": [str(output_vcf)]},
                )
            return
        except Exception as exc:  # pragma: no cover - exercised via integration
            _record_failure(ctx=ctx, stage="medaka", command=command, error=exc, group_name=group_name)
            _update_status(
                ctx,
                _task_id("medaka", group_name),
                {"state": "failed", "finished_at": _now_utc(), "error": str(exc)},
            )
            raise
    effective_tr_bed = _resolve_targeted_tr_bed(ctx)
    if not ctx.dry_run and not _text_file_has_records(effective_tr_bed):
        command = ["skip-empty-tr-bed", str(effective_tr_bed), str(output_vcf)]
        task = _task_id("medaka", group_name)
        command_path, _, _ = _task_paths(ctx, "medaka", group_name)
        _write_command(command_path, command)
        _update_status(
            ctx,
            task,
            {
                "state": "running",
                "started_at": _now_utc(),
                "outputs": [str(output_vcf)],
                "targeted_tr_bed": str(effective_tr_bed),
            },
        )
        _write_empty_medaka_outputs(stage_dir, output_vcf)
        _record_done(ctx=ctx, stage="medaka", outputs=[output_vcf], command=command, group_name=group_name)
        _update_status(
            ctx,
            task,
            {
                "state": "completed_empty",
                "finished_at": _now_utc(),
                "outputs": [str(output_vcf)],
                "targeted_tr_bed": str(effective_tr_bed),
            },
        )
        return
    sample_name = _sample_name(str(ctx.params["medaka_sample_name_template"]), ctx.sample_id, group_name)
    cmd = [
        ctx.tool_paths["medaka"],
        "tandem",
        "--workers",
        str(ctx.params["threads"]),
        "--model",
        str(ctx.params["medaka_model"]),
        "--sample_name",
        sample_name,
        "--padding",
        str(ctx.params["medaka_padding"]),
        group.bam_path,
        str(ctx.reference),
        str(effective_tr_bed),
        ctx.sex,
        str(stage_dir),
    ]
    if ctx.params["medaka_phasing"]:
        cmd[2:2] = ["--phasing", str(ctx.params["medaka_phasing"])]
    task = _task_id("medaka", group_name)
    command_path, stdout_path, stderr_path = _task_paths(ctx, "medaka", group_name)
    _write_command(command_path, cmd)
    if _should_skip(ctx=ctx, stage="medaka", outputs=[output_vcf], group_name=group_name):
        _update_status(
            ctx,
            task,
            {"state": "skipped", "updated_at": _now_utc(), "outputs": [str(output_vcf)]},
        )
        return
    _update_status(
        ctx,
        task,
        {
            "state": "running",
            "started_at": _now_utc(),
            "outputs": [str(output_vcf)],
            "targeted_tr_bed": str(effective_tr_bed),
        },
    )
    try:
        _run_command(
            cmd=cmd,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dry_run=ctx.dry_run,
        )
        if not ctx.dry_run:
            if not output_vcf.exists():
                _write_empty_medaka_outputs(stage_dir, output_vcf)
            _record_done(ctx=ctx, stage="medaka", outputs=[output_vcf], command=cmd, group_name=group_name)
            _update_status(
                ctx,
                task,
                {
                    "state": ("completed" if _vcf_has_records(output_vcf) else "completed_empty"),
                    "finished_at": _now_utc(),
                    "outputs": [str(output_vcf)],
                    "targeted_tr_bed": str(effective_tr_bed),
                },
            )
        else:
            _update_status(
                ctx,
                task,
                {"state": "dry_run", "finished_at": _now_utc(), "outputs": [str(output_vcf)]},
            )
    except Exception as exc:  # pragma: no cover - exercised via integration
        _record_failure(ctx=ctx, stage="medaka", command=cmd, error=exc, group_name=group_name)
        _update_status(
            ctx,
            task,
            {"state": "failed", "finished_at": _now_utc(), "error": str(exc)},
        )
        raise


def _run_trgt(ctx: RunContext, group_name: str) -> None:
    """Run trgt genotype on the per-group split BAM, then sort/index the spanning BAM.

    Mirrors _run_medaka's bookkeeping (status/done/failed records and skip logic)
    but invokes PacBio's trgt instead of medaka tandem.
    """
    group = _group_lookup(ctx, group_name)
    sample_name = _sample_name(
        str(ctx.params["trgt_sample_name_template"]), ctx.sample_id, group_name
    )
    stage_dir = _trgt_stage_dir(ctx, group_name)
    stage_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = _trgt_output_prefix(ctx, group_name, sample_name)
    output_vcf = _trgt_output_vcf(ctx, group_name, sample_name)
    spanning_bam = output_prefix.parent / f"{output_prefix.name}.spanning.bam"
    spanning_sorted = _trgt_output_spanning_sorted_bam(ctx, group_name, sample_name)
    spanning_sorted_bai = spanning_sorted.with_suffix(spanning_sorted.suffix + ".bai")

    effective_tr_bed = _resolve_targeted_tr_bed(ctx)
    task = _task_id("trgt", group_name)
    command_path, stdout_path, stderr_path = _task_paths(ctx, "trgt", group_name)

    if not ctx.dry_run and not _text_file_has_records(effective_tr_bed):
        empty_cmd = ["skip-empty-tr-bed", str(effective_tr_bed), str(output_vcf)]
        _write_command(command_path, empty_cmd)
        _update_status(
            ctx,
            task,
            {
                "state": "running",
                "started_at": _now_utc(),
                "outputs": [str(output_vcf), str(spanning_sorted)],
                "targeted_tr_bed": str(effective_tr_bed),
            },
        )
        _write_empty_vcf(output_vcf)
        spanning_sorted.write_bytes(b"")
        spanning_sorted_bai.write_bytes(b"")
        _record_done(
            ctx=ctx,
            stage="trgt",
            outputs=[output_vcf, spanning_sorted],
            command=empty_cmd,
            group_name=group_name,
        )
        _update_status(
            ctx,
            task,
            {
                "state": "completed_empty",
                "finished_at": _now_utc(),
                "outputs": [str(output_vcf), str(spanning_sorted)],
                "targeted_tr_bed": str(effective_tr_bed),
            },
        )
        return

    karyotype = _karyotype_for_trgt(ctx)
    threads = str(ctx.params["threads"])
    cmd = [
        ctx.tool_paths["trgt"],
        "genotype",
        "--genome",
        str(ctx.reference),
        "--reads",
        group.bam_path,
        "--repeats",
        str(effective_tr_bed),
        "--output-prefix",
        str(output_prefix),
        "--karyotype",
        karyotype,
        "--sample-name",
        sample_name,
        "--threads",
        threads,
    ]
    _write_command(command_path, cmd)
    outputs = [output_vcf, spanning_sorted]
    if _should_skip(ctx=ctx, stage="trgt", outputs=outputs, group_name=group_name):
        _update_status(
            ctx,
            task,
            {"state": "skipped", "updated_at": _now_utc(), "outputs": [str(p) for p in outputs]},
        )
        return
    _update_status(
        ctx,
        task,
        {
            "state": "running",
            "started_at": _now_utc(),
            "outputs": [str(p) for p in outputs],
            "targeted_tr_bed": str(effective_tr_bed),
        },
    )
    try:
        _run_command(
            cmd=cmd,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            dry_run=ctx.dry_run,
        )
        if not ctx.dry_run:
            if spanning_bam.exists():
                samtools = ctx.tool_paths["samtools"]
                sort_cmd = [samtools, "sort", "-@", threads, "-o", str(spanning_sorted), str(spanning_bam)]
                _run_command(cmd=sort_cmd, stdout_path=stdout_path, stderr_path=stderr_path, dry_run=False)
                index_cmd = [samtools, "index", "-@", threads, str(spanning_sorted)]
                _run_command(cmd=index_cmd, stdout_path=stdout_path, stderr_path=stderr_path, dry_run=False)
            _record_done(
                ctx=ctx,
                stage="trgt",
                outputs=outputs,
                command=cmd,
                group_name=group_name,
            )
            _update_status(
                ctx,
                task,
                {
                    "state": "completed" if output_vcf.exists() else "completed_empty",
                    "finished_at": _now_utc(),
                    "outputs": [str(p) for p in outputs],
                    "targeted_tr_bed": str(effective_tr_bed),
                },
            )
        else:
            _update_status(
                ctx,
                task,
                {"state": "dry_run", "finished_at": _now_utc(), "outputs": [str(p) for p in outputs]},
            )
    except Exception as exc:  # pragma: no cover - exercised via integration
        _record_failure(ctx=ctx, stage="trgt", command=cmd, error=exc, group_name=group_name)
        _update_status(
            ctx,
            task,
            {"state": "failed", "finished_at": _now_utc(), "error": str(exc)},
        )
        raise


def _resolve_tdb_input(ctx: RunContext, group_name: str) -> Path:
    run_vcf = _medaka_output_vcf(ctx, group_name)
    if run_vcf.exists():
        return run_vcf
    if ctx.dry_run:
        return run_vcf
    existing_vcf = _existing_split_medaka_vcf(ctx, group_name)
    if existing_vcf.exists():
        return existing_vcf
    raise FileNotFoundError(f"No Medaka TR VCF found for group {group_name}")


def _run_tdb_create(ctx: RunContext, group_name: str) -> None:
    input_vcf = _resolve_tdb_input(ctx, group_name)
    output_tdb = _tdb_output_path(ctx, group_name)
    output_tdb.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        ctx.tool_paths["tdb"],
        "create",
        "-o",
        str(output_tdb),
        "--mem",
        str(ctx.params["tdb_create_mem"]),
    ]
    if ctx.params["tdb_create_force"]:
        cmd.append("--force")
    cmd.append(str(input_vcf))
    _run_task(ctx=ctx, stage="tdb_create", group_name=group_name, command=cmd, outputs=[output_tdb])


def _run_tdb_merge(ctx: RunContext) -> None:
    output_tdb = _merged_tdb_output_path(ctx)
    output_tdb.parent.mkdir(parents=True, exist_ok=True)
    input_a = _tdb_output_path(ctx, ctx.selected_groups[0])
    input_b = _tdb_output_path(ctx, ctx.selected_groups[1])
    medaka_a = _resolve_tdb_input(ctx, ctx.selected_groups[0])
    medaka_b = _resolve_tdb_input(ctx, ctx.selected_groups[1])
    if not ctx.dry_run and (not _vcf_has_records(medaka_a) or not _vcf_has_records(medaka_b)):
        command = ["skip-empty-tdb-merge", str(medaka_a), str(medaka_b)]
        task = _task_id("tdb_merge", None)
        command_path, _, _ = _task_paths(ctx, "tdb_merge", None)
        _write_command(command_path, command)
        _update_status(
            ctx,
            task,
            {
                "state": "running",
                "started_at": _now_utc(),
                "empty_medaka_inputs": [str(medaka_a), str(medaka_b)],
            },
        )
        _record_done(ctx=ctx, stage="tdb_merge", outputs=[], command=command)
        _update_status(
            ctx,
            task,
            {
                "state": "completed_empty",
                "finished_at": _now_utc(),
                "empty_medaka_inputs": [str(medaka_a), str(medaka_b)],
            },
        )
        _run_tr_post_processing(ctx)
        return
    cmd = [
        ctx.tool_paths["tdb"],
        "merge",
        "-o",
        str(output_tdb),
        "--threads",
        str(ctx.params["threads"]),
        str(input_a),
        str(input_b),
    ]
    _run_task(ctx=ctx, stage="tdb_merge", command=cmd, outputs=[output_tdb])
    _run_tr_post_processing(ctx)


def _run_tr_post_processing(ctx: RunContext) -> None:
    if len(ctx.selected_groups) != 2:
        raise ValueError("tr_post_processing requires exactly two selected groups")
    output_dir = _tr_post_stage_dir(ctx)
    summary_path = output_dir / "summary.json"
    group_a, group_b = ctx.selected_groups
    sample_a_label = _sample_name(str(ctx.params["medaka_sample_name_template"]), ctx.sample_id, group_a)
    sample_b_label = _sample_name(str(ctx.params["medaka_sample_name_template"]), ctx.sample_id, group_b)
    cmd = [
        sys.executable,
        "-m",
        "sniffcell.discover.tr_post_processing",
        "--split-dir",
        str(ctx.split_dir),
        "--groups",
        ",".join(ctx.selected_groups),
        "--output-dir",
        str(output_dir),
        "--sample-id",
        ctx.sample_id,
        "--sample-a-label",
        sample_a_label,
        "--sample-b-label",
        sample_b_label,
        "--group-a-fasta",
        str(_resolve_trimmed_reads_input(ctx, group_a)),
        "--group-b-fasta",
        str(_resolve_trimmed_reads_input(ctx, group_b)),
        "--margin-bp",
        str(ctx.params["tr_margin_bp"]),
        "--min-supporting-reads",
        str(ctx.params["tr_min_supporting_reads"]),
    ]
    _run_task(
        ctx=ctx,
        stage="tr_post_processing",
        command=cmd,
        outputs=[summary_path],
    )


def _run_snv_post_processing(ctx: RunContext) -> None:
    if len(ctx.selected_groups) != 2:
        raise ValueError("snv_post_processing requires exactly two selected groups")
    group_a, group_b = ctx.selected_groups
    group_a_gvcf = _clair3_pileup_output_path(ctx, group_a)
    group_b_gvcf = _clair3_pileup_output_path(ctx, group_b)
    if not ctx.dry_run:
        _require_paths([group_a_gvcf, group_b_gvcf])
    output_dir = _snv_post_stage_dir(ctx)
    summary_path = output_dir / "summary.json"
    cmd = [
        sys.executable,
        "-m",
        "sniffcell.discover.snv_post_processing",
        "--split-dir",
        str(ctx.split_dir),
        "--groups",
        ",".join(ctx.selected_groups),
        "--group-a-gvcf",
        str(group_a_gvcf),
        "--group-b-gvcf",
        str(group_b_gvcf),
        "--output-dir",
        str(output_dir),
        "--sample-id",
        ctx.sample_id,
        "--min-dp",
        str(DEFAULT_SNV_POST_MIN_DP),
        "--max-dp",
        str(DEFAULT_SNV_POST_MAX_DP),
        "--min-dp-absence",
        str(DEFAULT_SNV_POST_MIN_DP_ABSENCE),
        "--min-gq",
        str(DEFAULT_SNV_POST_MIN_GQ),
        "--min-other-af",
        str(DEFAULT_SNV_POST_MIN_OTHER_AF),
    ]
    _run_task(
        ctx=ctx,
        stage="snv_post_processing",
        command=cmd,
        outputs=[summary_path],
    )


def _run_modkit(ctx: RunContext, group_name: str) -> None:
    group = _group_lookup(ctx, group_name)
    stage_dir = _modkit_stage_dir(ctx, group_name)
    stage_dir.mkdir(parents=True, exist_ok=True)
    output_bed = stage_dir / f"{_sanitize_token(group_name)}.cpg.bedmethyl.gz"
    cmd = [
        ctx.tool_paths["modkit"],
        "pileup",
        group.bam_path,
        str(output_bed),
        "--cpg",
        "--reference",
        str(ctx.reference),
        "--modified-bases",
        "5mC",
        "5hmC",
        "--bgzf",
        "-t",
        str(ctx.params["threads"]),
    ]
    _run_task(ctx=ctx, stage="modkit", group_name=group_name, command=cmd, outputs=[output_bed])
    if not ctx.dry_run:
        subprocess.run(
            [ctx.tool_paths["tabix"], "-f", "-p", "bed", str(output_bed)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def _run_clair3(ctx: RunContext, group_name: str) -> None:
    group = _group_lookup(ctx, group_name)
    stage_dir = _clair3_stage_dir(ctx, group_name)
    stage_dir.mkdir(parents=True, exist_ok=True)
    output_vcf = _clair3_merge_output_path(ctx, group_name)
    pileup_vcf = _clair3_pileup_output_path(ctx, group_name)
    model_path = ctx.params.get("clair3_model_path")
    if not model_path and not ctx.dry_run:
        raise ValueError(
            "clair3 stage requires --clair3-model-path. "
            "Provide the path to a Clair3 model directory (e.g. /path/to/models/r1041_e82_400bps_sup_v500)."
        )
    cmd = [
        ctx.tool_paths["clair3"],
        f"--bam_fn={group.bam_path}",
        f"--ref_fn={str(ctx.reference)}",
        f"--output={str(stage_dir)}",
        f"--threads={ctx.params['threads']}",
        f"--platform={ctx.params['clair3_platform']}",
        f"--model_path={model_path or 'PLACEHOLDER'}",
        "--gvcf",
        # "--include_all_ctgs",
        # "--remove_intermediate_dir",
    ]
    region_bed = _deconv_region_bed(ctx)
    if region_bed is not None:
        cmd.append(f"--bed_fn={region_bed}")
    _run_task(
        ctx=ctx,
        stage="clair3",
        group_name=group_name,
        command=cmd,
        outputs=[output_vcf, pileup_vcf],
        env=_clair3_env(ctx),
    )


def _write_harmonized_variant_summary(ctx: RunContext) -> Path | None:
    if len(ctx.selected_groups) != 2:
        return None

    sv_summary = _sv_post_stage_dir(ctx) / "summary.json"
    tr_summary = _tr_post_stage_dir(ctx) / "summary.json"
    snv_summary = _snv_post_stage_dir(ctx) / "summary.json"

    sv_payload = _read_json(sv_summary) if sv_summary.exists() else {}
    tr_payload = _read_json(tr_summary) if tr_summary.exists() else {}
    snv_payload = _read_json(snv_summary) if snv_summary.exists() else {}

    sv_bed_text = str(sv_payload.get("sv_bed_tsv", "")).strip() if isinstance(sv_payload, dict) else ""
    tr_bed_text = str(tr_payload.get("tr_bed_tsv", "")).strip() if isinstance(tr_payload, dict) else ""
    snv_tsv_text = str(snv_payload.get("merged_tsv", "")).strip() if isinstance(snv_payload, dict) else ""

    sv_bed = Path(sv_bed_text) if sv_bed_text else None
    tr_bed = Path(tr_bed_text) if tr_bed_text else None
    snv_tsv = Path(snv_tsv_text) if snv_tsv_text else None
    if (sv_bed is None or not sv_bed.exists()) and (tr_bed is None or not tr_bed.exists()):
        return None

    group_a, group_b = ctx.selected_groups
    group_a_label = _sample_name(str(ctx.params["medaka_sample_name_template"]), ctx.sample_id, group_a)
    group_b_label = _sample_name(str(ctx.params["medaka_sample_name_template"]), ctx.sample_id, group_b)
    include_all_tr_candidates = True
    output_path = ctx.run_root / "harmonized_variants.tsv"
    stats = write_harmonized_variants(
        output=output_path,
        group_a_label=group_a_label,
        group_b_label=group_b_label,
        tr_bed=tr_bed,
        sv_bed=sv_bed,
        snv_tsv=snv_tsv,
        include_all_tr_candidates=include_all_tr_candidates,
    )

    _write_json(
        ctx.run_root / "harmonized_variants_manifest.json",
        {
            "created_at": _now_utc(),
            "output": str(output_path),
            "group_a_label": group_a_label,
            "group_b_label": group_b_label,
            "tr_bed_tsv": (str(tr_bed) if tr_bed is not None and tr_bed.exists() else ""),
            "sv_bed_tsv": (str(sv_bed) if sv_bed is not None and sv_bed.exists() else ""),
            "snv_tsv": (str(snv_tsv) if snv_tsv is not None and snv_tsv.exists() else ""),
            "include_all_tr_candidates": include_all_tr_candidates,
            "counts": stats,
        },
    )
    return output_path


def _write_finalize_summary(ctx: RunContext) -> None:
    harmonized_path = _write_harmonized_variant_summary(ctx)
    summary = {
        "completed_at": _now_utc(),
        "sample_id": ctx.sample_id,
        "run_id": ctx.run_id,
        "run_root": str(ctx.run_root),
        "stages": list(ctx.stages),
        "selected_groups": ctx.selected_groups,
        "scheduler": ctx.scheduler,
        "dry_run": ctx.dry_run,
        "harmonized_variants_tsv": (str(harmonized_path) if harmonized_path is not None else ""),
    }
    _write_json(ctx.run_root / "run_summary.json", summary)


def _execute_local(ctx: RunContext) -> None:
    _clear_force_rerun_state(ctx)
    task_rows: list[dict[str, str]] = []
    for group_name in ctx.selected_groups:
        task_rows.append(
            {
                "task_id": _task_id("discover", group_name),
                "stage": "discover",
                "group_name": group_name,
                "scope": "group",
                "inputs": str(ctx.split_dir / "requested_group_splits.tsv"),
                "outputs": str(Path(_group_lookup(ctx, group_name).bam_path)),
            }
        )
    _write_task_manifest(ctx, task_rows)

    stage_handlers = {
        "sniffles": _run_sniffles,
        "sniffles_filter": _run_sniffles_filter,
        "kanpig": _run_kanpig,
        "medaka": _run_medaka,
        "trgt": _run_trgt,
        "tdb_create": _run_tdb_create,
        "clair3": _run_clair3,
        "modkit": _run_modkit,
    }
    for stage in ctx.stages:
        if stage in GROUP_SCOPED_STAGES:
            for group_name in ctx.selected_groups:
                stage_handlers[stage](ctx, group_name)
            if stage == "clair3" and len(ctx.selected_groups) == 2:
                _run_snv_post_processing(ctx)
        elif stage == "collapse":
            if len(ctx.selected_groups) != 2:
                raise ValueError("collapse requires exactly two selected groups")
            _run_collapse(ctx)
        elif stage == "tdb_merge":
            if len(ctx.selected_groups) != 2:
                raise ValueError("tdb_merge requires exactly two selected groups")
            _run_tdb_merge(ctx)
        else:
            raise ValueError(f"Unsupported stage: {stage}")
    _write_finalize_summary(ctx)


def _python_entrypoint() -> list[str]:
    return [sys.executable, "-m", "sniffcell.main", "discover"]


def _build_recursive_cli(
    ctx: RunContext,
    stage_list: str,
    groups: str | None = None,
) -> list[str]:
    cmd = _python_entrypoint() + [
        "--deconv-dir",
        str(ctx.deconv_dir),
        "--reference",
        str(ctx.reference),
        "--tr-bed",
        str(ctx.tr_bed),
        "--sex",
        ctx.sex,
        "--scheduler",
        "local",
        "--run-id",
        ctx.run_id,
        "--stages",
        stage_list,
        "--sniffles-bin",
        ctx.tool_paths["sniffles"],
        "--bcftools-bin",
        ctx.tool_paths["bcftools"],
        "--bgzip-bin",
        ctx.tool_paths["bgzip"],
        "--kanpig-bin",
        ctx.tool_paths["kanpig"],
        "--truvari-bin",
        ctx.tool_paths["truvari"],
        "--medaka-bin",
        ctx.tool_paths["medaka"],
        "--tdb-bin",
        ctx.tool_paths["tdb"],
        "--modkit-bin",
        ctx.tool_paths["modkit"],
        "--tabix-bin",
        ctx.tool_paths["tabix"],
        "--threads",
        str(ctx.params["threads"]),
        "--kanpig-seqsim",
        str(ctx.params["kanpig_seqsim"]),
        "--kanpig-sizesim",
        str(ctx.params["kanpig_sizesim"]),
        "--sniffles-mosaic-filter-expression",
        str(ctx.params["mosaic_filter_expression"]),
        "--sniffles-cluster-merge-len",
        str(ctx.params["sniffles_cluster_merge_len"]),
        "--truvari-refdist",
        str(ctx.params["truvari_refdist"]),
        "--truvari-pctseq",
        str(ctx.params["truvari_pctseq"]),
        "--truvari-pctsize",
        str(ctx.params["truvari_pctsize"]),
        "--medaka-model",
        str(ctx.params["medaka_model"]),
        "--medaka-padding",
        str(ctx.params["medaka_padding"]),
        "--tdb-create-mem",
        str(ctx.params["tdb_create_mem"]),
        "--kanpig-sample-name-template",
        str(ctx.params["kanpig_sample_name_template"]),
        "--medaka-sample-name-template",
        str(ctx.params["medaka_sample_name_template"]),
        "--mods-mode",
        str(ctx.params["mods_mode"]),
        "--clair3-bin",
        ctx.tool_paths["clair3"],
        "--clair3-platform",
        str(ctx.params["clair3_platform"]),
        "--trgt-bin",
        ctx.tool_paths.get("trgt") or "trgt",
        "--samtools-bin",
        ctx.tool_paths.get("samtools") or "samtools",
        "--trgt-sample-name-template",
        str(ctx.params["trgt_sample_name_template"]),
    ]
    if ctx.params.get("trgt_karyotype"):
        cmd.extend(["--trgt-karyotype", str(ctx.params["trgt_karyotype"])])
    if ctx.params.get("medaka_phasing"):
        cmd.extend(["--medaka-phasing", str(ctx.params["medaka_phasing"])])
    if ctx.params["collapse_use"] != "kanpig":
        cmd.extend(["--collapse-use", str(ctx.params["collapse_use"])])
    if ctx.params.get("tdb_create_force"):
        cmd.append("--tdb-create-force")
    if ctx.params.get("clair3_model_path"):
        cmd.extend(["--clair3-model-path", str(ctx.params["clair3_model_path"])])
    if groups:
        cmd.extend(["--groups", groups])
    if ctx.force:
        cmd.append("--force")
    if ctx.rerun_failed:
        cmd.append("--rerun-failed")
    return cmd


def _write_slurm_script(
    *,
    script_path: Path,
    job_name: str,
    cpus: int,
    time_limit: str,
    body_lines: Sequence[str],
    array: str | None,
) -> None:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={script_path.parent.parent / 'logs' / (job_name + '_%A_%a.out' if array else job_name + '_%j.out')}",
        f"#SBATCH --error={script_path.parent.parent / 'logs' / (job_name + '_%A_%a.err' if array else job_name + '_%j.err')}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --time={time_limit}",
    ]
    if array:
        lines.append(f"#SBATCH --array={array}")
    lines.extend(["set -euo pipefail", f"cd {shlex.quote(str(REPO_ROOT))}"])
    lines.extend(body_lines)
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text("\n".join(lines) + "\n")
    script_path.chmod(0o755)


def _cmd_with_group_var(ctx: RunContext, stage_list: str) -> str:
    placeholder = "__GROUP_NAME_PLACEHOLDER__"
    cmd = _build_recursive_cli(ctx, stage_list, placeholder)
    cmd_str = shlex.join(cmd)
    return cmd_str.replace(shlex.quote(placeholder), '"${GROUP_NAME}"')


def _render_slurm(ctx: RunContext) -> None:
    ctx.slurm_dir.mkdir(parents=True, exist_ok=True)
    groups = ctx.selected_groups
    threads = int(ctx.params["threads"])
    group_manifest = ctx.slurm_dir / "groups.tsv"
    with group_manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["index", "group_name"])
        for idx, group_name in enumerate(groups, start=1):
            writer.writerow([idx, group_name])

    # Per-stage: (job_name, time_limit)
    stage_time_map = {
        "sniffles": ("pp_sniffles", "12:00:00"),
        "sniffles_filter": ("pp_snifflt", "01:00:00"),
        "kanpig": ("pp_kanpig", "12:00:00"),
        "medaka": ("pp_medaka", "24:00:00"),
        "trgt": ("pp_trgt", "06:00:00"),
        "tdb_create": ("pp_tdbc", "04:00:00"),
        "clair3": ("pp_clair3", "24:00:00"),
        "snv_post_processing": ("pp_snvpost", "04:00:00"),
        "modkit": ("pp_modkit", "12:00:00"),
        "collapse": ("pp_collapse", "02:00:00"),
        "tdb_merge": ("pp_tdbmerge", "04:00:00"),
    }

    for stage in ctx.stages:
        if stage not in stage_time_map:
            continue
        job_name, time_limit = stage_time_map[stage]
        if stage in GROUP_SCOPED_STAGES:
            script_path = ctx.slurm_dir / f"{stage}.array.sbatch.sh"
            body = [
                f"GROUP_NAME=$(awk -F '\\t' 'NR>1 && $1==ENVIRON[\"SLURM_ARRAY_TASK_ID\"] {{print $2; exit}}' {shlex.quote(str(group_manifest))})",
                f"export PYTHONPATH={shlex.quote(str(REPO_ROOT / 'src'))}:${{PYTHONPATH:-}}",
                _cmd_with_group_var(ctx, stage),
            ]
            _write_slurm_script(
                script_path=script_path,
                job_name=job_name,
                cpus=threads,
                time_limit=time_limit,
                body_lines=body,
                array=f"1-{len(groups)}",
            )
            continue
        elif stage == "collapse" and len(groups) == 2:
            script_path = ctx.slurm_dir / "collapse.sbatch.sh"
            body = [
                f"export PYTHONPATH={shlex.quote(str(REPO_ROOT / 'src'))}:${{PYTHONPATH:-}}",
                shlex.join(_build_recursive_cli(ctx, "collapse", ",".join(groups))),
            ]
            _write_slurm_script(
                script_path=script_path, job_name=job_name, cpus=threads,
                time_limit=time_limit, body_lines=body, array=None,
            )
        elif stage == "tdb_merge" and len(groups) == 2:
            script_path = ctx.slurm_dir / "tdb_merge.sbatch.sh"
            body = [
                f"export PYTHONPATH={shlex.quote(str(REPO_ROOT / 'src'))}:${{PYTHONPATH:-}}",
                shlex.join(_build_recursive_cli(ctx, "tdb_merge", ",".join(groups))),
            ]
            _write_slurm_script(
                script_path=script_path, job_name=job_name, cpus=threads,
                time_limit=time_limit, body_lines=body, array=None,
            )

    if "clair3" in ctx.stages and len(groups) == 2:
        script_path = ctx.slurm_dir / "snv_post_processing.sbatch.sh"
        snv_post_cmd = [
            sys.executable,
            "-m",
            "sniffcell.main",
            "discover",
            "ctprocessing",
            "snv",
            "--split-dir",
            str(ctx.split_dir),
            "--groups",
            ",".join(groups),
            "--group-a-gvcf",
            str(_clair3_pileup_output_path(ctx, groups[0])),
            "--group-b-gvcf",
            str(_clair3_pileup_output_path(ctx, groups[1])),
            "--output-dir",
            str(_snv_post_stage_dir(ctx)),
            "--sample-id",
            ctx.sample_id,
            "--min-dp",
            str(DEFAULT_SNV_POST_MIN_DP),
            "--max-dp",
            str(DEFAULT_SNV_POST_MAX_DP),
            "--min-dp-absence",
            str(DEFAULT_SNV_POST_MIN_DP_ABSENCE),
            "--min-gq",
            str(DEFAULT_SNV_POST_MIN_GQ),
            "--min-other-af",
            str(DEFAULT_SNV_POST_MIN_OTHER_AF),
        ]
        body = [
            f"export PYTHONPATH={shlex.quote(str(REPO_ROOT / 'src'))}:${{PYTHONPATH:-}}",
            shlex.join(snv_post_cmd),
        ]
        _write_slurm_script(
            script_path=script_path,
            job_name="pp_snvpost",
            cpus=threads,
            time_limit=stage_time_map["snv_post_processing"][1],
            body_lines=body,
            array=None,
        )

    _render_submit_script(ctx)


def _render_submit_script(ctx: RunContext) -> Path:
    """Generate a bash submission script the user edits and runs to submit all jobs."""
    slurm_dir = ctx.slurm_dir
    stages_set = set(ctx.stages)
    account_default = ctx.params.get("slurm_account") or ""

    lines = [
        "#!/bin/bash",
        "# ===================================================================",
        f"# Pipeline submission script  —  {ctx.sample_id}  (run: {ctx.run_id})",
        f"# Generated: {_now_utc()}",
        "#",
        "# Instructions:",
        "#   1. Set PARTITION below (required)",
        "#   2. Set ACCOUNT if your cluster requires it",
        "#   3. Run:  bash submit_pipeline.sh",
        "# ===================================================================",
        "",
        'PARTITION=""    # Set your cluster partition, e.g. medium, high',
        f'ACCOUNT="{account_default}"    # Leave empty to omit --account',
        "",
        "# Helper: emit --account=VALUE only when ACCOUNT is non-empty",
        '_acct() { [[ -n "$ACCOUNT" ]] && echo "--account=$ACCOUNT" || true; }',
        "",
        "set -euo pipefail",
        f"D={shlex.quote(str(slurm_dir))}",
        "",
    ]

    jid_vars: dict[str, str] = {}

    groups = ctx.selected_groups

    # First-tier: independent jobs (only need BAMs)
    first_tier_stages = [s for s in ["sniffles", "clair3", "medaka", "modkit"] if s in stages_set]
    has_first_tier = bool(first_tier_stages)
    if has_first_tier:
        lines.append("# ---- Concurrent first-tier jobs (submit together, no upstream dependency) ---")
    for stage in first_tier_stages:
        var = stage.upper() + "_JID"
        jid_vars[stage] = var
        script_name = f"{stage}.array.sbatch.sh" if stage in GROUP_SCOPED_STAGES else f"{stage}.sbatch.sh"
        lines.append(f'{var}=$(sbatch --parsable --partition="$PARTITION" $(_acct) "$D/{script_name}")')
        lines.append(f'echo "  Submitted {stage}: ${{{var}}}"')
    if has_first_tier:
        lines.append("")

    # Dependent chain
    dep_chain = []
    if "clair3" in stages_set and len(groups) == 2:
        dep_chain.append(("snv_post_processing", "clair3"))
    dep_chain.extend([
        ("sniffles_filter", "sniffles"),
        ("kanpig", "sniffles_filter"),
        ("collapse", "kanpig"),
        ("tdb_create", "medaka"),
        ("tdb_merge", "tdb_create"),
    ])
    dep_section_started = False
    for stage, dep_stage in dep_chain:
        if stage not in stages_set and not (stage == "snv_post_processing" and "clair3" in stages_set and len(groups) == 2):
            continue
        if not dep_section_started:
            lines.append("# ---- Dependent stages -------------------------------------------------------")
            dep_section_started = True
        var = stage.upper() + "_JID"
        jid_vars[stage] = var
        script_name = f"{stage}.array.sbatch.sh" if stage in GROUP_SCOPED_STAGES else f"{stage}.sbatch.sh"
        sbatch_parts = ["sbatch", "--parsable", '--partition="$PARTITION"', "$(_acct)"]
        if dep_stage in jid_vars:
            sbatch_parts.append(f"--dependency=afterok:${{{jid_vars[dep_stage]}}}")
            dep_echo = f"(after {dep_stage})"
        else:
            dep_echo = "(upstream not in this run)"
        sbatch_parts.append(f'"$D/{script_name}"')
        lines.append(f'{var}=$({" ".join(sbatch_parts)})')
        lines.append(f'echo "  Submitted {stage}: ${{{var}}} {dep_echo}"')
    if dep_section_started:
        lines.append("")

    script_path = slurm_dir / "submit_pipeline.sh"
    script_path.write_text("\n".join(lines) + "\n")
    script_path.chmod(0o755)
    return script_path


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def discover_main(args) -> None:
    _setup_logging()
    ctx = _build_context(args)
    for path in (ctx.run_root, ctx.manifest_dir, ctx.status_dir, ctx.commands_dir, ctx.logs_dir, ctx.slurm_dir):
        path.mkdir(parents=True, exist_ok=True)
    _write_run_manifest(ctx)
    if ctx.scheduler == "local":
        _execute_local(ctx)
    elif ctx.scheduler == "slurm":
        _render_slurm(ctx)
    else:  # pragma: no cover
        raise ValueError(f"Unsupported scheduler: {ctx.scheduler}")
