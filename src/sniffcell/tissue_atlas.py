from __future__ import annotations

import json
import re
from importlib import resources
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


def default_tissue_atlas_path() -> str:
    return str(resources.files("sniffcell.data").joinpath("tissue_atlas.json"))


def normalize_tissue_token(value: Any) -> str:
    text = "" if pd.isna(value) else str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def atlas_tissue_records(atlas_payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = atlas_payload.get("__tissues__", [])
    if not isinstance(rows, list):
        return []
    return [row for row in rows if isinstance(row, dict)]


def resolve_tissue_key(celltypes_key: str, atlas_payload: Mapping[str, Any]) -> tuple[str, dict[str, Any] | None]:
    query = str(celltypes_key).strip()
    if query in atlas_payload and isinstance(atlas_payload[query], Mapping):
        return query, None

    query_norm = normalize_tissue_token(query)
    matches: list[dict[str, Any]] = []
    for row in atlas_tissue_records(atlas_payload):
        aliases = [
            row.get("code", ""),
            row.get("name", ""),
            row.get("key", ""),
        ]
        if any(normalize_tissue_token(alias) == query_norm for alias in aliases if str(alias).strip()):
            matches.append(row)

    if not matches:
        partial = [
            row
            for row in atlas_tissue_records(atlas_payload)
            if query_norm
            and (
                query_norm in normalize_tissue_token(row.get("name", ""))
                or query_norm in normalize_tissue_token(row.get("code", ""))
            )
        ]
        if len(partial) == 1:
            matches = partial
        elif len(partial) > 1:
            choices = ", ".join(f"{row.get('code')} ({row.get('name')})" for row in partial)
            raise ValueError(f"Tissue key '{query}' is ambiguous. Matches: {choices}")

    if not matches:
        available = ", ".join(
            f"{row.get('code')} ({row.get('name')})"
            for row in atlas_tissue_records(atlas_payload)
            if row.get("usable") and row.get("key")
        )
        raise KeyError(f"Tissue key '{query}' not found. Available tissue keys: {available}")

    if len(matches) > 1:
        choices = ", ".join(f"{row.get('code')} ({row.get('name')})" for row in matches)
        raise ValueError(f"Tissue key '{query}' is ambiguous. Matches: {choices}")

    row = matches[0]
    if not row.get("usable") or not row.get("key"):
        raise ValueError(
            f"Tissue '{row.get('code')} ({row.get('name')})' is listed in the atlas metadata "
            "but does not have enough reference cell types for sniffcell-lite find."
        )
    key = str(row["key"])
    if key not in atlas_payload or not isinstance(atlas_payload[key], Mapping):
        raise KeyError(f"Tissue '{row.get('code')} ({row.get('name')})' points to missing atlas key '{key}'.")
    return key, row


def write_tissue_catalog_manifest(
    *,
    output_path: Path,
    celltypes_file: str,
    requested_key: str,
    resolved_key: str,
    tissue_row: Mapping[str, Any] | None,
) -> None:
    manifest_path = output_path.with_suffix(output_path.suffix + ".catalog.json")
    payload: dict[str, Any] = {
        "catalog": str(output_path.resolve()),
        "igv_bed": str(output_path.with_suffix(output_path.suffix + ".igv.bed").resolve()),
        "celltypes_file": str(Path(celltypes_file).resolve()),
        "requested_key": requested_key,
        "resolved_key": resolved_key,
    }
    if tissue_row is not None:
        payload["tissue"] = {
            "code": tissue_row.get("code", ""),
            "name": tissue_row.get("name", ""),
            "confidence": tissue_row.get("confidence", ""),
            "notes": tissue_row.get("notes", ""),
        }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
