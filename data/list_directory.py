#!/usr/bin/env python3
"""Inventory downloaded ImmPort/GEO files and identify parser candidates.

Inputs are a project data root containing ``immport_cache`` and/or ``geo_cache``
plus a JSON classification-rules file. Outputs are deterministic TSV inventories
and a JSON provenance record. The module uses only the Python standard library,
performs no network requests, does not follow symlinks, and never invokes parsers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import platform
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

logger = logging.getLogger("downloaded_file_inventory")

INVENTORY_VERSION = "0.1.0"
INVENTORY_FILENAME = "downloaded_files.tsv"
CANDIDATES_FILENAME = "parser_candidates.tsv"
PROVENANCE_FILENAME = "file_inventory.provenance.json"
DEFAULT_DATA_ROOT = Path(__file__).resolve().parent
DEFAULT_RULES_PATH = DEFAULT_DATA_ROOT / "parser_candidate_rules.json"
DEFAULT_OUTPUT_DIRECTORY = "file_inventory"
REPOSITORY_ROOTS = {"ImmPort": "immport_cache", "GEO": "geo_cache"}
ACCESSION_PATTERNS = {
    "ImmPort": re.compile(r"^SDY\d+$", re.IGNORECASE),
    "GEO": re.compile(r"^GSE\d+$", re.IGNORECASE),
}
INVENTORY_COLUMNS = [
    "repository",
    "dataset_accession",
    "relative_path",
    "file_name",
    "size_bytes",
    "manifest_recorded",
    "artifact_type",
    "fetch_status",
    "source",
    "sha256",
    "candidate_type",
    "parser_support",
    "suggested_parser",
    "classification_reason",
]


class FileInventoryError(RuntimeError):
    """Raised when the inventory inputs or output contract are invalid."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_rules(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FileInventoryError(f"Could not read parser candidate rules {path}: {exc}") from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("rules"), list):
        raise FileInventoryError("Parser candidate rules must contain a top-level rules list")
    for index, rule in enumerate(payload["rules"]):
        required = {
            "id",
            "repositories",
            "filename_regex",
            "parser_support",
            "suggested_parser",
            "reason",
        }
        if not isinstance(rule, dict) or not required.issubset(rule):
            raise FileInventoryError(f"Parser candidate rule {index} is missing required fields")
        if rule["parser_support"] not in {"supported", "potential", "unsupported"}:
            raise FileInventoryError(
                f"Parser candidate rule {rule['id']!r} has invalid parser_support"
            )
        try:
            re.compile(str(rule["filename_regex"]))
        except re.error as exc:
            raise FileInventoryError(
                f"Parser candidate rule {rule['id']!r} has invalid regex: {exc}"
            ) from exc
    for pattern in payload.get("exclude_filename_regexes", []):
        try:
            re.compile(str(pattern))
        except re.error as exc:
            raise FileInventoryError(f"Invalid exclusion regex {pattern!r}: {exc}") from exc
    return payload


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise FileInventoryError(f"Could not read fetch manifest {path}: {exc}") from exc


def _fetch_records(cache_root: Path) -> Iterable[dict[str, Any]]:
    provenance_path = cache_root / "provenance_log.jsonl"
    if provenance_path.is_file():
        try:
            lines = provenance_path.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeError) as exc:
            raise FileInventoryError(f"Could not read {provenance_path}: {exc}") from exc
        for line_number, line in enumerate(lines, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise FileInventoryError(
                    f"Invalid JSON in {provenance_path} line {line_number}: {exc}"
                ) from exc
            if isinstance(record, dict):
                yield record
    manifest_path = cache_root / "manifest.json"
    if manifest_path.is_file():
        payload = _read_json(manifest_path)
        if not isinstance(payload, list):
            raise FileInventoryError(f"Fetch manifest must contain a list: {manifest_path}")
        yield from (record for record in payload if isinstance(record, dict))


def _artifact_key(
    repository: str,
    accession: str,
    local_file: str,
) -> tuple[str, str, str] | None:
    parts = Path(local_file).parts
    upper_parts = [part.upper() for part in parts]
    try:
        accession_index = upper_parts.index(accession.upper())
    except ValueError:
        return None
    suffix = Path(*parts[accession_index + 1 :]).as_posix()
    if not suffix:
        return None
    return repository, accession.upper(), suffix


def _manifest_index(cache_root: Path, repository: str) -> dict[tuple[str, str, str], dict]:
    index: dict[tuple[str, str, str], dict] = {}
    accession_field = "sdy_id" if repository == "ImmPort" else "gse_id"
    source_field = "source_path" if repository == "ImmPort" else "source_url"
    for record in _fetch_records(cache_root):
        accession = str(record.get(accession_field, "")).upper()
        if not ACCESSION_PATTERNS[repository].fullmatch(accession):
            continue
        for artifact in record.get("artifacts", []):
            if not isinstance(artifact, dict) or not artifact.get("local_file"):
                continue
            key = _artifact_key(repository, accession, str(artifact["local_file"]))
            if key is None:
                continue
            index[key] = {
                "artifact_type": str(artifact.get("artifact_type", "")),
                "fetch_status": str(artifact.get("status", "")),
                "source": str(artifact.get(source_field, "")),
                "sha256": str(artifact.get("sha256") or ""),
            }
    return index


def _classify(repository: str, filename: str, rules: dict[str, Any]) -> dict[str, str]:
    for rule in rules["rules"]:
        repositories = {str(value) for value in rule["repositories"]}
        if repository not in repositories and "*" not in repositories:
            continue
        if re.search(str(rule["filename_regex"]), filename):
            return {
                "candidate_type": str(rule["id"]),
                "parser_support": str(rule["parser_support"]),
                "suggested_parser": str(rule["suggested_parser"]),
                "classification_reason": str(rule["reason"]),
            }
    return {
        "candidate_type": "unclassified",
        "parser_support": "unsupported",
        "suggested_parser": "",
        "classification_reason": "No configured parser candidate rule matched.",
    }


def _is_excluded_filename(filename: str, rules: dict[str, Any]) -> bool:
    return any(
        re.search(str(pattern), filename)
        for pattern in rules.get("exclude_filename_regexes", [])
    )


def _write_tsv(records: list[dict[str, Any]], path: Path) -> None:
    try:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=INVENTORY_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise FileInventoryError(f"Could not write inventory {path}: {exc}") from exc


def _scan_repository(
    repository: str,
    cache_root: Path,
    rules: dict[str, Any],
) -> tuple[list[dict[str, Any]], Counter[str]]:
    records: list[dict[str, Any]] = []
    skipped: Counter[str] = Counter()
    if not cache_root.exists():
        return records, skipped
    if not cache_root.is_dir():
        raise FileInventoryError(f"Cache root is not a directory: {cache_root}")
    manifest = _manifest_index(cache_root, repository)
    excluded_directories = set(rules.get("exclude_directories", []))
    for directory, directory_names, filenames in os.walk(cache_root, followlinks=False):
        current = Path(directory)
        retained_directories = []
        for name in sorted(directory_names):
            path = current / name
            if name in excluded_directories or path.is_symlink():
                skipped["directory"] += 1
            else:
                retained_directories.append(name)
        directory_names[:] = retained_directories
        for filename in sorted(filenames):
            path = current / filename
            if path.is_symlink():
                skipped["symlink"] += 1
                continue
            if _is_excluded_filename(filename, rules):
                skipped["filename"] += 1
                continue
            try:
                relative = path.relative_to(cache_root)
                resolved = path.resolve(strict=True)
            except (OSError, ValueError) as exc:
                raise FileInventoryError(f"Could not resolve inventory file {path}: {exc}") from exc
            if not resolved.is_relative_to(cache_root.resolve()):
                skipped["outside_root"] += 1
                continue
            parts = relative.parts
            accession = parts[0].upper() if parts else ""
            if not ACCESSION_PATTERNS[repository].fullmatch(accession):
                accession = ""
            suffix = Path(*parts[1:]).as_posix() if accession else relative.as_posix()
            key = (repository, accession, suffix) if accession and suffix else None
            artifact = manifest.get(key, {}) if key else {}
            classification = _classify(repository, filename, rules)
            records.append(
                {
                    "repository": repository,
                    "dataset_accession": accession,
                    "relative_path": relative.as_posix(),
                    "file_name": filename,
                    "size_bytes": resolved.stat().st_size,
                    "manifest_recorded": "true" if artifact else "false",
                    "artifact_type": artifact.get("artifact_type", ""),
                    "fetch_status": artifact.get("fetch_status", ""),
                    "source": artifact.get("source", ""),
                    "sha256": artifact.get("sha256", ""),
                    **classification,
                }
            )
    return records, skipped


def build_file_inventory(
    data_root: str | Path = DEFAULT_DATA_ROOT,
    output_dir: str | Path | None = None,
    rules_path: str | Path = DEFAULT_RULES_PATH,
) -> dict[str, Any]:
    """Write a downloaded-file inventory and return its provenance record."""
    started = time.monotonic()
    root = Path(data_root).expanduser().resolve()
    rules_file = Path(rules_path).expanduser().resolve()
    destination = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else root / DEFAULT_OUTPUT_DIRECTORY
    )
    if not root.is_dir():
        raise FileInventoryError(f"Data root does not exist or is not a directory: {root}")
    rules = _load_rules(rules_file)
    records: list[dict[str, Any]] = []
    skipped: Counter[str] = Counter()
    scanned_roots = []
    for repository, directory_name in REPOSITORY_ROOTS.items():
        cache_root = root / directory_name
        repository_records, repository_skipped = _scan_repository(
            repository,
            cache_root,
            rules,
        )
        records.extend(repository_records)
        skipped.update(repository_skipped)
        scanned_roots.append(str(cache_root))
    records.sort(key=lambda item: (item["repository"], item["relative_path"]))
    candidates = [item for item in records if item["parser_support"] != "unsupported"]
    destination.mkdir(parents=True, exist_ok=True)
    inventory_path = destination / INVENTORY_FILENAME
    candidate_path = destination / CANDIDATES_FILENAME
    provenance_path = destination / PROVENANCE_FILENAME
    _write_tsv(records, inventory_path)
    _write_tsv(candidates, candidate_path)
    support_counts = Counter(item["parser_support"] for item in records)
    candidate_type_counts = Counter(item["candidate_type"] for item in candidates)
    provenance = {
        "status": "success",
        "run_started_at_utc": _utc_now(),
        "inventory_tool": "list_directory",
        "inventory_version": INVENTORY_VERSION,
        "python_version": platform.python_version(),
        "data_root": str(root),
        "scanned_roots": scanned_roots,
        "rules": {
            "path": str(rules_file),
            "version": rules.get("version"),
            "sha256": _sha256(rules_file),
        },
        "counts": {
            "files": len(records),
            "parser_candidates": len(candidates),
            "manifest_recorded": sum(
                item["manifest_recorded"] == "true" for item in records
            ),
            "parser_support": dict(sorted(support_counts.items())),
            "candidate_types": dict(sorted(candidate_type_counts.items())),
            "skipped": dict(sorted(skipped.items())),
        },
        "outputs": [
            {
                "type": "downloaded_file_inventory",
                "path": str(inventory_path),
                "rows": len(records),
                "size_bytes": inventory_path.stat().st_size,
                "sha256": _sha256(inventory_path),
            },
            {
                "type": "parser_candidates",
                "path": str(candidate_path),
                "rows": len(candidates),
                "size_bytes": candidate_path.stat().st_size,
                "sha256": _sha256(candidate_path),
            },
        ],
        "duration_sec": round(time.monotonic() - started, 3),
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    provenance["provenance_path"] = str(provenance_path)
    return provenance


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inventory downloaded files and identify configured parser candidates."
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help="Project data directory containing ImmPort/GEO caches.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory; defaults to <data-root>/file_inventory.",
    )
    parser.add_argument(
        "--rules",
        type=Path,
        default=DEFAULT_RULES_PATH,
        help="JSON parser-candidate classification rules.",
    )
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        result = build_file_inventory(args.data_root, args.output_dir, args.rules)
    except FileInventoryError as exc:
        logger.error("%s", exc)
        return 1
    logger.info(
        "Inventoried %s files and exposed %s parser candidates under %s",
        result["counts"]["files"],
        result["counts"]["parser_candidates"],
        Path(result["provenance_path"]).parent,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
