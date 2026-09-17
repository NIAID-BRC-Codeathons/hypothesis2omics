"""Export ImmPort neutralizing-antibody results as a Galaxy-ready TSV.

Inputs:
    Highest-release ImmPort Tab ZIP or extracted Tab directory for every study
    discovered below ``--cache-root``.

Outputs:
    ``ImmPort_neut_ab_titer_results.tsv`` and a JSON provenance sidecar below
    ``--output-dir``. Source columns are normalized to lowercase and values are
    preserved as submitted after surrounding whitespace is removed.

External dependencies:
    The local ``immport_batch_parse`` source-discovery module; no network access
    or credentials are required.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import logging
import os
import platform
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, TextIO

try:
    from data.immport_batch_parse import discover_tab_sources
except ModuleNotFoundError:
    from immport_batch_parse import discover_tab_sources  # type: ignore[no-redef]

logger = logging.getLogger("immport_neut_ab_export")

EXPORTER_VERSION = "0.1.0"
TABLE_FILENAME = "neut_ab_titer_result.txt"
DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent / "immport_cache"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "galaxy_file_input"
OUTPUT_FILENAME = "ImmPort_neut_ab_titer_results.tsv"
PROVENANCE_FILENAME = "ImmPort_neut_ab_titer_results.provenance.json"
REQUIRED_COLUMNS = {
    "result_id",
    "study_accession",
    "experiment_accession",
    "expsample_accession",
    "biosample_accession",
    "subject_accession",
    "study_time_collected",
    "study_time_collected_unit",
    "value_preferred",
    "value_reported",
    "unit_preferred",
    "unit_reported",
    "virus_strain_preferred",
    "virus_strain_reported",
}


class NeutralizingAntibodyExportError(RuntimeError):
    """Raised when neutralizing-antibody inputs violate the export contract."""


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalize_columns(columns: list[str], source: str) -> list[str]:
    normalized = [column.strip().lower() for column in columns]
    if any(not column for column in normalized):
        raise NeutralizingAntibodyExportError(f"{source} has an empty column name")
    if len(normalized) != len(set(normalized)):
        raise NeutralizingAntibodyExportError(
            f"{source} has duplicate columns after normalization"
        )
    missing = sorted(REQUIRED_COLUMNS - set(normalized))
    if missing:
        raise NeutralizingAntibodyExportError(
            f"{source} is missing required columns: {', '.join(missing)}"
        )
    return normalized


def _read_rows(handle: TextIO, source: str) -> tuple[list[str], list[dict[str, str]]]:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise NeutralizingAntibodyExportError(f"{source} has no header")
    source_columns = list(reader.fieldnames)
    columns = _normalize_columns(source_columns, source)
    rows: list[dict[str, str]] = []
    for line_number, source_row in enumerate(reader, start=2):
        if None in source_row:
            raise NeutralizingAntibodyExportError(
                f"{source} line {line_number} has more fields than its header"
            )
        row = {
            normalized: (source_row[source_column] or "").strip()
            for source_column, normalized in zip(source_columns, columns, strict=True)
        }
        rows.append(row)
    if not rows:
        raise NeutralizingAntibodyExportError(f"{source} has no data rows")
    return columns, rows


def _zip_member(archive: zipfile.ZipFile, source: Path) -> zipfile.ZipInfo:
    matches = [
        info
        for info in archive.infolist()
        if not info.is_dir()
        and PurePosixPath(info.filename).name == TABLE_FILENAME
        and "Tab" in PurePosixPath(info.filename).parts[:-1]
    ]
    if len(matches) != 1:
        raise NeutralizingAntibodyExportError(
            f"Expected one {TABLE_FILENAME} member in {source}; found {len(matches)}"
        )
    return matches[0]


def _read_source(
    accession: str,
    release: int,
    source: Path,
) -> tuple[list[str], list[dict[str, str]], dict[str, Any]]:
    if source.is_dir():
        table_path = source / TABLE_FILENAME
        if not table_path.is_file():
            raise NeutralizingAntibodyExportError(
                f"Required neutralization table is missing: {table_path}"
            )
        with table_path.open(encoding="utf-8", errors="strict", newline="") as handle:
            columns, rows = _read_rows(handle, str(table_path))
        input_record = {
            "study_accession": accession,
            "data_release": release,
            "source_type": "directory_table",
            "path": str(table_path),
            "size_bytes": table_path.stat().st_size,
            "sha256": _sha256(table_path),
            "rows": len(rows),
        }
    else:
        try:
            with zipfile.ZipFile(source) as archive:
                member = _zip_member(archive, source)
                with archive.open(member) as raw_handle:
                    with io.TextIOWrapper(
                        raw_handle,
                        encoding="utf-8",
                        errors="strict",
                        newline="",
                    ) as handle:
                        columns, rows = _read_rows(handle, f"{source}!{member.filename}")
                member_name = member.filename
                member_size = member.file_size
                member_crc32 = f"{member.CRC:08x}"
        except (OSError, UnicodeError, zipfile.BadZipFile) as exc:
            raise NeutralizingAntibodyExportError(
                f"Could not read ImmPort Tab ZIP {source}: {exc}"
            ) from exc
        input_record = {
            "study_accession": accession,
            "data_release": release,
            "source_type": "zip_member",
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "sha256": _sha256(source),
            "member": member_name,
            "member_size_bytes": member_size,
            "member_crc32": member_crc32,
            "rows": len(rows),
        }

    for row_number, row in enumerate(rows, start=2):
        if row["study_accession"].upper() != accession:
            raise NeutralizingAntibodyExportError(
                f"{accession} row {row_number} has study accession "
                f"{row['study_accession']!r}"
            )
        for column in ("result_id", "expsample_accession", "biosample_accession"):
            if not row[column]:
                raise NeutralizingAntibodyExportError(
                    f"{accession} row {row_number} has an empty {column}"
                )
    return columns, rows, input_record


def _write_tsv(columns: list[str], rows: list[dict[str, str]], output: Path) -> None:
    temporary = output.with_name(f"{output.name}.part")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()


def run_export(
    cache_root: str | Path = DEFAULT_CACHE_ROOT,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Export all discovered neutralization tables and return provenance."""
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    root = Path(cache_root).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / OUTPUT_FILENAME
    provenance_path = destination / PROVENANCE_FILENAME

    selected = discover_tab_sources(root)
    columns: list[str] | None = None
    combined_rows: list[dict[str, str]] = []
    inputs: list[dict[str, Any]] = []
    for accession, (release, source) in sorted(selected.items()):
        source_columns, rows, input_record = _read_source(accession, release, source)
        if columns is None:
            columns = source_columns
        elif source_columns != columns:
            raise NeutralizingAntibodyExportError(
                f"{source} columns do not match the first neutralization table"
            )
        combined_rows.extend(rows)
        inputs.append(input_record)

    if columns is None or not combined_rows:
        raise NeutralizingAntibodyExportError(
            f"No neutralizing-antibody results were found under {root}"
        )

    keys = [(row["study_accession"], row["result_id"]) for row in combined_rows]
    if len(keys) != len(set(keys)):
        raise NeutralizingAntibodyExportError(
            "Neutralizing-antibody results contain duplicate study/result identifiers"
        )
    combined_rows.sort(
        key=lambda row: (
            row["study_accession"],
            row["expsample_accession"],
            row["result_id"],
        )
    )
    _write_tsv(columns, combined_rows, output)

    preferred_count = sum(bool(row["value_preferred"]) for row in combined_rows)
    provenance = {
        "status": "success",
        "run_started_at_utc": started_at,
        "exporter": "immport_neut_ab_export",
        "exporter_version": EXPORTER_VERSION,
        "python_version": platform.python_version(),
        "selection": {
            "cache_root": str(root),
            "mode": "highest_release_tab_source_per_study",
            "table": TABLE_FILENAME,
        },
        "inputs": inputs,
        "counts": {
            "studies": len(inputs),
            "rows": len(combined_rows),
            "value_preferred_populated": preferred_count,
            "value_preferred_blank": len(combined_rows) - preferred_count,
        },
        "output": {
            "path": str(output),
            "rows": len(combined_rows),
            "columns": columns,
            "size_bytes": output.stat().st_size,
            "sha256": _sha256(output),
        },
        "duration_sec": round(time.monotonic() - started, 3),
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return provenance


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export ImmPort neutralizing-antibody results for Galaxy."
    )
    parser.add_argument(
        "--cache-root",
        default=str(DEFAULT_CACHE_ROOT),
        help="Root containing ImmPort study directories",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for the Galaxy TSV and provenance JSON",
    )
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        provenance = run_export(args.cache_root, args.output_dir)
    except (NeutralizingAntibodyExportError, OSError, ValueError) as exc:
        logger.error("Neutralizing-antibody export failed: %s", exc)
        return 1
    logger.info(
        "Wrote %s rows from %s studies to %s",
        provenance["counts"]["rows"],
        provenance["counts"]["studies"],
        provenance["output"]["path"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
