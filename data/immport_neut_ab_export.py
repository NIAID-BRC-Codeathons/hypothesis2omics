"""Export ImmPort neutralizing-antibody results as a Galaxy-ready TSV.

Inputs:
    Highest-release ImmPort Tab ZIP or extracted Tab directory for every study
    discovered below ``--cache-root``, plus the combined ImmPort sample manifest.
    Studies without a neutralizing-antibody result table are skipped and recorded
    in provenance.

Outputs:
    ``ImmPort_neut_ab_titer_results.tsv`` and a JSON provenance sidecar below
    ``--output-dir``. Each antibody result is repeated for every GEO sample from
    the same study and subject, with separate columns for the GEO sample timepoint.

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
import re
import time
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, TextIO

try:
    from data.immport_batch_parse import discover_tab_sources
except ModuleNotFoundError:
    from immport_batch_parse import discover_tab_sources  # type: ignore[no-redef]

logger = logging.getLogger("immport_neut_ab_export")

EXPORTER_VERSION = "0.3.0"
TABLE_FILENAME = "neut_ab_titer_result.txt"
DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent / "immport_cache"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "galaxy_file_input"
DEFAULT_MANIFEST_RELATIVE_PATH = Path("parsed") / "sample_manifest.tsv"
OUTPUT_FILENAME = "ImmPort_neut_ab_titer_results.tsv"
PROVENANCE_FILENAME = "ImmPort_neut_ab_titer_results.provenance.json"
REPOSITORY_TIME_COLUMNS = [
    "repository_study_time_collected",
    "repository_study_time_collected_unit",
]
REQUIRED_COLUMNS = {
    "result_id",
    "study_accession",
    "experiment_accession",
    "expsample_accession",
    "biosample_accession",
    "repository_accession",
    "repository_name",
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
MANIFEST_REQUIRED_COLUMNS = {
    "repository_accession",
    "repository_name",
    "study_accession",
    "study_time_collected",
    "study_time_collected_unit",
    "subject_accession",
}


class NeutralizingAntibodyExportError(RuntimeError):
    """Raised when neutralizing-antibody inputs violate the export contract."""


class MissingNeutralizingAntibodyTableError(NeutralizingAntibodyExportError):
    """Raised when a valid ImmPort Tab source has no neutralization table."""


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


def _read_gsm_links(
    manifest_path: Path,
) -> tuple[dict[tuple[str, str], list[dict[str, str]]], dict[str, Any]]:
    """Read GEO sample links keyed by ImmPort study and subject."""
    try:
        handle = manifest_path.open(encoding="utf-8", errors="strict", newline="")
    except OSError as exc:
        raise NeutralizingAntibodyExportError(
            f"Could not read sample manifest {manifest_path}: {exc}"
        ) from exc

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise NeutralizingAntibodyExportError(
                f"Sample manifest {manifest_path} has no header"
            )
        source_columns = list(reader.fieldnames)
        columns = [column.strip().lower() for column in source_columns]
        if len(columns) != len(set(columns)):
            raise NeutralizingAntibodyExportError(
                f"Sample manifest {manifest_path} has duplicate normalized columns"
            )
        missing = sorted(MANIFEST_REQUIRED_COLUMNS - set(columns))
        if missing:
            raise NeutralizingAntibodyExportError(
                f"Sample manifest {manifest_path} is missing required columns: "
                f"{', '.join(missing)}"
            )

        links: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        seen_accessions: dict[str, tuple[str, str]] = {}
        for line_number, source_row in enumerate(reader, start=2):
            if None in source_row:
                raise NeutralizingAntibodyExportError(
                    f"Sample manifest {manifest_path} line {line_number} has more "
                    "fields than its header"
                )
            row = {
                normalized: (source_row[source_column] or "").strip()
                for source_column, normalized in zip(
                    source_columns,
                    columns,
                    strict=True,
                )
            }
            accession = row["repository_accession"]
            if row["repository_name"].casefold() != "geo" or not re.fullmatch(
                r"(?i)GSM\d+",
                accession,
            ):
                continue
            key = (row["study_accession"].upper(), row["subject_accession"].upper())
            if not all(key):
                raise NeutralizingAntibodyExportError(
                    f"Sample manifest {manifest_path} line {line_number} has a GEO "
                    "sample without study_accession and subject_accession"
                )
            previous_key = seen_accessions.get(accession.upper())
            if previous_key is not None:
                raise NeutralizingAntibodyExportError(
                    f"Sample manifest {manifest_path} contains duplicate GEO accession "
                    f"{accession}"
                )
            seen_accessions[accession.upper()] = key
            links[key].append(
                {
                    "repository_name": "GEO",
                    "repository_accession": accession.upper(),
                    "repository_study_time_collected": row["study_time_collected"],
                    "repository_study_time_collected_unit": row[
                        "study_time_collected_unit"
                    ],
                }
            )

    for records in links.values():
        records.sort(key=lambda record: record["repository_accession"])
    return dict(links), {
        "path": str(manifest_path),
        "size_bytes": manifest_path.stat().st_size,
        "sha256": _sha256(manifest_path),
        "geo_samples": len(seen_accessions),
    }


def _repository_columns(columns: list[str]) -> list[str]:
    output_columns = list(columns)
    insertion_index = output_columns.index("repository_name") + 1
    output_columns[insertion_index:insertion_index] = REPOSITORY_TIME_COLUMNS
    return output_columns


def _link_gsm_rows(
    rows: list[dict[str, str]],
    links: dict[tuple[str, str], list[dict[str, str]]],
) -> list[dict[str, str]]:
    """Expand each antibody result to every GEO sample for its study and subject."""
    expanded: list[dict[str, str]] = []
    for row in rows:
        if row["repository_name"] or row["repository_accession"]:
            raise NeutralizingAntibodyExportError(
                "Cannot replace repository metadata already present on neutralizing-"
                f"antibody result {row['study_accession']}/{row['result_id']}"
            )
        key = (row["study_accession"].upper(), row["subject_accession"].upper())
        matches = links.get(key, [])
        if not matches:
            raise NeutralizingAntibodyExportError(
                "No GEO sample matched neutralizing-antibody result "
                f"{row['study_accession']}/{row['result_id']} for subject "
                f"{row['subject_accession']}"
            )
        for match in matches:
            linked_row = dict(row)
            linked_row.update(match)
            expanded.append(linked_row)
    return expanded


def _zip_member(archive: zipfile.ZipFile, source: Path) -> zipfile.ZipInfo:
    matches = [
        info
        for info in archive.infolist()
        if not info.is_dir()
        and PurePosixPath(info.filename).name == TABLE_FILENAME
        and "Tab" in PurePosixPath(info.filename).parts[:-1]
    ]
    if not matches:
        raise MissingNeutralizingAntibodyTableError(
            f"No {TABLE_FILENAME} member in {source}"
        )
    if len(matches) > 1:
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
            raise MissingNeutralizingAntibodyTableError(
                f"No {TABLE_FILENAME} in {source}"
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
    sample_manifest: str | Path | None = None,
) -> dict[str, Any]:
    """Export all discovered neutralization tables and return provenance."""
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    root = Path(cache_root).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    manifest_path = (
        Path(sample_manifest).expanduser().resolve()
        if sample_manifest is not None
        else root / DEFAULT_MANIFEST_RELATIVE_PATH
    )
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / OUTPUT_FILENAME
    provenance_path = destination / PROVENANCE_FILENAME

    selected = discover_tab_sources(root)
    columns: list[str] | None = None
    combined_rows: list[dict[str, str]] = []
    inputs: list[dict[str, Any]] = []
    skipped_inputs: list[dict[str, Any]] = []
    for accession, (release, source) in sorted(selected.items()):
        try:
            source_columns, rows, input_record = _read_source(accession, release, source)
        except MissingNeutralizingAntibodyTableError as exc:
            logger.info("Skipping %s: %s", accession, exc)
            skipped_inputs.append(
                {
                    "study_accession": accession,
                    "data_release": release,
                    "source_type": "directory" if source.is_dir() else "zip",
                    "path": str(source),
                    "reason": "neutralizing_antibody_table_not_present",
                }
            )
            continue
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
    gsm_links, manifest_record = _read_gsm_links(manifest_path)
    source_result_rows = len(combined_rows)
    combined_rows = _link_gsm_rows(combined_rows, gsm_links)
    columns = _repository_columns(columns)
    combined_rows.sort(
        key=lambda row: (
            row["study_accession"],
            row["expsample_accession"],
            row["repository_accession"],
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
            "gsm_linkage": "study_accession_and_subject_accession",
        },
        "inputs": inputs,
        "sample_manifest_input": manifest_record,
        "skipped_inputs": skipped_inputs,
        "counts": {
            "studies": len(inputs),
            "studies_discovered": len(selected),
            "studies_exported": len(inputs),
            "studies_skipped": len(skipped_inputs),
            "rows": len(combined_rows),
            "source_result_rows": source_result_rows,
            "linked_result_rows": len(combined_rows),
            "distinct_gsm_accessions": len(
                {row["repository_accession"] for row in combined_rows}
            ),
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
    parser.add_argument(
        "--sample-manifest",
        help=(
            "Combined sample manifest used to link results to all GEO samples; "
            "defaults to <cache-root>/parsed/sample_manifest.tsv"
        ),
    )
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        provenance = run_export(
            args.cache_root,
            args.output_dir,
            args.sample_manifest,
        )
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
