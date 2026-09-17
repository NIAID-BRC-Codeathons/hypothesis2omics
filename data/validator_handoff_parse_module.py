"""Create bounded GEO and ImmPort tables for scientific-validator handoff.

Inputs are a JSON configuration, one normalized ImmPort sample manifest, one
cached GEO series matrix, and one ImmPort Tab ZIP. Outputs are a validator-input
bundle containing three TSV files, per-file provenance, and a bundle manifest.
The configuration explicitly identifies the GEO feature and ImmPort outcome;
this module does not infer scientific eligibility or normalize measurements.

External dependencies: none beyond the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import logging
import platform
import re
import shutil
import time
import zipfile
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath
from typing import Any

logger = logging.getLogger("validator_handoff_parser")

PARSER_VERSION = "0.3.0"
SAMPLE_MANIFEST_FILENAME = "sample_manifest.tsv"
FEATURE_EXPRESSION_FILENAME = "feature_expression.tsv"
QUANTITATIVE_OUTCOME_FILENAME = "quantitative_outcome.tsv"
BUNDLE_MANIFEST_FILENAME = "validator_input_manifest.json"
GSM_PATTERN = re.compile(r"^GSM\d+$", re.IGNORECASE)
GSE_PATTERN = re.compile(r"^GSE\d+$", re.IGNORECASE)
GPL_PATTERN = re.compile(r"^GPL\d+$", re.IGNORECASE)
SDY_PATTERN = re.compile(r"^SDY\d+$", re.IGNORECASE)

GEO_COLUMNS = [
    "study_accession",
    "gse_accession",
    "gpl_accession",
    "gene",
    "feature_id",
    "sample_accession",
    "expression_value",
]
IMMPORT_COLUMNS = [
    "study_accession",
    "subject_accession",
    "biosample_accession",
    "timepoint_day",
    "outcome_name",
    "result_value",
    "result_unit",
]


class ValidatorHandoffParseError(RuntimeError):
    """Raised when configured source data violate the export contract."""


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require_text(mapping: dict[str, Any], key: str, section: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValidatorHandoffParseError(
            f"Configuration {section}.{key} must be a non-empty string"
        )
    return value.strip()


def _validate_accession(value: str, pattern: re.Pattern[str], label: str) -> str:
    normalized = value.upper()
    if not pattern.fullmatch(normalized):
        raise ValidatorHandoffParseError(f"Invalid {label}: {value!r}")
    return normalized


def _validate_numeric(value: str, label: str) -> None:
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValidatorHandoffParseError(f"{label} is not numeric: {value!r}") from exc
    if not number.is_finite():
        raise ValidatorHandoffParseError(f"{label} is not finite: {value!r}")


def _read_geo_feature(config: dict[str, Any], source: Path) -> list[dict[str, str]]:
    study = _validate_accession(
        _require_text(config, "study_accession", "geo"), SDY_PATTERN, "study accession"
    )
    expected_gse = _validate_accession(
        _require_text(config, "gse_accession", "geo"), GSE_PATTERN, "GSE accession"
    )
    expected_gpl = _validate_accession(
        _require_text(config, "gpl_accession", "geo"), GPL_PATTERN, "GPL accession"
    )
    gene = _require_text(config, "gene", "geo")
    feature = _require_text(config, "feature_id", "geo")

    gse = ""
    samples: list[str] = []
    platforms: list[str] = []
    table_found = False
    table_lines: list[str] = []
    try:
        with gzip.open(source, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            for line in handle:
                fields = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
                if fields and fields[0] == "!series_matrix_table_begin":
                    table_found = True
                    break
                if fields and fields[0] == "!Series_geo_accession" and len(fields) > 1:
                    gse = fields[1].strip().upper()
                elif fields and fields[0] == "!Sample_geo_accession":
                    samples = [value.strip().upper() for value in fields[1:]]
                elif fields and fields[0] == "!Sample_platform_id":
                    platforms = [value.strip().upper() for value in fields[1:]]
            if table_found:
                table_lines = list(handle)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValidatorHandoffParseError(f"Could not read GEO matrix {source}: {exc}") from exc

    if not table_found:
        raise ValidatorHandoffParseError(f"GEO matrix has no table-begin marker: {source}")
    if gse != expected_gse:
        raise ValidatorHandoffParseError(
            f"GEO matrix series {gse!r} does not match configured {expected_gse}"
        )
    if not samples or any(not GSM_PATTERN.fullmatch(value) for value in samples):
        raise ValidatorHandoffParseError("GEO matrix has invalid or missing sample accessions")
    if len(samples) != len(set(samples)):
        raise ValidatorHandoffParseError("GEO matrix has duplicate sample accessions")
    if len(platforms) != len(samples) or set(platforms) != {expected_gpl}:
        raise ValidatorHandoffParseError(
            f"GEO sample platforms do not all match configured {expected_gpl}"
        )

    reader = csv.reader(table_lines, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration as exc:
        raise ValidatorHandoffParseError("GEO expression table is empty") from exc
    if header != ["ID_REF", *samples]:
        raise ValidatorHandoffParseError("GEO expression header does not match sample metadata")

    matches: list[list[str]] = []
    end_markers = 0
    for fields in reader:
        if fields == ["!series_matrix_table_end"]:
            end_markers += 1
            continue
        if fields and fields[0] == feature:
            matches.append(fields)
    if end_markers != 1:
        raise ValidatorHandoffParseError("GEO expression table must have one end marker")
    if len(matches) != 1:
        raise ValidatorHandoffParseError(
            f"Expected one GEO row for feature {feature!r}; found {len(matches)}"
        )
    values = matches[0][1:]
    if len(values) != len(samples):
        raise ValidatorHandoffParseError(f"GEO feature {feature!r} has the wrong value count")

    records = []
    for sample, value in zip(samples, values, strict=True):
        value = value.strip()
        _validate_numeric(value, f"Expression value for {sample}")
        records.append(
            {
                "study_accession": study,
                "gse_accession": expected_gse,
                "gpl_accession": expected_gpl,
                "gene": gene,
                "feature_id": feature,
                "sample_accession": sample,
                "expression_value": value,
            }
        )
    return records


def _zip_member(archive: zipfile.ZipFile, filename: str) -> zipfile.ZipInfo:
    matches = [
        info
        for info in archive.infolist()
        if not info.is_dir()
        and PurePosixPath(info.filename).name == filename
        and "Tab" in PurePosixPath(info.filename).parts[:-1]
    ]
    if len(matches) != 1:
        raise ValidatorHandoffParseError(
            f"Expected one {filename} member in {archive.filename}; found {len(matches)}"
        )
    return matches[0]


def _read_zip_table(archive: zipfile.ZipFile, filename: str) -> list[dict[str, str]]:
    member = _zip_member(archive, filename)
    with archive.open(member) as raw_handle:
        with io.TextIOWrapper(raw_handle, encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValidatorHandoffParseError(f"ImmPort table {filename} has no header")
            normalized = [field.strip().upper() for field in reader.fieldnames]
            if len(normalized) != len(set(normalized)):
                raise ValidatorHandoffParseError(
                    f"ImmPort table {filename} has duplicate normalized columns"
                )
            records = []
            for row in reader:
                records.append(
                    {
                        key.strip().upper(): (value or "").strip()
                        for key, value in row.items()
                    }
                )
    return records


def _require_columns(records: list[dict[str, str]], columns: set[str], label: str) -> None:
    if not records:
        raise ValidatorHandoffParseError(f"ImmPort table {label} has no data rows")
    missing = sorted(columns - set(records[0]))
    if missing:
        raise ValidatorHandoffParseError(
            f"ImmPort table {label} is missing columns: {', '.join(missing)}"
        )


def _preferred(row: dict[str, str], preferred: str, reported: str) -> str:
    return row[preferred] or row[reported]


def _read_immport_outcomes(
    config: dict[str, Any], source: Path
) -> list[dict[str, str]]:
    study = _validate_accession(
        _require_text(config, "study_accession", "immport"),
        SDY_PATTERN,
        "study accession",
    )
    names = config.get("outcome_names")
    if not isinstance(names, list) or not names or any(
        not isinstance(name, str) or not name.strip() for name in names
    ):
        raise ValidatorHandoffParseError(
            "Configuration immport.outcome_names must be a non-empty string list"
        )
    selected_names = {name.strip().casefold() for name in names}
    expected_unit = _require_text(config, "timepoint_source_unit", "immport")
    precedence = config.get("field_precedence")
    if precedence != "preferred_then_reported":
        raise ValidatorHandoffParseError(
            "Configuration immport.field_precedence must be preferred_then_reported"
        )

    try:
        with zipfile.ZipFile(source) as archive:
            lab_rows = _read_zip_table(archive, "lab_test.txt")
            biosample_rows = _read_zip_table(archive, "biosample.txt")
    except (OSError, UnicodeError, zipfile.BadZipFile) as exc:
        raise ValidatorHandoffParseError(f"Could not read ImmPort ZIP {source}: {exc}") from exc

    _require_columns(
        lab_rows,
        {
            "BIOSAMPLE_ACCESSION",
            "NAME_PREFERRED",
            "NAME_REPORTED",
            "RESULT_UNIT_PREFERRED",
            "RESULT_UNIT_REPORTED",
            "RESULT_VALUE_PREFERRED",
            "RESULT_VALUE_REPORTED",
        },
        "lab_test.txt",
    )
    _require_columns(
        biosample_rows,
        {
            "BIOSAMPLE_ACCESSION",
            "STUDY_ACCESSION",
            "STUDY_TIME_COLLECTED",
            "STUDY_TIME_COLLECTED_UNIT",
            "SUBJECT_ACCESSION",
        },
        "biosample.txt",
    )

    biosamples: dict[str, dict[str, str]] = {}
    for row in biosample_rows:
        accession = row["BIOSAMPLE_ACCESSION"]
        if not accession or accession in biosamples:
            raise ValidatorHandoffParseError(
                f"biosample.txt has an empty or duplicate accession: {accession!r}"
            )
        biosamples[accession] = row

    records = []
    for lab_row in lab_rows:
        outcome_name = _preferred(lab_row, "NAME_PREFERRED", "NAME_REPORTED")
        if outcome_name.casefold() not in selected_names:
            continue
        biosample_accession = lab_row["BIOSAMPLE_ACCESSION"]
        if biosample_accession not in biosamples:
            raise ValidatorHandoffParseError(
                f"Outcome references missing biosample {biosample_accession!r}"
            )
        biosample = biosamples[biosample_accession]
        if biosample["STUDY_ACCESSION"].upper() != study:
            raise ValidatorHandoffParseError(
                f"Biosample {biosample_accession} does not belong to {study}"
            )
        if biosample["STUDY_TIME_COLLECTED_UNIT"].casefold() != expected_unit.casefold():
            raise ValidatorHandoffParseError(
                f"Biosample {biosample_accession} timepoint unit is not {expected_unit!r}"
            )
        subject = biosample["SUBJECT_ACCESSION"]
        timepoint = biosample["STUDY_TIME_COLLECTED"]
        value = _preferred(lab_row, "RESULT_VALUE_PREFERRED", "RESULT_VALUE_REPORTED")
        unit = _preferred(lab_row, "RESULT_UNIT_PREFERRED", "RESULT_UNIT_REPORTED")
        if not subject or not timepoint or not outcome_name or not unit:
            raise ValidatorHandoffParseError(
                f"Outcome for {biosample_accession} has an empty required field"
            )
        _validate_numeric(timepoint, f"Timepoint for {biosample_accession}")
        _validate_numeric(value, f"Outcome value for {biosample_accession}")
        records.append(
            {
                "study_accession": study,
                "subject_accession": subject,
                "biosample_accession": biosample_accession,
                "timepoint_day": timepoint,
                "outcome_name": outcome_name,
                "result_value": value,
                "result_unit": unit,
            }
        )
    if not records:
        raise ValidatorHandoffParseError("No ImmPort rows matched the configured outcome names")
    return records


def _write_tsv(records: list[dict[str, str]], columns: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _inspect_tsv(path: Path) -> tuple[int, list[str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader)
            rows = sum(1 for _ in reader)
    except StopIteration as exc:
        raise ValidatorHandoffParseError(f"Sample manifest is empty: {path}") from exc
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValidatorHandoffParseError(f"Could not read sample manifest {path}: {exc}") from exc
    columns = [column.strip() for column in header]
    if not columns or any(not column for column in columns):
        raise ValidatorHandoffParseError("Sample manifest has an empty column name")
    if len(columns) != len(set(columns)):
        raise ValidatorHandoffParseError("Sample manifest has duplicate columns")
    return rows, columns


def _provenance(
    source_label: str,
    source_path_text: str,
    source: Path,
    output_path_text: str,
    output: Path,
    rows: int,
    columns: list[str],
    selection: dict[str, Any],
    started_at: str,
    duration_sec: float,
) -> dict[str, Any]:
    return {
        "status": "success",
        "run_started_at_utc": started_at,
        "parser": "validator_handoff_parse_module",
        "parser_version": PARSER_VERSION,
        "python_version": platform.python_version(),
        "input": {
            "type": source_label,
            "path": source_path_text,
            "size_bytes": source.stat().st_size,
            "sha256": _sha256(source),
        },
        "selection": selection,
        "output": {
            "path": output_path_text,
            "rows": rows,
            "columns": columns,
            "size_bytes": output.stat().st_size,
            "sha256": _sha256(output),
        },
        "duration_sec": round(duration_sec, 3),
    }

def _configured_dataset_jobs(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Return validator handoff jobs in a backward-compatible form.

    Supported configuration styles:

    Legacy single-dataset mode::

        {
          "manifest": {...},
          "geo": {...},
          "immport": {...}
        }

    Multi-dataset mode::

        {
          "manifest": {...},
          "datasets": [
            {"name": "...", "geo": {...}, "immport": {...}},
            ...
          ]
        }

    Each dataset entry may contain GEO, ImmPort, or both. At least one section
    must be present. Exact duplicate output records are removed when jobs share
    an ImmPort study/outcome or otherwise overlap.
    """
    datasets = config.get("datasets")
    legacy_geo = config.get("geo")
    legacy_immport = config.get("immport")

    if datasets is None:
        if not isinstance(legacy_geo, dict) or not isinstance(legacy_immport, dict):
            raise ValidatorHandoffParseError(
                "Configuration requires either a non-empty datasets list or legacy geo and immport objects"
            )
        return [
            {
                "name": "legacy_single_dataset",
                "geo": legacy_geo,
                "immport": legacy_immport,
            }
        ]

    if legacy_geo is not None or legacy_immport is not None:
        raise ValidatorHandoffParseError(
            "Configuration must use either datasets or top-level geo/immport, not both"
        )
    if not isinstance(datasets, list) or not datasets:
        raise ValidatorHandoffParseError(
            "Configuration datasets must be a non-empty list"
        )

    jobs: list[dict[str, Any]] = []
    for index, dataset in enumerate(datasets):
        section = f"datasets[{index}]"
        if not isinstance(dataset, dict):
            raise ValidatorHandoffParseError(f"Configuration {section} must be an object")
        geo = dataset.get("geo")
        immport = dataset.get("immport")
        if geo is not None and not isinstance(geo, dict):
            raise ValidatorHandoffParseError(f"Configuration {section}.geo must be an object")
        if immport is not None and not isinstance(immport, dict):
            raise ValidatorHandoffParseError(f"Configuration {section}.immport must be an object")
        if geo is None and immport is None:
            raise ValidatorHandoffParseError(
                f"Configuration {section} must contain geo, immport, or both"
            )
        name = dataset.get("name")
        if name is None:
            name = f"dataset_{index + 1}"
        elif not isinstance(name, str) or not name.strip():
            raise ValidatorHandoffParseError(
                f"Configuration {section}.name must be a non-empty string when provided"
            )
        jobs.append({"name": str(name).strip(), "geo": geo, "immport": immport})
    return jobs


def _deduplicate_records(
    records: list[dict[str, str]], columns: list[str]
) -> list[dict[str, str]]:
    """Remove exact duplicate output rows while preserving first-seen order."""
    seen: set[tuple[str, ...]] = set()
    unique: list[dict[str, str]] = []
    for record in records:
        key = tuple(record[column] for column in columns)
        if key in seen:
            continue
        seen.add(key)
        unique.append(record)
    return unique


def _combined_provenance(
    source_label: str,
    sources: list[dict[str, Any]],
    output_path_text: str,
    output: Path,
    rows: int,
    columns: list[str],
    started_at: str,
    duration_sec: float,
) -> dict[str, Any]:
    """Create provenance for one or many configured source files."""
    if len(sources) == 1:
        item = sources[0]
        return _provenance(
            source_label,
            item["path_text"],
            item["path"],
            output_path_text,
            output,
            rows,
            columns,
            item["selection"],
            started_at,
            duration_sec,
        )

    return {
        "status": "success",
        "run_started_at_utc": started_at,
        "parser": "validator_handoff_parse_module",
        "parser_version": PARSER_VERSION,
        "python_version": platform.python_version(),
        "inputs": [
            {
                "type": source_label,
                "dataset_name": item["dataset_name"],
                "path": item["path_text"],
                "size_bytes": item["path"].stat().st_size,
                "sha256": _sha256(item["path"]),
                "selection": item["selection"],
            }
            for item in sources
        ],
        "output": {
            "path": output_path_text,
            "rows": rows,
            "columns": columns,
            "size_bytes": output.stat().st_size,
            "sha256": _sha256(output),
        },
        "duration_sec": round(duration_sec, 3),
    }



def run_parser(config_path: str | Path) -> dict[str, Any]:
    """Write one aggregated validator-input bundle and return provenance records.

    The output contract remains the same regardless of whether the configuration
    contains one dataset or many: sample_manifest.tsv, feature_expression.tsv,
    quantitative_outcome.tsv, their provenance records, and a bundle manifest.
    """
    run_started_at = _utc_now()
    config_source = Path(config_path)
    try:
        config = json.loads(config_source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValidatorHandoffParseError(f"Could not read configuration: {exc}") from exc
    if not isinstance(config, dict):
        raise ValidatorHandoffParseError("Configuration root must be an object")

    manifest = config.get("manifest")
    if not isinstance(manifest, dict):
        raise ValidatorHandoffParseError("Configuration requires a manifest object")
    jobs = _configured_dataset_jobs(config)

    output_dir_text = _require_text(config, "output_dir", "root")
    output_dir = Path(output_dir_text)

    # The normalized sample manifest is already study-aware, so it is copied once
    # and can contain rows for every study represented in the configured jobs.
    manifest_source_text = _require_text(manifest, "source", "manifest")
    manifest_source = Path(manifest_source_text)
    manifest_output = output_dir / SAMPLE_MANIFEST_FILENAME
    manifest_provenance_path = output_dir / "sample_manifest.provenance.json"
    manifest_started_at = _utc_now()
    started = time.monotonic()
    manifest_rows, manifest_columns = _inspect_tsv(manifest_source)
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    if manifest_source.resolve() != manifest_output.resolve():
        shutil.copyfile(manifest_source, manifest_output)
    manifest_provenance = _provenance(
        "normalized_immport_sample_manifest",
        manifest_source_text,
        manifest_source,
        str(manifest_output),
        manifest_output,
        manifest_rows,
        manifest_columns,
        {"mode": "verbatim_copy"},
        manifest_started_at,
        time.monotonic() - started,
    )
    _write_json(manifest_provenance, manifest_provenance_path)

    # Aggregate feature-expression rows across every configured GEO job.
    geo_output = output_dir / FEATURE_EXPRESSION_FILENAME
    geo_provenance_path = output_dir / "feature_expression.provenance.json"
    geo_started_at = _utc_now()
    started = time.monotonic()
    geo_records: list[dict[str, str]] = []
    geo_sources: list[dict[str, Any]] = []
    for job in jobs:
        geo = job.get("geo")
        if geo is None:
            continue
        geo_source_text = _require_text(geo, "source_matrix", f"datasets[{job['name']}].geo")
        geo_source = Path(geo_source_text)
        geo_records.extend(_read_geo_feature(geo, geo_source))
        geo_sources.append(
            {
                "dataset_name": job["name"],
                "path_text": geo_source_text,
                "path": geo_source,
                "selection": {
                    key: geo[key]
                    for key in [
                        "study_accession",
                        "gse_accession",
                        "gpl_accession",
                        "gene",
                        "feature_id",
                    ]
                },
            }
        )
    if not geo_sources:
        raise ValidatorHandoffParseError("No GEO jobs were configured")
    geo_records = _deduplicate_records(geo_records, GEO_COLUMNS)
    _write_tsv(geo_records, GEO_COLUMNS, geo_output)
    geo_provenance = _combined_provenance(
        "geo_series_matrix",
        geo_sources,
        str(geo_output),
        geo_output,
        len(geo_records),
        GEO_COLUMNS,
        geo_started_at,
        time.monotonic() - started,
    )
    _write_json(geo_provenance, geo_provenance_path)

    # Aggregate quantitative outcomes across every configured ImmPort job.
    # Exact duplicates are removed so the same study/outcome can safely be
    # referenced by more than one GEO dataset entry.
    immport_output = output_dir / QUANTITATIVE_OUTCOME_FILENAME
    immport_provenance_path = output_dir / "quantitative_outcome.provenance.json"
    immport_started_at = _utc_now()
    started = time.monotonic()
    immport_records: list[dict[str, str]] = []
    immport_sources: list[dict[str, Any]] = []
    for job in jobs:
        immport = job.get("immport")
        if immport is None:
            continue
        immport_source_text = _require_text(
            immport, "source_tab_zip", f"datasets[{job['name']}].immport"
        )
        immport_source = Path(immport_source_text)
        immport_records.extend(_read_immport_outcomes(immport, immport_source))
        immport_sources.append(
            {
                "dataset_name": job["name"],
                "path_text": immport_source_text,
                "path": immport_source,
                "selection": {
                    key: immport[key]
                    for key in [
                        "study_accession",
                        "outcome_names",
                        "timepoint_source_unit",
                        "field_precedence",
                    ]
                },
            }
        )
    if not immport_sources:
        raise ValidatorHandoffParseError("No ImmPort jobs were configured")
    immport_records = _deduplicate_records(immport_records, IMMPORT_COLUMNS)
    _write_tsv(immport_records, IMMPORT_COLUMNS, immport_output)
    immport_provenance = _combined_provenance(
        "immport_tab_zip",
        immport_sources,
        str(immport_output),
        immport_output,
        len(immport_records),
        IMMPORT_COLUMNS,
        immport_started_at,
        time.monotonic() - started,
    )
    _write_json(immport_provenance, immport_provenance_path)

    provenance_items = {
        "sample_manifest": (manifest_provenance, manifest_provenance_path),
        "feature_expression": (geo_provenance, geo_provenance_path),
        "quantitative_outcome": (immport_provenance, immport_provenance_path),
    }
    artifacts = {}
    for name, (provenance, provenance_path) in provenance_items.items():
        artifacts[name] = {
            **provenance["output"],
            "provenance": {
                "path": str(provenance_path),
                "size_bytes": provenance_path.stat().st_size,
                "sha256": _sha256(provenance_path),
            },
        }

    bundle_manifest = {
        "status": "success",
        "run_started_at_utc": run_started_at,
        "run_completed_at_utc": _utc_now(),
        "parser": "validator_handoff_parse_module",
        "parser_version": PARSER_VERSION,
        "python_version": platform.python_version(),
        "configuration": {
            "path": str(config_source),
            "size_bytes": config_source.stat().st_size,
            "sha256": _sha256(config_source),
        },
        "configuration_mode": "multi_dataset" if "datasets" in config else "legacy_single_dataset",
        "dataset_jobs": [job["name"] for job in jobs],
        "output_dir": output_dir_text,
        "artifacts": artifacts,
    }
    bundle_manifest_path = output_dir / BUNDLE_MANIFEST_FILENAME
    _write_json(bundle_manifest, bundle_manifest_path)
    return {
        "manifest": manifest_provenance,
        "geo": geo_provenance,
        "immport": immport_provenance,
        "bundle": bundle_manifest,
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a configured scientific-validator input bundle."
    )
    parser.add_argument("config", help="Path to validator handoff JSON configuration")
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        result = run_parser(args.config)
    except (ValidatorHandoffParseError, OSError, ValueError) as exc:
        logger.error("Validator handoff parsing failed: %s", exc)
        return 1
    logger.info(
        "Wrote %s manifest, %s GEO, and %s ImmPort rows",
        result["manifest"]["output"]["rows"],
        result["geo"]["output"]["rows"],
        result["immport"]["output"]["rows"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
