"""Create bounded GEO and ImmPort tables for scientific-validator handoff.

Inputs are a JSON configuration, one cached GEO series matrix, and one ImmPort
Tab ZIP. Outputs are validator-ready TSV files and provenance JSON records. The
configuration explicitly identifies the GEO feature and ImmPort outcome; this
module does not infer scientific eligibility or normalize source measurements.

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
import time
import zipfile
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path, PurePosixPath
from typing import Any

logger = logging.getLogger("validator_handoff_parser")

PARSER_VERSION = "0.1.0"
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
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


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


def run_parser(config_path: str | Path) -> dict[str, Any]:
    """Write both validator handoff tables and return their provenance records."""
    config_source = Path(config_path)
    try:
        config = json.loads(config_source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValidatorHandoffParseError(f"Could not read configuration: {exc}") from exc
    if not isinstance(config, dict):
        raise ValidatorHandoffParseError("Configuration root must be an object")
    geo = config.get("geo")
    immport = config.get("immport")
    if not isinstance(geo, dict) or not isinstance(immport, dict):
        raise ValidatorHandoffParseError("Configuration requires geo and immport objects")

    geo_source_text = _require_text(geo, "source_matrix", "geo")
    geo_output_text = _require_text(geo, "output", "geo")
    geo_provenance_text = _require_text(geo, "provenance", "geo")
    geo_source = Path(geo_source_text)
    geo_output = Path(geo_output_text)
    geo_started_at = _utc_now()
    started = time.monotonic()
    geo_records = _read_geo_feature(geo, geo_source)
    _write_tsv(geo_records, GEO_COLUMNS, geo_output)
    geo_provenance = _provenance(
        "geo_series_matrix",
        geo_source_text,
        geo_source,
        geo_output_text,
        geo_output,
        len(geo_records),
        GEO_COLUMNS,
        {
            key: geo[key]
            for key in [
                "study_accession",
                "gse_accession",
                "gpl_accession",
                "gene",
                "feature_id",
            ]
        },
        geo_started_at,
        time.monotonic() - started,
    )
    geo_provenance_path = Path(geo_provenance_text)
    geo_provenance_path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(geo_provenance, geo_provenance_path)

    immport_source_text = _require_text(immport, "source_tab_zip", "immport")
    immport_output_text = _require_text(immport, "output", "immport")
    immport_provenance_text = _require_text(immport, "provenance", "immport")
    immport_source = Path(immport_source_text)
    immport_output = Path(immport_output_text)
    immport_started_at = _utc_now()
    started = time.monotonic()
    immport_records = _read_immport_outcomes(immport, immport_source)
    _write_tsv(immport_records, IMMPORT_COLUMNS, immport_output)
    immport_provenance = _provenance(
        "immport_tab_zip",
        immport_source_text,
        immport_source,
        immport_output_text,
        immport_output,
        len(immport_records),
        IMMPORT_COLUMNS,
        {
            key: immport[key]
            for key in [
                "study_accession",
                "outcome_names",
                "timepoint_source_unit",
                "field_precedence",
            ]
        },
        immport_started_at,
        time.monotonic() - started,
    )
    immport_provenance_path = Path(immport_provenance_text)
    immport_provenance_path.parent.mkdir(parents=True, exist_ok=True)
    _write_json(immport_provenance, immport_provenance_path)
    return {"geo": geo_provenance, "immport": immport_provenance}


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create configured GEO feature and ImmPort outcome handoff tables."
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
        "Wrote %s GEO rows and %s ImmPort rows",
        result["geo"]["output"]["rows"],
        result["immport"]["output"]["rows"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
