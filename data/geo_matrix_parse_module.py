"""Parse planned GEO series matrices into bounded analysis units.

This module reads plan rows marked ``download``, finds the cached series matrix that
contains each requested GSM/GPL set, and writes a filtered probe-by-sample matrix.
It preserves submitted expression values and performs no normalization, annotation,
eligibility filtering, or interpretation of sample characteristics.

The parser requires pandas but no network access. Install the project's dependencies
from the repository root with ``uv sync``.

------------------------------------------------------------------------------
How to run this program
------------------------------------------------------------------------------

Command-line usage:
    uv run python data/geo_matrix_parse_module.py \
        data/geo_cache/plan/geo_download_plan.tsv \
        --immport-manifest data/immport_cache/parsed/sample_manifest.tsv \
        --geo-cache-root data/geo_cache \
        --output-dir data/geo_cache/parsed

Python usage:
    from data.geo_matrix_parse_module import parse_geo_matrices

    provenance = parse_geo_matrices(
        "data/geo_cache/plan/geo_download_plan.tsv",
        "data/immport_cache/parsed/sample_manifest.tsv",
        "data/geo_cache",
        "data/geo_cache/parsed",
    )

Input:
    ``geo_download_plan.tsv`` is written by ``geo_plan_module.py``. The ImmPort
    manifest supplies subject, biosample, cohort, and timepoint fields for each GSM.
    The GEO cache must contain downloaded ``*_series_matrix.txt.gz`` files.

Output:
    Each study/experiment/GSE/GPL unit receives ``expression.tsv.gz``, ``samples.tsv``,
    ``geo_sample_metadata_long.tsv``, and ``provenance.json``. The output root also
    receives ``geo_matrix_parse_manifest.json`` summarizing all units and failures.

Failure behavior:
    Missing or ambiguous matrices, missing GSMs, platform mismatches, malformed tables,
    and duplicate identifiers are recorded as failures. Other units continue, but the
    command exits 1 if any unit fails.
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
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger("geo_matrix_parser")

PARSER_VERSION = "0.1.0"
EXPRESSION_FILENAME = "expression.tsv.gz"
SAMPLES_FILENAME = "samples.tsv"
GEO_METADATA_FILENAME = "geo_sample_metadata_long.tsv"
UNIT_PROVENANCE_FILENAME = "provenance.json"
RUN_MANIFEST_FILENAME = "geo_matrix_parse_manifest.json"

GSM_PATTERN = re.compile(r"^GSM\d+$", re.IGNORECASE)
GSE_PATTERN = re.compile(r"^GSE\d+$", re.IGNORECASE)
GPL_PATTERN = re.compile(r"^GPL\d+$", re.IGNORECASE)
SDY_PATTERN = re.compile(r"^SDY\d+$", re.IGNORECASE)
EXP_PATTERN = re.compile(r"^EXP\d+$", re.IGNORECASE)

REQUIRED_PLAN_COLUMNS = [
    "study_accession",
    "experiment_accession",
    "gse_accession",
    "gpl_accession",
    "download_recommendation",
    "requested_sample_count",
    "gsm_accessions",
]
REQUIRED_MANIFEST_COLUMNS = [
    "study_accession",
    "experiment_accession",
    "repository_name",
    "repository_accession",
]
UNIT_KEY_COLUMNS = [
    "study_accession",
    "experiment_accession",
    "gse_accession",
    "gpl_accession",
]


class GeoMatrixParseError(RuntimeError):
    """Raised when a matrix or linkage table violates the parsing contract."""


@dataclass
class MatrixHeader:
    """GEO sample metadata required to select and validate a series matrix."""

    gse_accession: str
    gsm_accessions: list[str]
    gpl_accessions: set[str]
    sample_metadata: pd.DataFrame


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_tsv(path: Path, required: list[str], label: str) -> pd.DataFrame:
    if not path.is_file():
        raise GeoMatrixParseError(f"{label} does not exist: {path}")
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise GeoMatrixParseError(f"Could not read {label} {path}: {exc}") from exc
    normalized = [str(column).strip().lower() for column in frame.columns]
    if len(normalized) != len(set(normalized)):
        raise GeoMatrixParseError(f"{label} contains duplicate normalized columns: {path}")
    frame.columns = normalized
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise GeoMatrixParseError(f"{label} is missing columns: {', '.join(missing)}")
    for column in frame.columns:
        frame[column] = frame[column].str.strip()
    return frame


def _line_fields(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\r\n")], delimiter="\t", quotechar='"'))


def _read_matrix_header(path: Path) -> MatrixHeader:
    sample_rows: list[tuple[str, list[str]]] = []
    gse_accession = ""
    table_found = False
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            for line in handle:
                if line.startswith("!series_matrix_table_begin"):
                    table_found = True
                    break
                fields = _line_fields(line)
                if not fields:
                    continue
                if fields[0] == "!Series_geo_accession" and len(fields) > 1:
                    gse_accession = fields[1].strip().upper()
                elif fields[0].startswith("!Sample_"):
                    sample_rows.append((fields[0][len("!Sample_") :], fields[1:]))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise GeoMatrixParseError(f"Could not read GEO matrix header {path}: {exc}") from exc

    if not table_found:
        raise GeoMatrixParseError(f"GEO matrix has no table-begin marker: {path}")
    accession_rows = [values for name, values in sample_rows if name == "geo_accession"]
    if len(accession_rows) != 1:
        raise GeoMatrixParseError(f"GEO matrix must have one Sample_geo_accession row: {path}")
    gsm_accessions = [value.strip().upper() for value in accession_rows[0]]
    if not gsm_accessions or any(not GSM_PATTERN.fullmatch(value) for value in gsm_accessions):
        raise GeoMatrixParseError(f"GEO matrix contains an invalid GSM accession: {path}")
    if len(gsm_accessions) != len(set(gsm_accessions)):
        raise GeoMatrixParseError(f"GEO matrix contains duplicate GSM accessions: {path}")

    metadata_records = []
    occurrences: Counter[str] = Counter()
    for attribute_order, (attribute, values) in enumerate(sample_rows):
        if len(values) != len(gsm_accessions):
            raise GeoMatrixParseError(
                f"Sample metadata width differs from GSM count for {attribute!r}: {path}"
            )
        occurrence = occurrences[attribute]
        occurrences[attribute] += 1
        for gsm_accession, value in zip(gsm_accessions, values, strict=True):
            metadata_records.append(
                {
                    "gsm_accession": gsm_accession,
                    "attribute": attribute,
                    "occurrence": occurrence,
                    "attribute_order": attribute_order,
                    "value": value,
                }
            )
    metadata = pd.DataFrame(metadata_records)
    platform_values = metadata.loc[
        metadata["attribute"].eq("platform_id"), "value"
    ].str.upper()
    gpl_accessions = {value for value in platform_values if value}
    if not gpl_accessions or any(not GPL_PATTERN.fullmatch(value) for value in gpl_accessions):
        raise GeoMatrixParseError(f"GEO matrix has invalid or missing platform metadata: {path}")
    if not GSE_PATTERN.fullmatch(gse_accession):
        raise GeoMatrixParseError(f"GEO matrix has invalid or missing series accession: {path}")
    return MatrixHeader(gse_accession, gsm_accessions, gpl_accessions, metadata)


def _read_expression_table(path: Path, expected_samples: list[str]) -> pd.DataFrame:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            for line in handle:
                if line.startswith("!series_matrix_table_begin"):
                    break
            else:
                raise GeoMatrixParseError(f"GEO matrix has no table-begin marker: {path}")
            frame = pd.read_csv(
                handle,
                sep="\t",
                dtype=str,
                keep_default_na=False,
                quotechar='"',
            )
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise GeoMatrixParseError(f"Could not read GEO expression table {path}: {exc}") from exc

    frame.columns = [str(column).strip() for column in frame.columns]
    if not len(frame.columns) or frame.columns[0] != "ID_REF":
        raise GeoMatrixParseError(f"GEO expression table must begin with ID_REF: {path}")
    end_marker = frame["ID_REF"].eq("!series_matrix_table_end")
    if int(end_marker.sum()) != 1 or not end_marker.iloc[-1]:
        raise GeoMatrixParseError(f"GEO expression table has an invalid end marker: {path}")
    frame = frame.loc[~end_marker].copy()
    if frame["ID_REF"].eq("").any() or frame["ID_REF"].duplicated().any():
        raise GeoMatrixParseError(f"GEO expression table has empty or duplicate probe IDs: {path}")
    matrix_samples = list(frame.columns[1:])
    if matrix_samples != expected_samples:
        raise GeoMatrixParseError(f"Expression columns do not match sample metadata order: {path}")
    return frame


def _requested_gsms(value: str, expected_count: str) -> list[str]:
    accessions = [item.strip().upper() for item in value.split(";") if item.strip()]
    if not accessions or any(not GSM_PATTERN.fullmatch(item) for item in accessions):
        raise GeoMatrixParseError("Plan row contains an invalid or empty GSM list")
    if len(accessions) != len(set(accessions)):
        raise GeoMatrixParseError("Plan row contains duplicate GSM accessions")
    try:
        count = int(expected_count)
    except ValueError as exc:
        raise GeoMatrixParseError(f"Invalid requested_sample_count: {expected_count!r}") from exc
    if count != len(accessions):
        raise GeoMatrixParseError(
            f"Plan count {count} does not match its {len(accessions)} GSM accessions"
        )
    return accessions


def _find_matrix(
    geo_cache_root: Path,
    gse_accession: str,
    gpl_accession: str,
    requested: list[str],
) -> tuple[Path, MatrixHeader]:
    candidates = sorted((geo_cache_root / gse_accession).glob("*_series_matrix.txt.gz"))
    matches = []
    requested_set = set(requested)
    for candidate in candidates:
        header = _read_matrix_header(candidate)
        if (
            header.gse_accession == gse_accession
            and gpl_accession in header.gpl_accessions
            and requested_set.issubset(header.gsm_accessions)
        ):
            matches.append((candidate, header))
    if not matches:
        raise GeoMatrixParseError(
            f"No matrix contains all requested {gse_accession}/{gpl_accession} GSMs"
        )
    if len(matches) > 1:
        paths = ", ".join(str(path) for path, _ in matches)
        raise GeoMatrixParseError(
            f"Multiple matrices contain requested {gse_accession}/{gpl_accession} GSMs: {paths}"
        )
    return matches[0]


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def _write_gzip_tsv(frame: pd.DataFrame, path: Path) -> None:
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(fileobj=raw_handle, mode="wb", filename="", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_handle:
                frame.to_csv(text_handle, sep="\t", index=False, lineterminator="\n")


def _unit_directory(output_root: Path, row: dict[str, str]) -> Path:
    return (
        output_root
        / row["study_accession"]
        / row["experiment_accession"]
        / f"{row['gse_accession']}_{row['gpl_accession']}"
    )


def _validate_unit_accessions(row: dict[str, str]) -> None:
    contracts = {
        "study_accession": SDY_PATTERN,
        "experiment_accession": EXP_PATTERN,
        "gse_accession": GSE_PATTERN,
        "gpl_accession": GPL_PATTERN,
    }
    for field, pattern in contracts.items():
        if not pattern.fullmatch(row[field]):
            raise GeoMatrixParseError(f"Plan row has invalid {field}: {row[field]!r}")


def _linked_samples(
    manifest: pd.DataFrame,
    row: dict[str, str],
    requested: list[str],
    source_matrix: Path,
) -> pd.DataFrame:
    repository_accessions = manifest["repository_accession"].str.upper()
    selected = manifest.loc[
        manifest["study_accession"].str.upper().eq(row["study_accession"])
        & manifest["experiment_accession"].str.upper().eq(row["experiment_accession"])
        & manifest["repository_name"].str.casefold().eq("geo")
        & repository_accessions.isin(requested)
    ].copy()
    selected["repository_accession"] = selected["repository_accession"].str.upper()
    counts = selected["repository_accession"].value_counts()
    missing = sorted(set(requested) - set(counts.index))
    duplicated = sorted(counts[counts.ne(1)].index)
    if missing or duplicated:
        raise GeoMatrixParseError(
            f"ImmPort linkage mismatch; missing={missing[:5]}, duplicated={duplicated[:5]}"
        )
    selected["gse_accession"] = row["gse_accession"]
    selected["gpl_accession"] = row["gpl_accession"]
    selected["source_matrix"] = str(source_matrix)
    order = {accession: index for index, accession in enumerate(requested)}
    selected["_sample_order"] = selected["repository_accession"].map(order)
    return selected.sort_values("_sample_order").drop(columns="_sample_order")


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_unit(
    row: dict[str, str],
    manifest: pd.DataFrame,
    geo_cache_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    started = time.monotonic()
    provenance_path: Path | None = None
    base_record: dict[str, Any] = {
        "study_accession": row["study_accession"],
        "experiment_accession": row["experiment_accession"],
        "gse_accession": row["gse_accession"],
        "gpl_accession": row["gpl_accession"],
        "requested_samples": 0,
        "run_started_at_utc": _utc_now(),
    }
    try:
        _validate_unit_accessions(row)
        requested = _requested_gsms(row["gsm_accessions"], row["requested_sample_count"])
        unit_dir = _unit_directory(output_root, row)
        unit_dir.mkdir(parents=True, exist_ok=True)
        provenance_path = unit_dir / UNIT_PROVENANCE_FILENAME
        base_record.update(
            {
                "requested_samples": len(requested),
                "output_directory": str(unit_dir),
            }
        )
        matrix_path, header = _find_matrix(
            geo_cache_root,
            row["gse_accession"],
            row["gpl_accession"],
            requested,
        )
        full_expression = _read_expression_table(matrix_path, header.gsm_accessions)
        expression = full_expression.loc[:, ["ID_REF", *requested]]
        samples = _linked_samples(manifest, row, requested, matrix_path)
        metadata = header.sample_metadata.loc[
            header.sample_metadata["gsm_accession"].isin(requested)
        ].copy()
        metadata["_sample_order"] = metadata["gsm_accession"].map(
            {accession: index for index, accession in enumerate(requested)}
        )
        metadata = metadata.sort_values(
            ["_sample_order", "attribute_order", "occurrence"]
        ).drop(columns="_sample_order")

        expression_path = unit_dir / EXPRESSION_FILENAME
        samples_path = unit_dir / SAMPLES_FILENAME
        metadata_path = unit_dir / GEO_METADATA_FILENAME
        _write_gzip_tsv(expression, expression_path)
        _write_tsv(samples, samples_path)
        _write_tsv(metadata, metadata_path)
        record = {
            **base_record,
            "status": "success",
            "source_matrix": {
                "path": str(matrix_path),
                "size_bytes": matrix_path.stat().st_size,
                "sha256": _sha256(matrix_path),
            },
            "counts": {
                "probes": len(expression),
                "samples": len(requested),
                "geo_metadata_rows": len(metadata),
            },
            "outputs": [
                {
                    "type": output_type,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
                for output_type, path in [
                    ("expression", expression_path),
                    ("samples", samples_path),
                    ("geo_sample_metadata", metadata_path),
                ]
            ],
            "duration_sec": round(time.monotonic() - started, 3),
        }
    except Exception as exc:
        record = {
            **base_record,
            "status": "failed",
            "error": str(exc),
            "duration_sec": round(time.monotonic() - started, 3),
        }
        if provenance_path is not None:
            _write_json(record, provenance_path)
        return record

    assert provenance_path is not None
    _write_json(record, provenance_path)
    return record


def parse_geo_matrices(
    plan_path: str | Path,
    immport_manifest_path: str | Path,
    geo_cache_root: str | Path,
    output_dir: str | Path,
) -> dict[str, Any]:
    """Parse all plan rows marked download and return run-level provenance."""
    started = time.monotonic()
    started_at = _utc_now()
    plan_source = Path(plan_path).expanduser().resolve()
    manifest_source = Path(immport_manifest_path).expanduser().resolve()
    cache_root = Path(geo_cache_root).expanduser().resolve()
    output_root = Path(output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    plan = _read_tsv(plan_source, REQUIRED_PLAN_COLUMNS, "GEO download plan")
    manifest = _read_tsv(
        manifest_source,
        REQUIRED_MANIFEST_COLUMNS,
        "ImmPort sample manifest",
    )
    selected = plan.loc[plan["download_recommendation"].eq("download")].copy()
    if selected.empty:
        raise GeoMatrixParseError("GEO download plan has no rows marked download")
    if selected.duplicated(UNIT_KEY_COLUMNS).any():
        raise GeoMatrixParseError("GEO download plan contains duplicate analysis units")

    units = []
    for row in selected.to_dict(orient="records"):
        record = _parse_unit(row, manifest, cache_root, output_root)
        if record["status"] == "failed":
            logger.error(
                "Failed to parse %s/%s/%s/%s: %s",
                row["study_accession"],
                row["experiment_accession"],
                row["gse_accession"],
                row["gpl_accession"],
                record["error"],
            )
        units.append(record)

    failed = sum(unit["status"] == "failed" for unit in units)
    succeeded = len(units) - failed
    run_status = "success" if not failed else ("failed" if not succeeded else "partial")
    provenance = {
        "status": run_status,
        "run_started_at_utc": started_at,
        "parser": "geo_matrix_parse_module",
        "parser_version": PARSER_VERSION,
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "inputs": [
            {
                "type": input_type,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for input_type, path in [
                ("geo_download_plan", plan_source),
                ("immport_sample_manifest", manifest_source),
            ]
        ],
        "geo_cache_root": str(cache_root),
        "output_directory": str(output_root),
        "counts": {
            "plan_rows": len(plan),
            "download_units": len(units),
            "skipped_plan_rows": len(plan) - len(selected),
            "units_succeeded": succeeded,
            "units_failed": failed,
            "requested_samples": sum(unit["requested_samples"] for unit in units),
            "parsed_samples": sum(
                unit.get("counts", {}).get("samples", 0) for unit in units
            ),
        },
        "units": units,
        "duration_sec": round(time.monotonic() - started, 3),
    }
    _write_json(provenance, output_root / RUN_MANIFEST_FILENAME)
    return provenance


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse planned cached GEO matrices without normalization."
    )
    parser.add_argument("plan", help="Path to geo_download_plan.tsv")
    parser.add_argument("--immport-manifest", required=True)
    parser.add_argument("--geo-cache-root", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        provenance = parse_geo_matrices(
            args.plan,
            args.immport_manifest,
            args.geo_cache_root,
            args.output_dir,
        )
    except (GeoMatrixParseError, OSError, ValueError) as exc:
        logger.error("GEO matrix parsing failed: %s", exc)
        return 1

    counts = provenance["counts"]
    logger.info(
        "Parsed %s/%s units and %s samples",
        counts["units_succeeded"],
        counts["download_units"],
        counts["parsed_samples"],
    )
    return 0 if provenance["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
