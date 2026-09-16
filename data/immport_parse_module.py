"""Build a normalized ImmPort experimental-sample manifest.

This module joins ImmPort study, experiment, biosample, subject, cohort, and public
repository tables. It produces one row per experimental sample without selecting
assays, subjects, timepoints, or repository-linked samples.

The parser requires pandas but no network access or ImmPort credentials. Install the
project's pinned dependencies from the repository root with ``uv sync``.

------------------------------------------------------------------------------
How to run this program
------------------------------------------------------------------------------

Command-line usage:
    Parse an extracted ImmPort Tab directory:

    uv run python data/immport_parse_module.py \
        data/immport_cache/SDY1529/SDY1529-DR58_Tab/Tab \
        --output-dir data/immport_cache/SDY1529/parsed

    The input may instead be the unextracted Tab ZIP:

    uv run python data/immport_parse_module.py \
        data/immport_cache/SDY1529/SDY1529-DR58_Tab.zip \
        --output-dir data/immport_cache/SDY1529/parsed

    Run the script with ``--help`` for argument details.

Python usage:
    Use ``parse_immport_sample_links`` to return a DataFrame without writing files:

    from data.immport_parse_module import parse_immport_sample_links

    samples = parse_immport_sample_links(
        "data/immport_cache/SDY1529/SDY1529-DR58_Tab.zip"
    )

    Use ``run_parser`` to write the manifest and provenance and return the provenance
    dictionary:

    from data.immport_parse_module import run_parser

    provenance = run_parser(
        "data/immport_cache/SDY1529/SDY1529-DR58_Tab.zip",
        "data/immport_cache/SDY1529/parsed",
    )

Input:
    The source must contain ``experiment.txt``, ``expsample.txt``,
    ``expsample_2_biosample.txt``, ``biosample.txt``, ``subject.txt``,
    ``arm_2_subject.txt``, ``arm_or_cohort.txt``, and
    ``expsample_public_repository.txt``. ZIP members are read directly without
    extracting the archive.

Output:
    ``sample_manifest.tsv`` contains normalized linkage and study-design fields.
    ``sample_manifest.provenance.json`` records source members, SHA-256 hashes,
    versions, row counts, output paths, status, errors, and timing.

Validation and failures:
    Required tables and columns, unique join keys, join cardinality, and study
    accessions are validated. Missing public-repository links remain blank. Invalid
    input raises ``ImmportParseError``; command-line failures log an error and exit 1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

import pandas as pd

logger = logging.getLogger("immport_parser")

PARSER_VERSION = "0.1.0"
OUTPUT_FILENAME = "sample_manifest.tsv"
PROVENANCE_FILENAME = "sample_manifest.provenance.json"

TABLE_COLUMNS = {
    "experiment.txt": [
        "experiment_accession",
        "measurement_technique",
        "name",
        "study_accession",
    ],
    "expsample.txt": [
        "expsample_accession",
        "experiment_accession",
        "result_schema",
        "upload_result_status",
    ],
    "expsample_2_biosample.txt": [
        "expsample_accession",
        "biosample_accession",
    ],
    "biosample.txt": [
        "biosample_accession",
        "planned_visit_accession",
        "study_accession",
        "study_time_collected",
        "study_time_collected_unit",
        "study_time_t0_event",
        "subject_accession",
        "type",
        "subtype",
    ],
    "subject.txt": [
        "subject_accession",
        "ancestral_population",
        "ethnicity",
        "gender",
        "race",
        "race_specify",
        "species",
        "strain",
    ],
    "arm_2_subject.txt": [
        "arm_accession",
        "subject_accession",
        "age_event",
        "age_unit",
        "max_subject_age",
        "min_subject_age",
        "subject_phenotype",
        "subject_location",
        "max_subject_age_in_years",
        "min_subject_age_in_years",
    ],
    "arm_or_cohort.txt": [
        "arm_accession",
        "name",
        "study_accession",
        "type_reported",
        "type_preferred",
    ],
    "expsample_public_repository.txt": [
        "expsample_accession",
        "repository_accession",
        "repository_name",
    ],
}

OUTPUT_COLUMNS = [
    "study_accession",
    "experiment_accession",
    "experiment_name",
    "measurement_technique",
    "expsample_accession",
    "result_schema",
    "upload_result_status",
    "repository_name",
    "repository_accession",
    "biosample_accession",
    "planned_visit_accession",
    "subject_accession",
    "study_time_collected",
    "study_time_collected_unit",
    "study_time_t0_event",
    "biosample_type",
    "biosample_subtype",
    "arm_accession",
    "arm_name",
    "arm_type_reported",
    "arm_type_preferred",
    "age_event",
    "age_unit",
    "min_subject_age",
    "max_subject_age",
    "min_subject_age_in_years",
    "max_subject_age_in_years",
    "subject_phenotype",
    "subject_location",
    "ancestral_population",
    "ethnicity",
    "gender",
    "race",
    "race_specify",
    "species",
    "strain",
]


class ImmportParseError(RuntimeError):
    """Raised when required ImmPort tables do not satisfy the join contract."""


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalize_table(frame: pd.DataFrame, source: str, filename: str) -> pd.DataFrame:
    """Normalize and validate one table loaded from a directory or ZIP archive."""
    normalized = [str(column).strip().lower() for column in frame.columns]
    if len(normalized) != len(set(normalized)):
        raise ImmportParseError(f"Duplicate normalized columns in {source}")
    frame.columns = normalized

    required = TABLE_COLUMNS[filename]
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ImmportParseError(f"{source} is missing required columns: {', '.join(missing)}")

    selected = frame.loc[:, required].copy()
    for column in selected.columns:
        selected[column] = selected[column].str.strip()
    return selected


def _read_stream(handle: Any, source: str, filename: str) -> pd.DataFrame:
    try:
        frame = pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise ImmportParseError(f"Could not read {source}: {exc}") from exc
    return _normalize_table(frame, source, filename)


def _zip_members(archive: zipfile.ZipFile) -> dict[str, zipfile.ZipInfo]:
    """Resolve each required table to exactly one member below a Tab directory."""
    members = {}
    for filename in TABLE_COLUMNS:
        matches = [
            info
            for info in archive.infolist()
            if not info.is_dir()
            and PurePosixPath(info.filename).name == filename
            and "Tab" in PurePosixPath(info.filename).parts[:-1]
        ]
        if not matches:
            raise ImmportParseError(
                f"Required ImmPort table is missing from {archive.filename}: {filename}"
            )
        if len(matches) > 1:
            paths = ", ".join(sorted(info.filename for info in matches))
            raise ImmportParseError(
                f"Multiple {filename} members found in {archive.filename}: {paths}"
            )
        members[filename] = matches[0]
    return members


def _load_tables(tab_source: Path) -> dict[str, pd.DataFrame]:
    if tab_source.is_dir():
        tables = {}
        for filename in TABLE_COLUMNS:
            path = tab_source / filename
            if not path.is_file():
                raise ImmportParseError(f"Required ImmPort table is missing: {path}")
            with path.open("rb") as handle:
                tables[filename] = _read_stream(handle, str(path), filename)
        return tables

    if not tab_source.is_file():
        raise ImmportParseError(f"ImmPort Tab source does not exist: {tab_source}")
    if tab_source.suffix.lower() != ".zip":
        raise ImmportParseError(f"ImmPort Tab source must be a directory or ZIP: {tab_source}")

    try:
        with zipfile.ZipFile(tab_source) as archive:
            members = _zip_members(archive)
            tables = {}
            for filename, member in members.items():
                source = f"{tab_source}!{member.filename}"
                with archive.open(member) as handle:
                    tables[filename] = _read_stream(handle, source, filename)
            return tables
    except zipfile.BadZipFile as exc:
        raise ImmportParseError(f"Invalid ZIP archive {tab_source}: {exc}") from exc


def _validate_key(frame: pd.DataFrame, key: str, filename: str) -> None:
    if frame[key].eq("").any():
        raise ImmportParseError(f"{filename} contains an empty {key}")
    duplicates = frame.loc[frame[key].duplicated(keep=False), key].unique()
    if len(duplicates):
        examples = ", ".join(sorted(duplicates)[:5])
        raise ImmportParseError(f"{filename} contains duplicate {key} values: {examples}")


def _validate_study_accessions(frame: pd.DataFrame) -> None:
    columns = [
        "study_accession",
        "biosample_study_accession",
        "arm_study_accession",
    ]
    mismatched = pd.Series(False, index=frame.index)
    for column in columns[1:]:
        mismatched |= frame[column].ne("") & frame["study_accession"].ne(frame[column])
    if mismatched.any():
        samples = frame.loc[mismatched, "expsample_accession"].head(5).tolist()
        raise ImmportParseError(
            "Study accessions disagree across joined tables for experimental samples: "
            + ", ".join(samples)
        )


def parse_immport_sample_links(tab_source: str | Path) -> pd.DataFrame:
    """Return one normalized row per ImmPort experimental sample.

    All joins are structural and left-preserving. Missing optional links remain blank;
    missing required tables, columns, keys, or one-to-one relationships raise
    ``ImmportParseError``.
    """
    source_path = Path(tab_source).expanduser().resolve()
    tables = _load_tables(source_path)

    key_contracts = {
        "experiment.txt": "experiment_accession",
        "expsample.txt": "expsample_accession",
        "expsample_2_biosample.txt": "expsample_accession",
        "biosample.txt": "biosample_accession",
        "subject.txt": "subject_accession",
        "arm_2_subject.txt": "subject_accession",
        "arm_or_cohort.txt": "arm_accession",
        "expsample_public_repository.txt": "expsample_accession",
    }
    for filename, key in key_contracts.items():
        _validate_key(tables[filename], key, filename)

    experiment = tables["experiment.txt"].rename(columns={"name": "experiment_name"})
    biosample = tables["biosample.txt"].rename(
        columns={
            "study_accession": "biosample_study_accession",
            "type": "biosample_type",
            "subtype": "biosample_subtype",
        }
    )
    arm = tables["arm_or_cohort.txt"].rename(
        columns={
            "name": "arm_name",
            "study_accession": "arm_study_accession",
            "type_reported": "arm_type_reported",
            "type_preferred": "arm_type_preferred",
        }
    )

    frame = tables["expsample.txt"].merge(
        experiment,
        on="experiment_accession",
        how="left",
        validate="many_to_one",
    )
    frame = frame.merge(
        tables["expsample_2_biosample.txt"],
        on="expsample_accession",
        how="left",
        validate="one_to_one",
    )
    frame = frame.merge(
        biosample,
        on="biosample_accession",
        how="left",
        validate="many_to_one",
    )
    frame = frame.merge(
        tables["subject.txt"],
        on="subject_accession",
        how="left",
        validate="many_to_one",
    )
    frame = frame.merge(
        tables["arm_2_subject.txt"],
        on="subject_accession",
        how="left",
        validate="many_to_one",
    )
    frame = frame.merge(
        arm,
        on="arm_accession",
        how="left",
        validate="many_to_one",
    )
    frame = frame.merge(
        tables["expsample_public_repository.txt"],
        on="expsample_accession",
        how="left",
        validate="one_to_one",
    ).fillna("")

    if len(frame) != len(tables["expsample.txt"]):
        raise ImmportParseError("A join changed the number of experimental samples")
    _validate_study_accessions(frame)

    return frame.loc[:, OUTPUT_COLUMNS].sort_values("expsample_accession").reset_index(drop=True)


def _input_provenance(tab_source: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if tab_source.is_file() and tab_source.suffix.lower() == ".zip":
        try:
            with zipfile.ZipFile(tab_source) as archive:
                members = _zip_members(archive)
                inputs = [
                    {
                        "filename": filename,
                        "member": member.filename,
                        "size_bytes": member.file_size,
                        "crc32": f"{member.CRC:08x}",
                    }
                    for filename, member in members.items()
                ]
        except zipfile.BadZipFile as exc:
            raise ImmportParseError(f"Invalid ZIP archive {tab_source}: {exc}") from exc
        source = {
            "type": "zip",
            "path": str(tab_source),
            "size_bytes": tab_source.stat().st_size,
            "sha256": _sha256(tab_source),
        }
        return source, inputs

    records = []
    for filename in TABLE_COLUMNS:
        path = tab_source / filename
        if not path.is_file():
            continue
        records.append(
            {
                "filename": filename,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    return {"type": "directory", "path": str(tab_source)}, records


def run_parser(tab_source: str | Path, output_dir: str | Path) -> dict[str, Any]:
    """Parse an ImmPort Tab directory or ZIP, write outputs, and return provenance."""
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    source_path = Path(tab_source).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    output_path = destination / OUTPUT_FILENAME
    provenance_path = destination / PROVENANCE_FILENAME

    source_provenance, inputs = _input_provenance(source_path)
    base_provenance: dict[str, Any] = {
        "run_started_at_utc": started_at,
        "parser": "immport_parse_module",
        "parser_version": PARSER_VERSION,
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "input_source": source_provenance,
        "inputs": inputs,
    }

    try:
        frame = parse_immport_sample_links(source_path)
        frame.to_csv(output_path, sep="\t", index=False, lineterminator="\n")
        provenance = {
            **base_provenance,
            "status": "success",
            "output": {
                "path": str(output_path),
                "size_bytes": output_path.stat().st_size,
                "sha256": _sha256(output_path),
                "rows": len(frame),
                "columns": list(frame.columns),
            },
            "counts": {
                "experimental_samples": len(frame),
                "repository_linked_samples": int(frame["repository_accession"].ne("").sum()),
            },
            "duration_sec": round(time.monotonic() - started, 3),
        }
    except Exception as exc:
        provenance = {
            **base_provenance,
            "status": "failed",
            "error": str(exc),
            "duration_sec": round(time.monotonic() - started, 3),
        }
        provenance_path.write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        raise

    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return provenance


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a normalized sample-linkage manifest from ImmPort Tab tables."
    )
    parser.add_argument("tab_source", help="Path to an ImmPort Tab directory or *_Tab.zip")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for sample_manifest.tsv and its provenance JSON",
    )
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        provenance = run_parser(args.tab_source, args.output_dir)
    except (ImmportParseError, OSError, ValueError) as exc:
        logger.error("ImmPort parsing failed: %s", exc)
        return 1

    logger.info(
        "Wrote %s experimental samples to %s",
        provenance["counts"]["experimental_samples"],
        provenance["output"]["path"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
