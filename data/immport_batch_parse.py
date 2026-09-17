"""Parse every current ImmPort Tab source discovered under a cache root.

This module discovers ImmPort study directories, selects one Tab source per study,
calls ``immport_parse_module``, and combines successful study manifests. Study
accessions are discovered from the filesystem and are not hardcoded.

The batch parser requires pandas and ``immport_parse_module`` but no network access
or ImmPort credentials. Install dependencies from the repository root with ``uv sync``.

------------------------------------------------------------------------------
How to run this program
------------------------------------------------------------------------------

Command-line usage:
    Parse all studies under the default ``data/immport_cache`` root:

    uv run python data/immport_batch_parse.py

    Parse a cache located elsewhere:

    uv run python data/immport_batch_parse.py \
        --cache-root path/to/immport_cache

    Run the script with ``--help`` for argument details.

Python usage:
    Import and call ``run_batch``. It writes outputs and returns the batch provenance
    dictionary:

    from data.immport_batch_parse import run_batch

    provenance = run_batch("data/immport_cache")

Input layout and source selection:
    Each direct child study directory may contain an extracted
    ``<SDY_ID>-DR<n>_Tab/Tab/`` directory, an unextracted
    ``<SDY_ID>-DR<n>_Tab.zip``, or both. The highest numbered data release is used.
    For the same release, extracted Tab tables are preferred. MySQL ZIPs and archives
    nested under other directories are ignored.

Output:
    Each successful study writes ``<SDY_ID>/parsed/sample_manifest.tsv`` and its
    provenance JSON. The cache root also receives ``parsed/sample_manifest.tsv``
    containing all successful studies, ``parsed/geo_series_links.tsv`` containing
    structured GSE links extracted from ``study_link.txt``, and
    ``parsed/batch_manifest.json`` containing selected sources, study statuses,
    counts, errors, output hashes, and timing.

Failure behavior:
    A study failure is recorded and does not stop later studies. The command exits 1
    if any study fails or a GEO-linked study has no structured GSE link; otherwise it
    exits 0. If no usable Tab source is discovered, ``ImmportBatchParseError`` is raised.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import re
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

import pandas as pd

try:
    from data.immport_parse_module import (
        OUTPUT_COLUMNS,
        OUTPUT_FILENAME,
        PARSER_VERSION,
        run_parser,
    )
except ModuleNotFoundError:
    from immport_parse_module import (  # type: ignore[no-redef]
        OUTPUT_COLUMNS,
        OUTPUT_FILENAME,
        PARSER_VERSION,
        run_parser,
    )

logger = logging.getLogger("immport_batch_parser")

BATCH_VERSION = "0.2.0"
DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent / "immport_cache"
SOURCE_PATTERN = re.compile(
    r"^(?P<accession>SDY\d+)-DR(?P<release>\d+)_Tab(?P<zip>\.zip)?$",
    re.IGNORECASE,
)
COMBINED_FILENAME = "sample_manifest.tsv"
GEO_SERIES_LINKS_FILENAME = "geo_series_links.tsv"
BATCH_PROVENANCE_FILENAME = "batch_manifest.json"
STUDY_LINK_FILENAME = "study_link.txt"
GSE_LINK_PATTERN = re.compile(r"\bGSE\d+\b", re.IGNORECASE)
GEO_SERIES_LINK_COLUMNS = [
    "study_accession",
    "gse_accession",
    "link_name",
    "link_value",
    "data_release",
    "input_source",
]


class ImmportBatchParseError(RuntimeError):
    """Raised when the cache root contains no usable ImmPort Tab sources."""


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _study_link_member(archive: zipfile.ZipFile, source: Path) -> zipfile.ZipInfo | None:
    matches = [
        info
        for info in archive.infolist()
        if not info.is_dir()
        and PurePosixPath(info.filename).name == STUDY_LINK_FILENAME
        and "Tab" in PurePosixPath(info.filename).parts[:-1]
    ]
    if len(matches) > 1:
        paths = ", ".join(sorted(info.filename for info in matches))
        raise ImmportBatchParseError(
            f"Multiple {STUDY_LINK_FILENAME} members found in {source}: {paths}"
        )
    return matches[0] if matches else None


def _normalize_study_links(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    frame.columns = [str(column).strip().lower() for column in frame.columns]
    if len(frame.columns) != len(set(frame.columns)):
        raise ImmportBatchParseError(f"Duplicate normalized columns in {source}")
    required = {"name", "study_accession", "value"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ImmportBatchParseError(
            f"{source} is missing required columns: {', '.join(missing)}"
        )
    selected = frame.loc[:, ["name", "study_accession", "value"]].copy()
    for column in selected.columns:
        selected[column] = selected[column].astype(str).str.strip()
    return selected


def _read_geo_series_links(
    accession: str,
    release: int,
    source: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Extract structured GSE accessions from one optional ImmPort study-link table."""
    if source.is_dir():
        table_path = source / STUDY_LINK_FILENAME
        if not table_path.is_file():
            return [], {"status": "missing", "path": str(table_path)}
        try:
            frame = pd.read_csv(table_path, sep="\t", dtype=str, keep_default_na=False)
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            raise ImmportBatchParseError(f"Could not read {table_path}: {exc}") from exc
        input_record = {
            "status": "success",
            "type": "file",
            "path": str(table_path),
            "size_bytes": table_path.stat().st_size,
            "sha256": _sha256(table_path),
        }
        source_reference = str(table_path)
    else:
        try:
            with zipfile.ZipFile(source) as archive:
                member = _study_link_member(archive, source)
                if member is None:
                    return [], {
                        "status": "missing",
                        "type": "zip_member",
                        "path": str(source),
                        "member": None,
                        "sha256": _sha256(source),
                    }
                with archive.open(member) as handle:
                    frame = pd.read_csv(
                        handle,
                        sep="\t",
                        dtype=str,
                        keep_default_na=False,
                    )
                member_name = member.filename
                member_size = member.file_size
                member_crc32 = f"{member.CRC:08x}"
        except (OSError, UnicodeError, zipfile.BadZipFile, pd.errors.ParserError) as exc:
            raise ImmportBatchParseError(f"Could not read {source}: {exc}") from exc
        input_record = {
            "status": "success",
            "type": "zip_member",
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "sha256": _sha256(source),
            "member": member_name,
            "member_size_bytes": member_size,
            "member_crc32": member_crc32,
        }
        source_reference = f"{source}!{member_name}"

    normalized = _normalize_study_links(frame, source_reference)
    records: list[dict[str, Any]] = []
    for row_number, row in normalized.iterrows():
        row_accession = row["study_accession"].upper()
        if row_accession != accession:
            raise ImmportBatchParseError(
                f"{source_reference} row {row_number + 2} belongs to "
                f"{row_accession!r}, not {accession}"
            )
        for gse_accession in sorted(set(GSE_LINK_PATTERN.findall(row["value"].upper()))):
            records.append(
                {
                    "study_accession": accession,
                    "gse_accession": gse_accession,
                    "link_name": row["name"],
                    "link_value": row["value"],
                    "data_release": release,
                    "input_source": source_reference,
                }
            )
    return records, input_record


def discover_tab_sources(cache_root: str | Path) -> dict[str, tuple[int, Path]]:
    """Return the highest-release Tab source for each discovered study directory."""
    root = Path(cache_root).expanduser().resolve()
    if not root.is_dir():
        raise ImmportBatchParseError(f"ImmPort cache root does not exist: {root}")

    selected: dict[str, tuple[int, Path]] = {}
    for study_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        candidates: list[tuple[Path, Path]] = []
        for candidate in study_dir.iterdir():
            if candidate.is_file() and candidate.suffix.lower() == ".zip":
                candidates.append((candidate, candidate))
            elif candidate.is_dir() and (candidate / "Tab").is_dir():
                candidates.append((candidate, candidate / "Tab"))

        for named_source, parser_source in sorted(candidates):
            match = SOURCE_PATTERN.fullmatch(named_source.name)
            if not match:
                continue
            accession = match.group("accession").upper()
            if study_dir.name.upper() != accession:
                continue
            release = int(match.group("release"))
            current = selected.get(accession)
            prefer_extracted = parser_source.is_dir()
            current_is_zip = current is not None and current[1].is_file()
            if (
                current is None
                or release > current[0]
                or (release == current[0] and prefer_extracted and current_is_zip)
            ):
                selected[accession] = (release, parser_source.resolve())

    if not selected:
        raise ImmportBatchParseError(f"No extracted Tab directories or Tab ZIPs found under {root}")
    return selected


def run_batch(cache_root: str | Path) -> dict[str, Any]:
    """Parse all discovered studies and write combined output and provenance."""
    started = time.monotonic()
    started_at = datetime.now(timezone.utc).isoformat()
    root = Path(cache_root).expanduser().resolve()
    selected = discover_tab_sources(root)
    combined_dir = root / "parsed"
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_path = combined_dir / COMBINED_FILENAME
    geo_series_links_path = combined_dir / GEO_SERIES_LINKS_FILENAME
    provenance_path = combined_dir / BATCH_PROVENANCE_FILENAME

    study_records = []
    frames = []
    geo_series_link_records: list[dict[str, Any]] = []
    for accession, (release, source) in sorted(selected.items()):
        output_dir = root / accession / "parsed"
        record: dict[str, Any] = {
            "study_accession": accession,
            "data_release": release,
            "input_source": str(source),
            "input_source_type": "directory" if source.is_dir() else "zip",
            "output_directory": str(output_dir),
        }
        try:
            study_provenance = run_parser(source, output_dir)
            frame = pd.read_csv(
                output_dir / OUTPUT_FILENAME,
                sep="\t",
                dtype=str,
                keep_default_na=False,
            )
            study_geo_links, study_link_input = _read_geo_series_links(
                accession,
                release,
                source,
            )
            frames.append(frame)
            geo_series_link_records.extend(study_geo_links)
            record.update(
                {
                    "status": "success",
                    "experimental_samples": len(frame),
                    "repository_linked_samples": study_provenance["counts"][
                        "repository_linked_samples"
                    ],
                    "geo_series_links": len(study_geo_links),
                    "study_link_input": study_link_input,
                }
            )
        except Exception as exc:
            logger.error("Failed to parse %s: %s", accession, exc)
            record.update({"status": "failed", "error": str(exc)})
        study_records.append(record)

    if frames:
        combined = pd.concat(frames, ignore_index=True)
        combined = combined.loc[:, OUTPUT_COLUMNS].sort_values(
            ["study_accession", "expsample_accession"]
        )
    else:
        combined = pd.DataFrame(columns=OUTPUT_COLUMNS)
    combined.to_csv(combined_path, sep="\t", index=False, lineterminator="\n")

    geo_series_links = pd.DataFrame(
        geo_series_link_records,
        columns=GEO_SERIES_LINK_COLUMNS,
    )
    if not geo_series_links.empty:
        geo_series_links = geo_series_links.drop_duplicates().sort_values(
            ["study_accession", "gse_accession", "link_value"]
        )
    geo_series_links.to_csv(
        geo_series_links_path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    geo_sample_mask = (
        combined["repository_name"].str.casefold().eq("geo")
        & combined["repository_accession"].str.fullmatch(r"(?i)GSM\d+")
    )
    geo_linked_studies = set(combined.loc[geo_sample_mask, "study_accession"])
    series_linked_studies = set(geo_series_links["study_accession"])
    missing_geo_series_links = sorted(geo_linked_studies - series_linked_studies)
    records_by_study = {record["study_accession"]: record for record in study_records}
    for accession in geo_linked_studies:
        record = records_by_study.get(accession)
        if record is not None:
            record["geo_series_link_status"] = (
                "missing" if accession in missing_geo_series_links else "resolved"
            )

    failures = sum(record["status"] == "failed" for record in study_records)
    status = "partial" if failures or missing_geo_series_links else "success"
    provenance = {
        "status": status,
        "run_started_at_utc": started_at,
        "batch_parser": "immport_batch_parse",
        "batch_parser_version": BATCH_VERSION,
        "sample_parser_version": PARSER_VERSION,
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "cache_root": str(root),
        "studies": study_records,
        "counts": {
            "studies_discovered": len(study_records),
            "studies_succeeded": len(study_records) - failures,
            "studies_failed": failures,
            "experimental_samples": len(combined),
            "repository_linked_samples": int(
                combined["repository_accession"].ne("").sum()
            ),
            "geo_linked_studies": len(geo_linked_studies),
            "geo_series_links": len(geo_series_links),
            "geo_series_accessions": geo_series_links["gse_accession"].nunique(),
            "geo_linked_studies_without_series_link": len(missing_geo_series_links),
        },
        "combined_output": {
            "path": str(combined_path),
            "size_bytes": combined_path.stat().st_size,
            "sha256": _sha256(combined_path),
            "rows": len(combined),
            "columns": list(combined.columns),
        },
        "geo_series_links_output": {
            "path": str(geo_series_links_path),
            "size_bytes": geo_series_links_path.stat().st_size,
            "sha256": _sha256(geo_series_links_path),
            "rows": len(geo_series_links),
            "columns": list(geo_series_links.columns),
        },
        "geo_linked_studies_without_series_link": missing_geo_series_links,
        "duration_sec": round(time.monotonic() - started, 3),
    }
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return provenance


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse the highest-release ImmPort Tab source for every cached study."
    )
    parser.add_argument(
        "--cache-root",
        default=str(DEFAULT_CACHE_ROOT),
        help="Root containing dynamically discovered ImmPort study directories",
    )
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        provenance = run_batch(args.cache_root)
    except (ImmportBatchParseError, OSError, ValueError) as exc:
        logger.error("ImmPort batch parsing failed: %s", exc)
        return 1

    counts = provenance["counts"]
    logger.info(
        "Parsed %s/%s studies into %s",
        counts["studies_succeeded"],
        counts["studies_discovered"],
        provenance["combined_output"]["path"],
    )
    return 0 if provenance["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
