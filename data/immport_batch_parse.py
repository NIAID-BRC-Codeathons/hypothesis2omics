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
    containing all successful studies and ``parsed/batch_manifest.json`` containing
    selected sources, study statuses, counts, errors, output hashes, and timing.

Failure behavior:
    A study failure is recorded and does not stop later studies. The command exits 1
    if any study fails; otherwise it exits 0. If no usable Tab source is discovered,
    ``ImmportBatchParseError`` is raised.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import re
import time
from datetime import datetime, timezone
from pathlib import Path
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

BATCH_VERSION = "0.1.0"
DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent / "immport_cache"
SOURCE_PATTERN = re.compile(
    r"^(?P<accession>SDY\d+)-DR(?P<release>\d+)_Tab(?P<zip>\.zip)?$",
    re.IGNORECASE,
)
COMBINED_FILENAME = "sample_manifest.tsv"
BATCH_PROVENANCE_FILENAME = "batch_manifest.json"


class ImmportBatchParseError(RuntimeError):
    """Raised when the cache root contains no usable ImmPort Tab sources."""


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    provenance_path = combined_dir / BATCH_PROVENANCE_FILENAME

    study_records = []
    frames = []
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
            frames.append(frame)
            record.update(
                {
                    "status": "success",
                    "experimental_samples": len(frame),
                    "repository_linked_samples": study_provenance["counts"][
                        "repository_linked_samples"
                    ],
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

    failures = sum(record["status"] == "failed" for record in study_records)
    provenance = {
        "status": "partial" if failures else "success",
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
        },
        "combined_output": {
            "path": str(combined_path),
            "size_bytes": combined_path.stat().st_size,
            "sha256": _sha256(combined_path),
            "rows": len(combined),
            "columns": list(combined.columns),
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
    return 1 if counts["studies_failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
