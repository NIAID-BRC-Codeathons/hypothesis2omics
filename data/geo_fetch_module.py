"""
Dataset acquisition for the Hypothesis2Omics pipeline (GEO).

For each GEO Series accession (GSExxxx), this module downloads the family SOFT file
and every published series-matrix file, validates and caches the compressed files
under data/geo_cache, and records per-artifact provenance. GEOparse reads the family
SOFT file to report platform and sample counts.

GEO downloads require internet access and the GEOparse package, but no account or
credentials. Install the project's dependencies with `uv sync` before running.

------------------------------------------------------------------------------
How to run this program
------------------------------------------------------------------------------

Command-line usage:
    Running the script without GSE arguments reads accessions from the structured
    ImmPort links produced by ``immport_batch_parse.py``:

    python data/geo_fetch_module.py

    To fetch only specific studies, list their accessions:

    python data/geo_fetch_module.py GSE13699 GSE125921

    Use --destdir to select another cache directory. Existing valid files are reused;
    use --force to download and validate them again.

Python usage:
    Import and call fetch_geo_datasets with the studies you want. The function returns
    one FetchRecord per study:

    from data.geo_fetch_module import fetch_geo_datasets

    results = fetch_geo_datasets(
        ["GSE13699", "GSE125921"],
        provenance_log_path="data/geo_cache/provenance_log.jsonl",
    )

Output:
    By default, files are saved under data/geo_cache/<GSE_ID>/. Each FetchRecord
    includes status, source URLs, local paths, file sizes, SHA-256 checksums, errors,
    timing, and parsed platform/sample counts. Command-line runs also write
    data/geo_cache/manifest.json and append to data/geo_cache/provenance_log.jsonl.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import logging
import os
import platform
import re
import shutil
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from typing import Optional
from urllib.parse import unquote, urljoin, urlparse
from urllib.request import Request, urlopen

try:
    import GEOparse
except ImportError:
    GEOparse = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("geo_fetcher")

DATA_DIR = Path(__file__).resolve().parent
DEFAULT_CACHE_DIR = DATA_DIR / "geo_cache"
DEFAULT_GEO_LINKS_PATH = DATA_DIR / "immport_cache" / "parsed" / "geo_series_links.tsv"
SELECTION_PROVENANCE_FILENAME = "geo_fetch_selection.provenance.json"
FETCHER_VERSION = "0.2.0"
NCBI_GEO_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series"
USER_AGENT = "Hypothesis2Omics/0.1 (GEO dataset retrieval)"
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
GSE_PATTERN = re.compile(r"^GSE\d+$", re.IGNORECASE)


class GeoSelectionError(RuntimeError):
    """Raised when GEO accessions cannot be selected from parsed ImmPort links."""


@dataclass
class ArtifactRecord:
    """Provenance for one downloaded or cached GEO artifact."""

    artifact_type: str
    status: str  # "success", "failed", or "cached"
    source_url: str
    local_file: str
    retrieval_tool: str
    retrieval_tool_version: str
    sha256: Optional[str] = None
    size_bytes: Optional[int] = None
    error: Optional[str] = None


@dataclass
class FetchRecord:
    """Dataset-level result containing per-artifact provenance."""

    gse_id: str
    status: str  # "success", "partial", "failed", or "cached"
    run_started_at_utc: str
    destdir: str
    geoparse_version: Optional[str]
    family_soft_url: Optional[str]
    matrix_index_url: Optional[str]
    artifacts: list[ArtifactRecord] = field(default_factory=list)
    n_platforms: Optional[int] = None
    n_samples: Optional[int] = None
    error: Optional[str] = None
    duration_sec: Optional[float] = None

    def to_dict(self) -> dict:
        return asdict(self)


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list) -> None:
        if tag.lower() != "a":
            return
        for name, value in attrs:
            if name.lower() == "href" and value:
                self.hrefs.append(value)


def _normalize_gse_id(gse_id: str) -> str:
    normalized_id = gse_id.strip().upper()
    if not GSE_PATTERN.fullmatch(normalized_id):
        raise ValueError(f"Invalid GEO Series accession: {gse_id!r}")
    return normalized_id


def load_gse_ids_from_parsed_links(path: str | Path) -> list[str]:
    """Load unique GSE accessions from an ImmPort ``geo_series_links.tsv`` file."""

    links_path = Path(path)
    if not links_path.is_file():
        raise GeoSelectionError(
            f"Parsed GEO links file does not exist: {links_path}. "
            "Run data/immport_batch_parse.py first."
        )

    try:
        with links_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise GeoSelectionError(f"Parsed GEO links file is empty: {links_path}")
            normalized_headers = [header.strip().lower() for header in reader.fieldnames]
            if len(normalized_headers) != len(set(normalized_headers)):
                raise GeoSelectionError(
                    f"Parsed GEO links file has duplicate normalized columns: {links_path}"
                )
            if "gse_accession" not in normalized_headers:
                raise GeoSelectionError(
                    f"Parsed GEO links file is missing gse_accession: {links_path}"
                )
            accession_header = reader.fieldnames[normalized_headers.index("gse_accession")]
            raw_accessions = [row.get(accession_header, "") for row in reader]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise GeoSelectionError(
            f"Could not read parsed GEO links file {links_path}: {exc}"
        ) from exc

    accessions: list[str] = []
    seen: set[str] = set()
    for raw_accession in raw_accessions:
        try:
            accession = _normalize_gse_id(raw_accession or "")
        except ValueError as exc:
            raise GeoSelectionError(f"Invalid gse_accession in {links_path}: {exc}") from exc
        if accession not in seen:
            seen.add(accession)
            accessions.append(accession)

    if not accessions:
        raise GeoSelectionError(f"Parsed GEO links file contains no GSE accessions: {links_path}")
    return accessions


def _gse_range_subdir(gse_id: str) -> str:
    return re.sub(r"\d{1,3}$", "nnn", gse_id)


def _gse_directory_url(gse_id: str) -> str:
    return f"{NCBI_GEO_BASE_URL}/{_gse_range_subdir(gse_id)}/{gse_id}"


def _gse_family_soft_url(gse_id: str) -> str:
    normalized_id = _normalize_gse_id(gse_id)
    return f"{_gse_directory_url(normalized_id)}/soft/{normalized_id}_family.soft.gz"


def _gse_matrix_directory_url(gse_id: str) -> str:
    normalized_id = _normalize_gse_id(gse_id)
    return f"{_gse_directory_url(normalized_id)}/matrix/"


def _installed_geoparse_version() -> Optional[str]:
    try:
        return importlib.metadata.version("GEOparse")
    except importlib.metadata.PackageNotFoundError:
        return None


def _open_url(url: str, timeout: int = 60):
    request = Request(url, headers={"User-Agent": USER_AGENT})
    return urlopen(request, timeout=timeout)


def _discover_series_matrix_urls(gse_id: str) -> list[str]:
    normalized_id = _normalize_gse_id(gse_id)
    directory_url = _gse_matrix_directory_url(normalized_id)
    with _open_url(directory_url) as response:
        listing = response.read().decode("utf-8", errors="replace")

    parser = _LinkParser()
    parser.feed(listing)
    filename_pattern = re.compile(
        rf"^{re.escape(normalized_id)}(?:-[^/]+)?_series_matrix\.txt\.gz$"
    )
    directory_path = urlparse(directory_url).path
    urls = set()

    for href in parser.hrefs:
        candidate_url = urljoin(directory_url, href)
        parsed = urlparse(candidate_url)
        filename = unquote(Path(parsed.path).name)
        if parsed.netloc != "ftp.ncbi.nlm.nih.gov":
            continue
        if not parsed.path.startswith(directory_path):
            continue
        if filename_pattern.fullmatch(filename):
            urls.add(candidate_url)

    return sorted(urls)


def _sha256_of_file(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_gzip(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Downloaded file is missing or empty: {path}")
    with gzip.open(path, "rb") as handle:
        for _ in iter(lambda: handle.read(DOWNLOAD_CHUNK_SIZE), b""):
            pass


def _artifact_record(
    artifact_type: str,
    status: str,
    source_url: str,
    destination: Path,
    error: Optional[str] = None,
) -> ArtifactRecord:
    exists = destination.is_file()
    return ArtifactRecord(
        artifact_type=artifact_type,
        status=status,
        source_url=source_url,
        local_file=str(destination),
        retrieval_tool="urllib.request",
        retrieval_tool_version=platform.python_version(),
        sha256=_sha256_of_file(destination) if exists else None,
        size_bytes=destination.stat().st_size if exists else None,
        error=error,
    )


def _fetch_gzip_artifact(
    source_url: str,
    destination: Path,
    artifact_type: str,
    force: bool,
) -> ArtifactRecord:
    if destination.exists() and not force:
        try:
            _validate_gzip(destination)
            return _artifact_record(artifact_type, "cached", source_url, destination)
        except (OSError, ValueError) as exc:
            logger.warning("Invalid cached file %s; downloading again: %s", destination, exc)

    temporary_path = destination.with_name(f"{destination.name}.part")
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with _open_url(source_url) as response, temporary_path.open("wb") as output:
            shutil.copyfileobj(response, output, length=DOWNLOAD_CHUNK_SIZE)
        _validate_gzip(temporary_path)
        os.replace(temporary_path, destination)
        return _artifact_record(artifact_type, "success", source_url, destination)
    except Exception as exc:
        if temporary_path.exists():
            temporary_path.unlink()
        return _artifact_record(
            artifact_type,
            "failed",
            source_url,
            destination,
            error=str(exc),
        )


def _summarize_status(artifacts: Sequence[ArtifactRecord]) -> str:
    statuses = {artifact.status for artifact in artifacts}
    if statuses == {"cached"}:
        return "cached"
    if "failed" not in statuses:
        return "success"
    if statuses == {"failed"}:
        return "failed"
    return "partial"


def fetch_one(
    gse_id: str,
    destdir: Path,
    force: bool = False,
    silent: bool = True,
) -> FetchRecord:
    """Fetch one GSE family SOFT file and all available series matrices."""

    start = time.time()
    now = datetime.now(timezone.utc).isoformat()
    artifacts: list[ArtifactRecord] = []
    errors: list[str] = []

    try:
        normalized_id = _normalize_gse_id(gse_id)
    except ValueError as exc:
        return FetchRecord(
            gse_id=gse_id,
            status="failed",
            run_started_at_utc=now,
            destdir=str(destdir),
            geoparse_version=_installed_geoparse_version(),
            family_soft_url=None,
            matrix_index_url=None,
            error=str(exc),
            duration_sec=round(time.time() - start, 3),
        )

    gse_dir = destdir / normalized_id
    gse_dir.mkdir(parents=True, exist_ok=True)
    geoparse_version = _installed_geoparse_version()
    soft_url = _gse_family_soft_url(normalized_id)
    matrix_index_url = _gse_matrix_directory_url(normalized_id)

    if GEOparse is None:
        error = "GEOparse is not installed. Install the project's pinned dependencies."
        return FetchRecord(
            gse_id=normalized_id,
            status="failed",
            run_started_at_utc=now,
            destdir=str(gse_dir),
            geoparse_version=geoparse_version,
            family_soft_url=soft_url,
            matrix_index_url=matrix_index_url,
            error=error,
            duration_sec=round(time.time() - start, 3),
        )

    soft_path = gse_dir / f"{normalized_id}_family.soft.gz"
    soft_artifact = _fetch_gzip_artifact(soft_url, soft_path, "family_soft", force)
    artifacts.append(soft_artifact)

    try:
        matrix_urls = _discover_series_matrix_urls(normalized_id)
    except Exception as exc:
        matrix_urls = []
        errors.append(f"Series-matrix discovery failed: {exc}")

    if not matrix_urls and not errors:
        errors.append("No series-matrix files were published for this accession.")

    for matrix_url in matrix_urls:
        filename = unquote(Path(urlparse(matrix_url).path).name)
        artifacts.append(
            _fetch_gzip_artifact(
                matrix_url,
                gse_dir / filename,
                "series_matrix",
                force,
            )
        )

    failed_artifacts = [artifact for artifact in artifacts if artifact.status == "failed"]
    errors.extend(
        f"{artifact.artifact_type} failed: {artifact.error}" for artifact in failed_artifacts
    )

    n_platforms = None
    n_samples = None
    if soft_artifact.status != "failed":
        try:
            gse = GEOparse.get_GEO(filepath=str(soft_path), silent=silent)
            n_platforms = len(gse.gpls)
            n_samples = len(gse.gsms)
        except Exception as exc:
            errors.append(f"GEOparse could not parse the family SOFT file: {exc}")

    status = _summarize_status(artifacts)
    if errors and status in {"success", "cached"}:
        status = "partial"

    return FetchRecord(
        gse_id=normalized_id,
        status=status,
        run_started_at_utc=now,
        destdir=str(gse_dir),
        geoparse_version=geoparse_version,
        family_soft_url=soft_url,
        matrix_index_url=matrix_index_url,
        artifacts=artifacts,
        n_platforms=n_platforms,
        n_samples=n_samples,
        error="; ".join(errors) or None,
        duration_sec=round(time.time() - start, 3),
    )


def fetch_geo_datasets(
    gse_ids: Sequence[str],
    destdir: str = str(DEFAULT_CACHE_DIR),
    force: bool = False,
    provenance_log_path: Optional[str] = None,
) -> list[FetchRecord]:
    """Fetch family SOFT and series-matrix artifacts for each GSE accession."""

    base = Path(destdir)
    base.mkdir(parents=True, exist_ok=True)
    records: list[FetchRecord] = []
    provenance_path = Path(provenance_log_path) if provenance_log_path else None
    if provenance_path:
        provenance_path.parent.mkdir(parents=True, exist_ok=True)

    for gse_id in gse_ids:
        record = fetch_one(gse_id, base, force=force)
        records.append(record)
        if provenance_path:
            with provenance_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record.to_dict()) + "\n")

    _print_summary(records)
    return records


def fetch_series_matrices(
    gse_ids: Sequence[str],
    destdir: str = str(DEFAULT_CACHE_DIR),
    force: bool = False,
    provenance_log_path: Optional[str] = None,
) -> list[FetchRecord]:
    """Backward-compatible alias; fetches both SOFT and series-matrix artifacts."""

    return fetch_geo_datasets(gse_ids, destdir, force, provenance_log_path)


def _print_summary(records: Sequence[FetchRecord]) -> None:
    logger.info("Fetch summary")
    for record in records:
        logger.info("[%s] %s (%ss)", record.status.upper(), record.gse_id, record.duration_sec)
        if record.error:
            logger.info("    %s", record.error)
    available = sum(1 for record in records if record.status in {"success", "cached"})
    logger.info("%s/%s datasets completely available in cache.", available, len(records))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch GEO family SOFT and series-matrix files with provenance."
    )
    parser.add_argument(
        "gse_ids",
        nargs="*",
        help=(
            "GSE accessions; when omitted, accessions are read from "
            "data/immport_cache/parsed/geo_series_links.tsv."
        ),
    )
    parser.add_argument(
        "--links-file",
        default=str(DEFAULT_GEO_LINKS_PATH),
        help=(
            "Parsed ImmPort GEO-series link table used when no GSE accessions are given "
            "(default: data/immport_cache/parsed/geo_series_links.tsv)."
        ),
    )
    parser.add_argument(
        "--destdir",
        default=str(DEFAULT_CACHE_DIR),
        help="Cache directory (default: data/geo_cache).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Download and validate artifacts again even when cached.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    cache_dir = Path(args.destdir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    if args.gse_ids:
        selected_gse_ids = [_normalize_gse_id(gse_id) for gse_id in args.gse_ids]
        selection_mode = "command_line"
        selection_input = {"gse_ids": args.gse_ids}
    else:
        links_path = Path(args.links_file)
        selected_gse_ids = load_gse_ids_from_parsed_links(links_path)
        selection_mode = "parsed_immport_geo_series_links"
        selection_input = {
            "path": str(links_path),
            "sha256": _sha256_of_file(links_path),
        }

    provenance_path = cache_dir / "provenance_log.jsonl"
    results = fetch_geo_datasets(
        selected_gse_ids,
        destdir=str(cache_dir),
        force=args.force,
        provenance_log_path=str(provenance_path),
    )

    manifest_path = cache_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump([record.to_dict() for record in results], handle, indent=2)
    logger.info("Manifest written to %s", manifest_path)

    status_counts: dict[str, int] = {}
    for record in results:
        status_counts[record.status] = status_counts.get(record.status, 0) + 1
    selection_provenance = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool": "data.geo_fetch_module",
        "tool_version": FETCHER_VERSION,
        "selection_mode": selection_mode,
        "input": selection_input,
        "selected_gse_ids": selected_gse_ids,
        "force": args.force,
        "outputs": {
            "cache_directory": str(cache_dir),
            "manifest": str(manifest_path),
            "provenance_log": str(provenance_path),
        },
        "result_status_counts": status_counts,
    }
    selection_provenance_path = cache_dir / SELECTION_PROVENANCE_FILENAME
    with selection_provenance_path.open("w", encoding="utf-8") as handle:
        json.dump(selection_provenance, handle, indent=2)
        handle.write("\n")
    logger.info("Selection provenance written to %s", selection_provenance_path)


if __name__ == "__main__":
    main()
