"""
Dataset acquisition for the Hypothesis2Omics pipeline (ImmPort).

For each ImmPort Study Accession (SDYxxxx), this module authenticates against the
ImmPort Shared Data API, retrieves the study manifest, resolves each file's DRS ID
to a signed download URL, caches the artifacts under data/immport_cache, and records
per-artifact provenance (mirrors the structure of geo_fetch_module.py).

Unlike GEO, ImmPort requires an API key to resolve authenticated download URLs.
See:
  https://docs.immport.org/apidocumentation/immport-auth-service/apikeys/
  https://docs.immport.org/download/guide/

------------------------------------------------------------------------------
CREDENTIALS — READ BEFORE RUNNING (team / codeathon usage)
------------------------------------------------------------------------------
Every team member should register their own free ImmPort account and use their
OWN credentials. Do NOT share a single username/password or API key across the
team, and do NOT commit credentials to the repo or hardcode them in this file.

Register:    https://www.immport.org
Get an API key (recommended): https://www.immport.org/auth/api/keys
    -> select scopes: "browse" and "download"
    -> download/copy the key and treat it like a password

Use either the downloaded key file (recommended) or set the raw key as an
environment variable (never hardcode it):

    python data/immport_fetch_module.py SDY1529 \
        --api-key-file data/immport_cache/immport-key-....json

    export IMMPORT_API_KEY="..."
    python data/immport_fetch_module.py SDY1529

For live demos, use a dedicated/throwaway demo account rather than a personal
one, and rotate its key afterward.

------------------------------------------------------------------------------
How to run this program
------------------------------------------------------------------------------

Command-line usage:
    CANDIDATE_SDY_IDS is the default study list. Running the script without SDY
    arguments fetches those candidates:

    python data/immport_fetch_module.py \
        --api-key-file data/immport_cache/immport-key-....json

    To fetch only specific studies, list their accessions:

    python data/immport_fetch_module.py SDY1529 SDY1264 \
        --api-key-file data/immport_cache/immport-key-....json

Python usage:
    Import and call fetch_immport_datasets with the studies you want. The function
    returns one FetchRecord per study:

    from data.immport_fetch_module import fetch_immport_datasets

    results = fetch_immport_datasets(
        ["SDY1529", "SDY1264"],
        api_key_file="data/immport_cache/immport-key-....json",
    )

Output:
    By default, files are saved under data/immport_cache/<SDY_ID>/. The fetcher
    preserves ImmPort subdirectories such as StudyFiles/ and ResultFiles/ and saves
    <SDY_ID>_filepath_manifest.json alongside them. Use --destdir to select another
    cache root.

    Each FetchRecord includes status, source paths, local paths, file sizes, SHA-256
    checksums, errors, timing, and listed/downloaded file counts. Command-line runs
    also write data/immport_cache/manifest.json and append to
    data/immport_cache/provenance_log.jsonl.

"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import platform
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger("immport_fetcher")

DATA_DIR = Path(__file__).resolve().parent
DEFAULT_CACHE_DIR = DATA_DIR / "immport_cache"

QUERY_BASE_URL = "https://www.immport.org/data/query"
STUDY_MANIFEST_ENDPOINT = f"{QUERY_BASE_URL}/api/study/manifest"
DRS_DOWNLOAD_ENDPOINT = f"{QUERY_BASE_URL}/drs/download"

USER_AGENT = "Hypothesis2Omics/0.1 (ImmPort dataset retrieval)"
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
MAX_ERROR_BODY_CHARS = 500
SOURCE_MISSING_MARKERS = ("nosuchkey", "specified key does not exist")
CANDIDATE_SDY_IDS = ["SDY1529", "SDY1264", "SDY1294", "SDY1289"]


# --------------------------------------------------------------------------- #
# Provenance records (same shape/spirit as geo_fetch_module.py)
# --------------------------------------------------------------------------- #


@dataclass
class ArtifactRecord:
    """Provenance for one downloaded or cached ImmPort file."""

    artifact_type: str  # e.g. "result_file", "study_file", "filepath_manifest"
    status: str  # "success", "failed", "source_missing", or "cached"
    source_path: str  # ImmPort-internal path, e.g. "/SDY208/ResultFiles/..."
    local_file: str
    retrieval_tool: str
    retrieval_tool_version: str
    file_detail: Optional[str] = None
    measurement_technique: Optional[str] = None
    sha256: Optional[str] = None
    size_bytes: Optional[int] = None
    error: Optional[str] = None


@dataclass
class FetchRecord:
    """Study-level result containing per-artifact provenance."""

    sdy_id: str
    status: str  # "success", "partial", "failed", or "cached"
    run_started_at_utc: str
    destdir: str
    api_base_url: str
    filepath_manifest_url: Optional[str]
    artifacts: list[ArtifactRecord] = field(default_factory=list)
    n_files_listed: Optional[int] = None
    n_files_downloaded: Optional[int] = None
    error: Optional[str] = None
    duration_sec: Optional[float] = None

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------- #
# Authentication
# --------------------------------------------------------------------------- #


class ImmportAuthError(RuntimeError):
    """Raised when an API key is missing or cannot be loaded."""


class ImmportSession:
    """Hold an ImmPort API key loaded directly, from a JSON key file, or the environment."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        api_key_file: Optional[str] = None,
    ) -> None:
        key_file = api_key_file or os.environ.get("IMMPORT_API_KEY_FILE")
        self.api_key = api_key or (self._read_api_key_file(key_file) if key_file else None)
        self.api_key = self.api_key or os.environ.get("IMMPORT_API_KEY")

        if not self.api_key:
            raise ImmportAuthError(
                "No ImmPort API key found. Pass --api-key-file, set "
                "IMMPORT_API_KEY_FILE, or set IMMPORT_API_KEY. See the README for setup."
            )

    @staticmethod
    def _read_api_key_file(path: str) -> str:
        try:
            payload = json.loads(Path(path).expanduser().read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ImmportAuthError(f"Could not read ImmPort API key file {path!r}: {exc}") from exc
        key = payload.get("api_key") if isinstance(payload, dict) else None
        if not isinstance(key, str) or not key.strip():
            raise ImmportAuthError(f"ImmPort API key file {path!r} has no 'api_key' value")
        return key.strip()

    def auth_header(self) -> dict:
        return {"Authorization": f"Bearer {self.api_key}"}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _normalize_sdy_id(sdy_id: str) -> str:
    normalized_id = sdy_id.strip().upper()
    if not normalized_id.startswith("SDY") or not normalized_id[3:].isdigit():
        raise ValueError(f"Invalid ImmPort study accession: {sdy_id!r}")
    return normalized_id


def _authed_get(url: str, session: ImmportSession, timeout: int = 60):
    request = Request(
        url,
        headers={"User-Agent": USER_AGENT, **session.auth_header()},
    )
    return urlopen(request, timeout=timeout)


def _sha256_of_file(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _http_error_detail(exc: HTTPError, *secrets: Optional[str]) -> str:
    """Return bounded HTTP diagnostics without exposing credentials or signed URLs."""
    detail = f"HTTP {exc.code}: {exc.reason}"
    try:
        body = exc.read().decode("utf-8", errors="replace")
    except Exception:
        body = ""

    if body:
        body = " ".join(body.split())
        for secret in secrets:
            if secret:
                body = body.replace(secret, "[REDACTED]")
        detail = f"{detail}; response={body[:MAX_ERROR_BODY_CHARS]!r}"
    return detail


def _exception_detail(exc: Exception, *secrets: Optional[str]) -> str:
    """Return an exception message with known secrets removed."""
    detail = str(exc)
    for secret in secrets:
        if secret:
            detail = detail.replace(secret, "[REDACTED]")
    return detail


def _signed_url_metadata(download_url: str) -> str:
    """Return non-secret SigV4 metadata for diagnostics."""
    parsed = urlparse(download_url)
    query = parse_qs(parsed.query)

    def first_value(name: str) -> str:
        values = query.get(name)
        return values[0] if values else "missing"

    return ", ".join(
        (
            f"host={parsed.hostname or 'missing'}",
            f"signed_headers={first_value('X-Amz-SignedHeaders')}",
            f"issued_at={first_value('X-Amz-Date')}",
            f"expires_sec={first_value('X-Amz-Expires')}",
        )
    )


def _is_source_missing(detail: str) -> bool:
    """Return whether an HTTP diagnostic explicitly identifies a missing object."""
    normalized = detail.lower()
    return any(marker in normalized for marker in SOURCE_MISSING_MARKERS)


def _probe_stream_for_missing_source(
    file_uuid: str,
    session: ImmportSession,
) -> Optional[str]:
    """Use the DRS stream method to confirm whether a denied object is missing."""
    resolver_url = f"{DRS_DOWNLOAD_ENDPOINT}/stream/{file_uuid}"
    stream_url: Optional[str] = None
    try:
        with _authed_get(resolver_url, session) as response:
            payload = json.loads(response.read().decode("utf-8"))
        if not isinstance(payload, dict):
            return None
        stream_url = payload.get("url")
        if not isinstance(stream_url, str) or not stream_url:
            return None
        with urlopen(Request(stream_url), timeout=120):
            return None
    except HTTPError as exc:
        detail = _http_error_detail(exc, session.api_key, stream_url)
        return detail if _is_source_missing(detail) else None
    except Exception:
        return None


def _fetch_filepath_manifest(sdy_id: str, session: ImmportSession) -> list[dict]:
    """Retrieve the current study manifest, including DRS IDs for every file."""
    url = f"{STUDY_MANIFEST_ENDPOINT}/{sdy_id}?fileType=all&format=json"
    with _authed_get(url, session) as response:
        payload = json.loads(response.read().decode("utf-8"))
    if isinstance(payload, dict) and "content" in payload:
        payload = payload["content"]
    if not isinstance(payload, list):
        raise ValueError(f"Unexpected study manifest response for {sdy_id}: {type(payload)}")

    seen = set()
    deduped = []
    for entry in payload:
        identity = entry.get("fileUUID") or entry.get("path")
        if identity and identity not in seen:
            seen.add(identity)
            deduped.append(entry)
    return deduped


def _download_file(
    file_path: str,
    file_uuid: str,
    destination: Path,
    session: ImmportSession,
    force: bool,
    artifact_type: str,
    expected_size: Optional[int] = None,
) -> ArtifactRecord:
    """Resolve a DRS ID to a signed S3 URL and download one ImmPort file."""
    cached_size_is_valid = (
        destination.exists()
        and destination.stat().st_size > 0
        and (expected_size is None or destination.stat().st_size == expected_size)
    )
    if cached_size_is_valid and not force:
        return ArtifactRecord(
            artifact_type=artifact_type,
            status="cached",
            source_path=file_path,
            local_file=str(destination),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            sha256=_sha256_of_file(destination),
            size_bytes=destination.stat().st_size,
        )

    temporary_path = destination.with_name(f"{destination.name}.part")

    resolver_url = f"{DRS_DOWNLOAD_ENDPOINT}/s3/{file_uuid}"
    try:
        with _authed_get(resolver_url, session) as response:
            resolver_payload = json.loads(response.read().decode("utf-8"))
        if not isinstance(resolver_payload, dict):
            raise ValueError("DRS resolver response is not a JSON object")
        download_url = resolver_payload.get("url")
        if not isinstance(download_url, str) or not download_url:
            raise ValueError("DRS response is missing a download URL")
    except HTTPError as exc:
        if temporary_path.exists():
            temporary_path.unlink()
        return ArtifactRecord(
            artifact_type=artifact_type,
            status="failed",
            source_path=file_path,
            local_file=str(destination),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            error=f"DRS resolver failed: {_http_error_detail(exc, session.api_key)}",
        )
    except Exception as exc:
        if temporary_path.exists():
            temporary_path.unlink()
        return ArtifactRecord(
            artifact_type=artifact_type,
            status="failed",
            source_path=file_path,
            local_file=str(destination),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            error=f"DRS resolver failed: {_exception_detail(exc, session.api_key)}",
        )

    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        request = Request(download_url, headers={"User-Agent": USER_AGENT})
        with urlopen(request, timeout=120) as response, temporary_path.open("wb") as output:
            while True:
                chunk = response.read(DOWNLOAD_CHUNK_SIZE)
                if not chunk:
                    break
                output.write(chunk)

        if not temporary_path.is_file() or temporary_path.stat().st_size == 0:
            raise ValueError("Downloaded file is missing or empty")
        if expected_size is not None and temporary_path.stat().st_size != expected_size:
            raise ValueError(
                f"Downloaded size {temporary_path.stat().st_size} does not match "
                f"manifest size {expected_size}"
            )

        os.replace(temporary_path, destination)
        return ArtifactRecord(
            artifact_type=artifact_type,
            status="success",
            source_path=file_path,
            local_file=str(destination),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            sha256=_sha256_of_file(destination),
            size_bytes=destination.stat().st_size,
        )
    except HTTPError as exc:
        if temporary_path.exists():
            temporary_path.unlink()
        detail = _http_error_detail(exc, session.api_key, download_url)
        status = "source_missing" if _is_source_missing(detail) else "failed"
        if status == "failed" and exc.code == 403:
            stream_detail = _probe_stream_for_missing_source(file_uuid, session)
            if stream_detail:
                status = "source_missing"
                detail = f"{detail}; stream probe={stream_detail}"
        return ArtifactRecord(
            artifact_type=artifact_type,
            status=status,
            source_path=file_path,
            local_file=str(destination),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            error=(
                f"Signed download failed: {detail}; "
                f"signed_url_metadata=({_signed_url_metadata(download_url)})"
            ),
        )
    except Exception as exc:
        if temporary_path.exists():
            temporary_path.unlink()
        return ArtifactRecord(
            artifact_type=artifact_type,
            status="failed",
            source_path=file_path,
            local_file=str(destination),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            error=(
                "Signed download failed: "
                f"{_exception_detail(exc, session.api_key, download_url)}"
            ),
        )


def _summarize_status(artifacts: Sequence[ArtifactRecord]) -> str:
    statuses = {artifact.status for artifact in artifacts}
    if not statuses:
        return "failed"
    if statuses == {"cached"}:
        return "cached"
    if statuses <= {"success", "cached"}:
        return "success"
    if statuses <= {"failed", "source_missing"}:
        return "failed"
    return "partial"


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #


def fetch_one(
    sdy_id: str,
    destdir: Path,
    session: ImmportSession,
    force: bool = False,
    max_files: Optional[int] = None,
) -> FetchRecord:
    """Fetch the file-path manifest and every listed file for one ImmPort study."""

    start = time.time()
    now = datetime.now(timezone.utc).isoformat()

    try:
        normalized_id = _normalize_sdy_id(sdy_id)
    except ValueError as exc:
        return FetchRecord(
            sdy_id=sdy_id,
            status="failed",
            run_started_at_utc=now,
            destdir=str(destdir),
            api_base_url=QUERY_BASE_URL,
            filepath_manifest_url=None,
            error=str(exc),
            duration_sec=round(time.time() - start, 3),
        )

    sdy_dir = destdir / normalized_id
    sdy_dir.mkdir(parents=True, exist_ok=True)
    manifest_url = f"{STUDY_MANIFEST_ENDPOINT}/{normalized_id}?fileType=all&format=json"

    try:
        manifest_entries = _fetch_filepath_manifest(normalized_id, session)
    except (HTTPError, URLError, ValueError, json.JSONDecodeError) as exc:
        return FetchRecord(
            sdy_id=normalized_id,
            status="failed",
            run_started_at_utc=now,
            destdir=str(sdy_dir),
            api_base_url=QUERY_BASE_URL,
            filepath_manifest_url=manifest_url,
            error=f"filePath manifest retrieval failed: {exc}",
            duration_sec=round(time.time() - start, 3),
        )

    # Persist the raw manifest for provenance/reproducibility.
    manifest_path = sdy_dir / f"{normalized_id}_filepath_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest_entries, handle, indent=2)

    artifacts: list[ArtifactRecord] = [
        ArtifactRecord(
            artifact_type="filepath_manifest",
            status="success",
            source_path=manifest_url,
            local_file=str(manifest_path),
            retrieval_tool="urllib.request",
            retrieval_tool_version=platform.python_version(),
            size_bytes=manifest_path.stat().st_size,
            sha256=_sha256_of_file(manifest_path),
        )
    ]

    entries_to_fetch = manifest_entries[:max_files] if max_files is not None else manifest_entries
    errors: list[str] = []

    for entry in entries_to_fetch:
        file_path = entry.get("path")
        file_uuid = entry.get("fileUUID")
        if not file_path or not file_uuid:
            errors.append(f"Manifest entry is missing path or fileUUID: {entry!r}")
            continue
        relative_parts = PurePosixPath(file_path.lstrip("/")).parts
        if relative_parts and relative_parts[0] == normalized_id:
            relative_parts = relative_parts[1:]
        if not relative_parts or ".." in relative_parts:
            errors.append(f"Manifest entry has an unsafe path: {file_path!r}")
            continue
        destination = sdy_dir.joinpath(*relative_parts)
        artifact = _download_file(
            file_path,
            file_uuid,
            destination,
            session,
            force,
            artifact_type=entry.get("fileType") or "study_file",
            expected_size=entry.get("filesizeBytes"),
        )
        artifacts.append(artifact)
        if artifact.status in {"failed", "source_missing"}:
            errors.append(f"{Path(file_path).name}: {artifact.error}")

    n_downloaded = sum(
        1
        for a in artifacts
        if a.artifact_type != "filepath_manifest" and a.status in {"success", "cached"}
    )
    status = _summarize_status(artifacts)
    if len(entries_to_fetch) < len(manifest_entries) and status in {"success", "cached"}:
        status = "partial"
    if errors and status in {"success", "cached"}:
        status = "partial"

    return FetchRecord(
        sdy_id=normalized_id,
        status=status,
        run_started_at_utc=now,
        destdir=str(sdy_dir),
        api_base_url=QUERY_BASE_URL,
        filepath_manifest_url=manifest_url,
        artifacts=artifacts,
        n_files_listed=len(manifest_entries),
        n_files_downloaded=n_downloaded,
        error="; ".join(errors) or None,
        duration_sec=round(time.time() - start, 3),
    )


def fetch_immport_datasets(
    sdy_ids: Sequence[str],
    destdir: str = str(DEFAULT_CACHE_DIR),
    api_key: Optional[str] = None,
    api_key_file: Optional[str] = None,
    force: bool = False,
    max_files_per_study: Optional[int] = None,
    provenance_log_path: Optional[str] = None,
) -> list[FetchRecord]:
    """
    Fetch file-path manifests and files for each SDY accession, with provenance logging.

    The API key resolves from arguments first, then IMMPORT_API_KEY_FILE or
    IMMPORT_API_KEY. Use your own ImmPort account and do not share the key.
    """

    session = ImmportSession(api_key=api_key, api_key_file=api_key_file)

    base = Path(destdir)
    base.mkdir(parents=True, exist_ok=True)
    records: list[FetchRecord] = []
    provenance_path = Path(provenance_log_path) if provenance_log_path else None
    if provenance_path:
        provenance_path.parent.mkdir(parents=True, exist_ok=True)

    for sdy_id in sdy_ids:
        record = fetch_one(sdy_id, base, session, force=force, max_files=max_files_per_study)
        records.append(record)
        if provenance_path:
            with provenance_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record.to_dict()) + "\n")

    _print_summary(records)
    return records


def _print_summary(records: Sequence[FetchRecord]) -> None:
    logger.info("Fetch summary")
    for record in records:
        logger.info(
            "[%s] %s — %s/%s files (%ss)",
            record.status.upper(),
            record.sdy_id,
            record.n_files_downloaded,
            record.n_files_listed,
            record.duration_sec,
        )
        if record.error:
            logger.info("    %s", record.error)
    available = sum(1 for record in records if record.status in {"success", "cached"})
    logger.info("%s/%s studies completely available in cache.", available, len(records))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch ImmPort study file manifests and result files with provenance. "
        "Requires your own ImmPort API key; see the README."
    )
    parser.add_argument(
        "sdy_ids",
        nargs="*",
        default=CANDIDATE_SDY_IDS,
        help="ImmPort study accessions; defaults to the project candidate list.",
    )
    parser.add_argument(
        "--destdir",
        default=str(DEFAULT_CACHE_DIR),
        help="Cache directory (default: data/immport_cache).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Download and validate artifacts again even when cached.",
    )
    parser.add_argument(
        "--max-files-per-study",
        type=int,
        default=None,
        help="Optional cap on number of files downloaded per study (manifest is always full).",
    )
    parser.add_argument(
        "--api-key-file",
        default=None,
        help="Path to the JSON key file downloaded from ImmPort (recommended).",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    cache_dir = Path(args.destdir)
    provenance_path = cache_dir / "provenance_log.jsonl"
    try:
        results = fetch_immport_datasets(
            args.sdy_ids,
            destdir=str(cache_dir),
            api_key_file=args.api_key_file,
            force=args.force,
            max_files_per_study=args.max_files_per_study,
            provenance_log_path=str(provenance_path),
        )
    except ImmportAuthError as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc

    manifest_path = cache_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump([record.to_dict() for record in results], handle, indent=2)
    logger.info("Manifest written to %s", manifest_path)


if __name__ == "__main__":
    main()
