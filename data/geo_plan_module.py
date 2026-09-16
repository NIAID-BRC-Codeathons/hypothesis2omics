"""Resolve ImmPort-linked GEO samples into a bounded GEO download plan.

This module maps GEO GSM accessions from an ImmPort sample manifest to their parent
GSE series and GPL platforms using NCBI E-utilities. It groups samples into candidate
downloads and marks explicit SuperSeries records as ``skip_superseries``. It retrieves
metadata only and never downloads expression matrices or raw archives.

The planner requires pandas and network access for uncached accessions. Install the
project's pinned dependencies from the repository root with ``uv sync``. No NCBI
account is required; an optional ``NCBI_API_KEY`` environment variable may be used.

------------------------------------------------------------------------------
How to run this program
------------------------------------------------------------------------------

Command-line usage:
    Build a plan from the combined ImmPort sample manifest:

    uv run python data/geo_plan_module.py \
        data/immport_cache/parsed/sample_manifest.tsv \
        --output-dir data/geo_cache/plan

    Existing GSM and GSE metadata caches are reused. To refresh them:

    uv run python data/geo_plan_module.py \
        data/immport_cache/parsed/sample_manifest.tsv \
        --output-dir data/geo_cache/plan \
        --force

    ``--batch-size`` controls accessions per NCBI query and ``--timeout`` controls
    each request timeout. Run the script with ``--help`` for all options.

Python usage:
    Import and call ``plan_geo_downloads``. It writes the plan and returns its
    provenance dictionary:

    from data.geo_plan_module import plan_geo_downloads

    provenance = plan_geo_downloads(
        "data/immport_cache/parsed/sample_manifest.tsv",
        "data/geo_cache/plan",
    )

Input:
    The input is the combined TSV written by ``immport_batch_parse.py``. Rows with
    ``repository_name`` equal to GEO and a nonempty ``repository_accession`` are
    considered. Valid ``GSM<number>`` accessions are queried; malformed accessions are
    retained with ``invalid_accession`` status.

Output:
    ``geo_download_plan.tsv`` groups requested GSMs by ImmPort experiment, GSE, and
    GPL. It includes series metadata, sample counts, resolution status, and a
    ``download``, ``skip_superseries``, or ``review`` recommendation.

    ``gsm_resolution_cache.json`` and ``gse_metadata_cache.json`` allow repeat runs
    without duplicate NCBI requests. ``geo_download_plan.provenance.json`` records
    source and output hashes, parameters, query batches, cache details, statuses,
    counts, versions, and timing. API keys are never written to provenance.

Failure behavior:
    Failed query batches are recorded as ``query_failed`` and processing continues.
    Unresolved metadata remains visible for review. The command exits 0 only when all
    GSM and GSE metadata resolves; partial or failed plans exit 1.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import platform
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

logger = logging.getLogger("geo_planner")

PLANNER_VERSION = "0.1.0"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
USER_AGENT = "Hypothesis2Omics/0.1 (GEO metadata planning)"
DEFAULT_BATCH_SIZE = 100
DEFAULT_TIMEOUT = 60
REQUEST_INTERVAL_SEC = 0.34
CACHE_SCHEMA_VERSION = 1
PLAN_FILENAME = "geo_download_plan.tsv"
CACHE_FILENAME = "gsm_resolution_cache.json"
SERIES_CACHE_FILENAME = "gse_metadata_cache.json"
PROVENANCE_FILENAME = "geo_download_plan.provenance.json"
GSM_PATTERN = re.compile(r"^GSM\d+$", re.IGNORECASE)

REQUIRED_INPUT_COLUMNS = [
    "study_accession",
    "experiment_accession",
    "experiment_name",
    "measurement_technique",
    "repository_name",
    "repository_accession",
]

PLAN_COLUMNS = [
    "study_accession",
    "experiment_accession",
    "experiment_name",
    "measurement_technique",
    "gse_accession",
    "gpl_accession",
    "series_title",
    "series_type",
    "series_sample_count",
    "series_metadata_status",
    "download_recommendation",
    "resolution_status",
    "requested_sample_count",
    "gsm_accessions",
]


class GeoPlanError(RuntimeError):
    """Raised when the input or cached GEO metadata is invalid."""


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _post_json(url: str, fields: dict[str, str], timeout: int) -> dict[str, Any]:
    request = Request(
        url,
        data=urlencode(fields).encode("utf-8"),
        headers={
            "User-Agent": USER_AGENT,
            "Content-Type": "application/x-www-form-urlencoded",
        },
        method="POST",
    )
    with urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read().decode("utf-8"))
    time.sleep(REQUEST_INTERVAL_SEC)
    return payload


def _accession_values(value: Any, prefix: str) -> list[str]:
    numbers = re.findall(rf"(?:{prefix})?(\d+)", str(value), flags=re.IGNORECASE)
    return sorted({f"{prefix}{number}" for number in numbers}, key=lambda item: int(item[3:]))


def _resolve_batch(
    accessions: list[str],
    timeout: int,
    api_key: str | None,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    query_started = _utc_now()
    terms = " OR ".join(f"{accession}[ACCN]" for accession in accessions)
    search_fields = {
        "db": "gds",
        "term": f"({terms}) AND gsm[ETYP]",
        "retmode": "json",
        "retmax": str(len(accessions)),
        "tool": "Hypothesis2Omics",
    }
    if api_key:
        search_fields["api_key"] = api_key
    search_payload = _post_json(ESEARCH_URL, search_fields, timeout)
    identifiers = search_payload.get("esearchresult", {}).get("idlist", [])
    if not isinstance(identifiers, list):
        raise GeoPlanError("NCBI ESearch response has no valid idlist")

    summaries: dict[str, Any] = {}
    if identifiers:
        summary_fields = {
            "db": "gds",
            "id": ",".join(str(identifier) for identifier in identifiers),
            "retmode": "json",
            "tool": "Hypothesis2Omics",
        }
        if api_key:
            summary_fields["api_key"] = api_key
        summary_payload = _post_json(ESUMMARY_URL, summary_fields, timeout)
        result = summary_payload.get("result", {})
        if not isinstance(result, dict):
            raise GeoPlanError("NCBI ESummary response has no valid result")
        for identifier in result.get("uids", []):
            summary = result.get(str(identifier), {})
            accession = str(summary.get("accession", "")).strip().upper()
            if str(summary.get("entrytype", "")).upper() == "GSM" and accession in accessions:
                summaries[accession] = summary

    resolved_at = _utc_now()
    records = {}
    for accession in accessions:
        summary = summaries.get(accession)
        gse_accessions = _accession_values(summary.get("gse", ""), "GSE") if summary else []
        gpl_accessions = _accession_values(summary.get("gpl", ""), "GPL") if summary else []
        records[accession] = {
            "status": "resolved" if gse_accessions else "unresolved",
            "gse_accessions": gse_accessions,
            "gpl_accessions": gpl_accessions,
            "sample_title": str(summary.get("title", "")) if summary else "",
            "resolved_at_utc": resolved_at,
            "source": "NCBI GDS ESearch/ESummary",
        }

    query_record = {
        "query_type": "gsm_resolution",
        "started_at_utc": query_started,
        "search_endpoint": ESEARCH_URL,
        "summary_endpoint": ESUMMARY_URL,
        "requested_accessions": accessions,
        "returned_uids": len(identifiers),
        "status": "success",
    }
    return records, query_record


def _resolve_series_batch(
    accessions: list[str],
    timeout: int,
    api_key: str | None,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Resolve GSE summaries and classify explicit NCBI SuperSeries records."""
    query_started = _utc_now()
    terms = " OR ".join(f"{accession}[ACCN]" for accession in accessions)
    search_fields = {
        "db": "gds",
        "term": f"({terms}) AND gse[ETYP]",
        "retmode": "json",
        "retmax": str(len(accessions)),
        "tool": "Hypothesis2Omics",
    }
    if api_key:
        search_fields["api_key"] = api_key
    search_payload = _post_json(ESEARCH_URL, search_fields, timeout)
    identifiers = search_payload.get("esearchresult", {}).get("idlist", [])
    if not isinstance(identifiers, list):
        raise GeoPlanError("NCBI GSE ESearch response has no valid idlist")

    summaries: dict[str, Any] = {}
    if identifiers:
        summary_fields = {
            "db": "gds",
            "id": ",".join(str(identifier) for identifier in identifiers),
            "retmode": "json",
            "tool": "Hypothesis2Omics",
        }
        if api_key:
            summary_fields["api_key"] = api_key
        summary_payload = _post_json(ESUMMARY_URL, summary_fields, timeout)
        result = summary_payload.get("result", {})
        if not isinstance(result, dict):
            raise GeoPlanError("NCBI GSE ESummary response has no valid result")
        for identifier in result.get("uids", []):
            summary = result.get(str(identifier), {})
            accession = str(summary.get("accession", "")).strip().upper()
            if str(summary.get("entrytype", "")).upper() == "GSE" and accession in accessions:
                summaries[accession] = summary

    resolved_at = _utc_now()
    records = {}
    for accession in accessions:
        summary = summaries.get(accession)
        if summary is None:
            records[accession] = {
                "status": "unresolved",
                "series_title": "",
                "series_type": "unknown",
                "series_sample_count": None,
                "resolved_at_utc": resolved_at,
                "source": "NCBI GDS ESearch/ESummary",
            }
            continue

        summary_text = " ".join(str(summary.get("summary", "")).casefold().split())
        series_type = (
            "superseries" if "this superseries is composed of" in summary_text else "series"
        )
        raw_sample_count = summary.get("n_samples")
        try:
            sample_count = int(raw_sample_count)
        except (TypeError, ValueError):
            sample_count = None
        records[accession] = {
            "status": "resolved",
            "series_title": str(summary.get("title", "")),
            "series_type": series_type,
            "series_sample_count": sample_count,
            "resolved_at_utc": resolved_at,
            "source": "NCBI GDS ESearch/ESummary",
        }

    query_record = {
        "query_type": "gse_metadata",
        "started_at_utc": query_started,
        "search_endpoint": ESEARCH_URL,
        "summary_endpoint": ESUMMARY_URL,
        "requested_accessions": accessions,
        "returned_uids": len(identifiers),
        "status": "success",
    }
    return records, query_record


def _load_cache(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GeoPlanError(f"Could not read GEO metadata cache {path}: {exc}") from exc
    if payload.get("schema_version") != CACHE_SCHEMA_VERSION:
        raise GeoPlanError(f"Unsupported GEO metadata cache schema in {path}")
    records = payload.get("records")
    if not isinstance(records, dict):
        raise GeoPlanError(f"GEO metadata cache has no records object: {path}")
    return records


def _write_cache(path: Path, records: dict[str, dict[str, Any]]) -> None:
    payload = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "updated_at_utc": _utc_now(),
        "records": dict(sorted(records.items())),
    }
    temporary = path.with_name(f"{path.name}.part")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _read_manifest(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise GeoPlanError(f"ImmPort sample manifest does not exist: {path}")
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise GeoPlanError(f"Could not read {path}: {exc}") from exc
    frame.columns = [str(column).strip().lower() for column in frame.columns]
    missing = sorted(set(REQUIRED_INPUT_COLUMNS) - set(frame.columns))
    if missing:
        raise GeoPlanError(f"{path} is missing required columns: {', '.join(missing)}")
    for column in REQUIRED_INPUT_COLUMNS:
        frame[column] = frame[column].str.strip()
    return frame


def _build_plan(
    source: pd.DataFrame,
    resolutions: dict[str, dict[str, Any]],
    series_metadata: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    rows = []
    for record in source.to_dict(orient="records"):
        gsm = record["repository_accession"].upper()
        resolution = resolutions[gsm]
        gse_accessions = resolution.get("gse_accessions") or [""]
        gpl_accession = ";".join(resolution.get("gpl_accessions", []))
        for gse_accession in gse_accessions:
            series = series_metadata.get(
                gse_accession,
                {
                    "status": "not_applicable" if not gse_accession else "unresolved",
                    "series_title": "",
                    "series_type": "unknown",
                    "series_sample_count": None,
                },
            )
            if series["series_type"] == "superseries":
                recommendation = "skip_superseries"
            elif series["status"] == "resolved":
                recommendation = "download"
            else:
                recommendation = "review"
            rows.append(
                {
                    "study_accession": record["study_accession"],
                    "experiment_accession": record["experiment_accession"],
                    "experiment_name": record["experiment_name"],
                    "measurement_technique": record["measurement_technique"],
                    "gse_accession": gse_accession,
                    "gpl_accession": gpl_accession,
                    "series_title": series["series_title"],
                    "series_type": series["series_type"],
                    "series_sample_count": series["series_sample_count"],
                    "series_metadata_status": series["status"],
                    "download_recommendation": recommendation,
                    "resolution_status": resolution["status"],
                    "gsm_accession": gsm,
                }
            )

    if not rows:
        return pd.DataFrame(columns=PLAN_COLUMNS)
    detailed = pd.DataFrame(rows)
    group_columns = PLAN_COLUMNS[:-2]
    grouped = (
        detailed.groupby(group_columns, dropna=False, sort=True)["gsm_accession"]
        .agg(
            requested_sample_count="nunique",
            gsm_accessions=lambda values: ";".join(sorted(set(values))),
        )
        .reset_index()
    )
    return grouped.loc[:, PLAN_COLUMNS]


def plan_geo_downloads(
    manifest_path: str | Path,
    output_dir: str | Path,
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    force: bool = False,
) -> dict[str, Any]:
    """Resolve GEO-linked GSM accessions and write a grouped download plan."""
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if timeout < 1:
        raise ValueError("timeout must be at least 1 second")

    started = time.monotonic()
    started_at = _utc_now()
    source_path = Path(manifest_path).expanduser().resolve()
    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    plan_path = destination / PLAN_FILENAME
    cache_path = destination / CACHE_FILENAME
    series_cache_path = destination / SERIES_CACHE_FILENAME
    provenance_path = destination / PROVENANCE_FILENAME

    manifest = _read_manifest(source_path)
    geo_rows = manifest.loc[
        manifest["repository_name"].str.casefold().eq("geo")
        & manifest["repository_accession"].ne(""),
        REQUIRED_INPUT_COLUMNS,
    ].copy()
    geo_rows["repository_accession"] = geo_rows["repository_accession"].str.upper()

    cache = _load_cache(cache_path)
    resolutions: dict[str, dict[str, Any]] = {}
    valid_accessions = sorted(
        accession
        for accession in geo_rows["repository_accession"].unique()
        if GSM_PATTERN.fullmatch(accession)
    )
    invalid_accessions = sorted(
        set(geo_rows["repository_accession"].unique()) - set(valid_accessions)
    )
    for accession in invalid_accessions:
        resolutions[accession] = {
            "status": "invalid_accession",
            "gse_accessions": [],
            "gpl_accessions": [],
        }

    pending = []
    for accession in valid_accessions:
        if not force and accession in cache:
            resolutions[accession] = cache[accession]
        else:
            pending.append(accession)

    query_records = []
    api_key = os.environ.get("NCBI_API_KEY")
    for offset in range(0, len(pending), batch_size):
        batch = pending[offset : offset + batch_size]
        try:
            batch_records, query_record = _resolve_batch(batch, timeout, api_key)
            resolutions.update(batch_records)
            cache.update(batch_records)
            query_records.append(query_record)
            _write_cache(cache_path, cache)
        except Exception as exc:
            logger.error("GEO metadata batch failed for %s accessions: %s", len(batch), exc)
            query_records.append(
                {
                    "query_type": "gsm_resolution",
                    "started_at_utc": _utc_now(),
                    "search_endpoint": ESEARCH_URL,
                    "summary_endpoint": ESUMMARY_URL,
                    "requested_accessions": batch,
                    "status": "failed",
                    "error": str(exc),
                }
            )
            for accession in batch:
                resolutions[accession] = {
                    "status": "query_failed",
                    "gse_accessions": [],
                    "gpl_accessions": [],
                }

    if not cache_path.exists():
        _write_cache(cache_path, cache)

    gse_accessions = sorted(
        {
            gse_accession
            for resolution in resolutions.values()
            for gse_accession in resolution.get("gse_accessions", [])
        }
    )
    series_cache = _load_cache(series_cache_path)
    series_metadata: dict[str, dict[str, Any]] = {}
    pending_series = []
    for accession in gse_accessions:
        if not force and accession in series_cache:
            series_metadata[accession] = series_cache[accession]
        else:
            pending_series.append(accession)

    for offset in range(0, len(pending_series), batch_size):
        batch = pending_series[offset : offset + batch_size]
        try:
            batch_records, query_record = _resolve_series_batch(batch, timeout, api_key)
            series_metadata.update(batch_records)
            series_cache.update(batch_records)
            query_records.append(query_record)
            _write_cache(series_cache_path, series_cache)
        except Exception as exc:
            logger.error("GEO series metadata batch failed for %s accessions: %s", len(batch), exc)
            query_records.append(
                {
                    "query_type": "gse_metadata",
                    "started_at_utc": _utc_now(),
                    "search_endpoint": ESEARCH_URL,
                    "summary_endpoint": ESUMMARY_URL,
                    "requested_accessions": batch,
                    "status": "failed",
                    "error": str(exc),
                }
            )
            for accession in batch:
                series_metadata[accession] = {
                    "status": "query_failed",
                    "series_title": "",
                    "series_type": "unknown",
                    "series_sample_count": None,
                }

    if not series_cache_path.exists():
        _write_cache(series_cache_path, series_cache)

    plan = _build_plan(geo_rows, resolutions, series_metadata)
    plan.to_csv(plan_path, sep="\t", index=False, lineterminator="\n")
    statuses = {accession: resolutions[accession]["status"] for accession in resolutions}
    resolved = sum(status == "resolved" for status in statuses.values())
    unresolved = len(statuses) - resolved
    series_statuses = {
        accession: series_metadata[accession]["status"] for accession in series_metadata
    }
    unresolved_series = sum(value != "resolved" for value in series_statuses.values())
    status = (
        "success"
        if unresolved == 0 and unresolved_series == 0
        else ("failed" if resolved == 0 else "partial")
    )
    recommended_gses = set(
        plan.loc[plan["download_recommendation"].eq("download"), "gse_accession"]
    )
    skipped_superseries = set(
        plan.loc[plan["download_recommendation"].eq("skip_superseries"), "gse_accession"]
    )

    provenance = {
        "status": status,
        "run_started_at_utc": started_at,
        "planner": "geo_plan_module",
        "planner_version": PLANNER_VERSION,
        "python_version": platform.python_version(),
        "pandas_version": pd.__version__,
        "input": {
            "path": str(source_path),
            "size_bytes": source_path.stat().st_size,
            "sha256": _sha256(source_path),
        },
        "parameters": {
            "batch_size": batch_size,
            "timeout": timeout,
            "force": force,
        },
        "queries": query_records,
        "cache": {
            "path": str(cache_path),
            "size_bytes": cache_path.stat().st_size,
            "sha256": _sha256(cache_path),
            "records": len(cache),
        },
        "series_cache": {
            "path": str(series_cache_path),
            "size_bytes": series_cache_path.stat().st_size,
            "sha256": _sha256(series_cache_path),
            "records": len(series_cache),
        },
        "counts": {
            "geo_linked_rows": len(geo_rows),
            "unique_gsm_accessions": len(statuses),
            "resolved_gsm_accessions": resolved,
            "unresolved_gsm_accessions": unresolved,
            "unique_gse_accessions": len(gse_accessions),
            "resolved_gse_accessions": len(series_statuses) - unresolved_series,
            "unresolved_gse_accessions": unresolved_series,
            "recommended_gse_accessions": len(recommended_gses),
            "skipped_superseries_accessions": len(skipped_superseries),
            "download_plan_rows": len(plan),
        },
        "resolution_status_counts": {
            key: int(value)
            for key, value in pd.Series(list(statuses.values()), dtype=str).value_counts().items()
        },
        "series_status_counts": {
            key: int(value)
            for key, value in pd.Series(
                list(series_statuses.values()), dtype=str
            ).value_counts().items()
        },
        "output": {
            "path": str(plan_path),
            "size_bytes": plan_path.stat().st_size,
            "sha256": _sha256(plan_path),
            "rows": len(plan),
            "columns": list(plan.columns),
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
        description="Resolve ImmPort GSM links into a metadata-only GEO download plan."
    )
    parser.add_argument("manifest", help="Combined ImmPort sample_manifest.tsv")
    parser.add_argument("--output-dir", required=True, help="Directory for plan and metadata cache")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--force", action="store_true", help="Refresh cached GSM metadata")
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()
    try:
        provenance = plan_geo_downloads(
            args.manifest,
            args.output_dir,
            batch_size=args.batch_size,
            timeout=args.timeout,
            force=args.force,
        )
    except (GeoPlanError, OSError, ValueError) as exc:
        logger.error("GEO planning failed: %s", exc)
        return 1

    counts = provenance["counts"]
    logger.info(
        "Resolved %s/%s GSM accessions into %s plan rows",
        counts["resolved_gsm_accessions"],
        counts["unique_gsm_accessions"],
        counts["download_plan_rows"],
    )
    return 0 if provenance["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
