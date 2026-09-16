"""Controlled MCP adapters for the deterministic modules under ``data``.

The adapters call the existing Python functions directly. All input and output
paths are derived from one configured data root; MCP callers cannot provide
arbitrary filesystem paths. Fetch operations reuse valid cached files and load
credentials only through the data modules' existing environment-variable logic.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from data.geo_fetch_module import fetch_geo_datasets
from data.geo_matrix_parse_module import parse_geo_matrices as run_geo_matrix_parser
from data.geo_plan_module import plan_geo_downloads
from data.immport_batch_parse import run_batch as run_immport_batch_parser
from data.immport_fetch_module import fetch_immport_datasets

DEFAULT_DATA_ROOT = Path(__file__).resolve().parents[1] / "data"
SDY_PATTERN = re.compile(r"^SDY\d+$", re.IGNORECASE)
GSE_PATTERN = re.compile(r"^GSE\d+$", re.IGNORECASE)
PLAN_FILENAME = "geo_download_plan.tsv"
IMMPORT_MANIFEST_FILENAME = "sample_manifest.tsv"
GEO_PARSE_MANIFEST_FILENAME = "geo_matrix_parse_manifest.json"


class PipelineToolError(RuntimeError):
    """Raised when an executable pipeline tool cannot satisfy its contract."""


def _normalized_accessions(
    values: Sequence[str],
    pattern: re.Pattern[str],
    label: str,
) -> list[str]:
    normalized = [str(value).strip().upper() for value in values]
    if not normalized or any(not pattern.fullmatch(value) for value in normalized):
        raise PipelineToolError(f"{label} must contain at least one valid accession")
    if len(normalized) != len(set(normalized)):
        raise PipelineToolError(f"{label} contains duplicate accessions")
    return normalized


def _overall_status(statuses: list[str]) -> str:
    completed = {"success", "cached"}
    if statuses and all(status in completed for status in statuses):
        return "success"
    if statuses and all(status == "failed" for status in statuses):
        return "failed"
    return "partial"


def _artifact_summary(artifact: Any) -> dict[str, Any]:
    return {
        "artifact_type": artifact.artifact_type,
        "status": artifact.status,
        "local_file": artifact.local_file,
        "size_bytes": artifact.size_bytes,
        "sha256": artifact.sha256,
        "error": artifact.error,
    }


class PipelineTools:
    """Execute the repository's ingestion stages within one fixed data root."""

    def __init__(self, data_root: str | Path = DEFAULT_DATA_ROOT) -> None:
        self.data_root = Path(data_root).expanduser().resolve()
        self.immport_cache = self.data_root / "immport_cache"
        self.immport_parsed = self.immport_cache / "parsed"
        self.immport_manifest = self.immport_parsed / IMMPORT_MANIFEST_FILENAME
        self.geo_cache = self.data_root / "geo_cache"
        self.geo_plan = self.geo_cache / "plan"
        self.geo_plan_path = self.geo_plan / PLAN_FILENAME
        self.geo_parsed = self.geo_cache / "parsed"

    def fetch_immport_studies(
        self,
        study_accessions: list[str],
        max_files_per_study: int | None = None,
    ) -> dict[str, Any]:
        """Fetch ImmPort studies into the configured cache with provenance."""
        accessions = _normalized_accessions(
            study_accessions,
            SDY_PATTERN,
            "study_accessions",
        )
        if max_files_per_study is not None and max_files_per_study < 1:
            raise PipelineToolError("max_files_per_study must be positive when provided")
        records = fetch_immport_datasets(
            accessions,
            destdir=str(self.immport_cache),
            force=False,
            max_files_per_study=max_files_per_study,
            provenance_log_path=str(self.immport_cache / "provenance_log.jsonl"),
        )
        results = []
        for record in records:
            results.append(
                {
                    "study_accession": record.sdy_id,
                    "status": record.status,
                    "files_listed": record.n_files_listed,
                    "files_downloaded": record.n_files_downloaded,
                    "duration_sec": record.duration_sec,
                    "error": record.error,
                    "artifacts": [_artifact_summary(item) for item in record.artifacts],
                }
            )
        return {
            "operation": "fetch_immport_studies",
            "status": _overall_status([item["status"] for item in results]),
            "cache_root": str(self.immport_cache),
            "results": results,
        }

    def parse_immport_studies(self) -> dict[str, Any]:
        """Parse every cached ImmPort study and combine its sample linkage."""
        provenance = run_immport_batch_parser(self.immport_cache)
        return {
            "operation": "parse_immport_studies",
            "status": provenance["status"],
            "counts": provenance["counts"],
            "studies": provenance["studies"],
            "combined_output": provenance["combined_output"],
            "manifest_path": str(self.immport_parsed / "batch_manifest.json"),
            "duration_sec": provenance["duration_sec"],
        }

    def plan_geo_retrieval(
        self,
        batch_size: int = 100,
        timeout: int = 60,
    ) -> dict[str, Any]:
        """Resolve linked GSM accessions and write the bounded GEO download plan."""
        provenance = plan_geo_downloads(
            self.immport_manifest,
            self.geo_plan,
            batch_size=batch_size,
            timeout=timeout,
            force=False,
        )
        return {
            "operation": "plan_geo_retrieval",
            "status": provenance["status"],
            "counts": provenance["counts"],
            "resolution_status_counts": provenance["resolution_status_counts"],
            "series_status_counts": provenance["series_status_counts"],
            "plan_path": str(self.geo_plan_path),
            "provenance_path": str(self.geo_plan / "geo_download_plan.provenance.json"),
            "duration_sec": provenance["duration_sec"],
        }

    def fetch_planned_geo(self) -> dict[str, Any]:
        """Fetch only GSE accessions marked ``download`` in the current plan."""
        if not self.geo_plan_path.is_file():
            raise PipelineToolError(f"GEO download plan does not exist: {self.geo_plan_path}")
        try:
            plan = pd.read_csv(
                self.geo_plan_path,
                sep="\t",
                dtype=str,
                keep_default_na=False,
            )
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            raise PipelineToolError(f"Could not read GEO download plan: {exc}") from exc
        required = {"gse_accession", "download_recommendation"}
        missing = sorted(required - set(plan.columns))
        if missing:
            raise PipelineToolError(f"GEO download plan is missing columns: {missing}")
        selected = plan.loc[
            plan["download_recommendation"].eq("download"),
            "gse_accession",
        ]
        accessions = _normalized_accessions(
            list(dict.fromkeys(selected)),
            GSE_PATTERN,
            "planned_gse_accessions",
        )
        records = fetch_geo_datasets(
            accessions,
            destdir=str(self.geo_cache),
            force=False,
            provenance_log_path=str(self.geo_cache / "provenance_log.jsonl"),
        )
        results = []
        for record in records:
            results.append(
                {
                    "gse_accession": record.gse_id,
                    "status": record.status,
                    "platforms": record.n_platforms,
                    "samples": record.n_samples,
                    "duration_sec": record.duration_sec,
                    "error": record.error,
                    "artifacts": [_artifact_summary(item) for item in record.artifacts],
                }
            )
        return {
            "operation": "fetch_planned_geo",
            "status": _overall_status([item["status"] for item in results]),
            "planned_accessions": accessions,
            "cache_root": str(self.geo_cache),
            "results": results,
        }

    def parse_geo_matrices(self) -> dict[str, Any]:
        """Parse all downloaded plan units into bounded expression artifacts."""
        provenance = run_geo_matrix_parser(
            self.geo_plan_path,
            self.immport_manifest,
            self.geo_cache,
            self.geo_parsed,
        )
        units = [
            {
                "analysis_unit_id": "/".join(
                    [
                        unit["study_accession"],
                        unit["experiment_accession"],
                        f"{unit['gse_accession']}_{unit['gpl_accession']}",
                    ]
                ),
                "status": unit["status"],
                "counts": unit.get("counts", {}),
                "error": unit.get("error"),
            }
            for unit in provenance["units"]
        ]
        return {
            "operation": "parse_geo_matrices",
            "status": provenance["status"],
            "counts": provenance["counts"],
            "units": units,
            "manifest_path": str(self.geo_parsed / GEO_PARSE_MANIFEST_FILENAME),
            "duration_sec": provenance["duration_sec"],
        }
