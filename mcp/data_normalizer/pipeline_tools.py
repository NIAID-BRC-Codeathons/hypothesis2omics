"""Controlled MCP adapters for the deterministic modules under ``data``.

The adapters call the existing Python functions directly. All input and output
paths are derived from one configured data root; MCP callers cannot provide
arbitrary filesystem paths. Fetch operations reuse valid cached files and load
credentials only through the data modules' existing environment-variable logic.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from data.geo_fetch_module import fetch_geo_datasets
from data.geo_matrix_parse_module import parse_geo_matrices as run_geo_matrix_parser
from data.geo_plan_module import plan_geo_downloads
from data.immport_batch_parse import run_batch as run_immport_batch_parser
from data.immport_fetch_module import fetch_immport_datasets
from data.list_directory import build_file_inventory

DEFAULT_DATA_ROOT = Path(__file__).resolve().parents[2] / "data"
PIPELINE_TOOLS_VERSION = "0.2.0"
SDY_PATTERN = re.compile(r"^SDY\d+$", re.IGNORECASE)
GSE_PATTERN = re.compile(r"^GSE\d+$", re.IGNORECASE)
PLAN_FILENAME = "geo_download_plan.tsv"
IMMPORT_MANIFEST_FILENAME = "sample_manifest.tsv"
GEO_SERIES_LINKS_FILENAME = "geo_series_links.tsv"
GEO_FETCH_SELECTION_FILENAME = "geo_fetch_selection.provenance.json"
GEO_PARSE_MANIFEST_FILENAME = "geo_matrix_parse_manifest.json"
MAX_INVENTORY_CANDIDATES = 100


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


def _sha256(path: Path, chunk_size: int = 65536) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path, required: set[str], label: str) -> pd.DataFrame:
    if not path.is_file():
        raise PipelineToolError(f"{label} does not exist: {path}")
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError) as exc:
        raise PipelineToolError(f"Could not read {label}: {exc}") from exc
    frame.columns = [str(column).strip().lower() for column in frame.columns]
    if len(frame.columns) != len(set(frame.columns)):
        raise PipelineToolError(f"{label} has duplicate normalized columns")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise PipelineToolError(f"{label} is missing columns: {missing}")
    return frame


class PipelineTools:
    """Execute the repository's ingestion stages within one fixed data root."""

    def __init__(self, data_root: str | Path = DEFAULT_DATA_ROOT) -> None:
        self.data_root = Path(data_root).expanduser().resolve()
        self.immport_cache = self.data_root / "immport_cache"
        self.immport_parsed = self.immport_cache / "parsed"
        self.immport_manifest = self.immport_parsed / IMMPORT_MANIFEST_FILENAME
        self.geo_series_links = self.immport_parsed / GEO_SERIES_LINKS_FILENAME
        self.geo_cache = self.data_root / "geo_cache"
        self.geo_plan = self.geo_cache / "plan"
        self.geo_plan_path = self.geo_plan / PLAN_FILENAME
        self.geo_parsed = self.geo_cache / "parsed"
        self.file_inventory = self.data_root / "file_inventory"

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

    def inventory_downloaded_files(self) -> dict[str, Any]:
        """Inventory cached downloads and expose configured parser candidates."""
        provenance = build_file_inventory(
            data_root=self.data_root,
            output_dir=self.file_inventory,
        )
        candidate_output = next(
            item for item in provenance["outputs"] if item["type"] == "parser_candidates"
        )
        try:
            candidates = pd.read_csv(
                candidate_output["path"],
                sep="\t",
                dtype=str,
                keep_default_na=False,
            )
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            raise PipelineToolError(f"Could not read parser candidate inventory: {exc}") from exc
        supported = candidates.loc[candidates["parser_support"].eq("supported")]
        exposed = supported.head(MAX_INVENTORY_CANDIDATES).to_dict(orient="records")
        return {
            "operation": "inventory_downloaded_files",
            "status": provenance["status"],
            "counts": provenance["counts"],
            "inventory_path": provenance["outputs"][0]["path"],
            "candidate_path": candidate_output["path"],
            "provenance_path": provenance["provenance_path"],
            "supported_candidates": exposed,
            "supported_candidates_truncated": len(supported) > len(exposed),
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

    def fetch_linked_geo(self) -> dict[str, Any]:
        """Fetch GSE accessions linked to GEO samples in parsed ImmPort outputs."""
        manifest = _read_tsv(
            self.immport_manifest,
            {"study_accession", "repository_name", "repository_accession"},
            "ImmPort sample manifest",
        )
        links = _read_tsv(
            self.geo_series_links,
            {"study_accession", "gse_accession"},
            "ImmPort GEO series links",
        )

        geo_rows = manifest.loc[
            manifest["repository_name"].str.strip().str.casefold().eq("geo")
            & manifest["repository_accession"].str.strip().ne("")
        ].copy()
        if geo_rows.empty:
            raise PipelineToolError("ImmPort sample manifest has no GEO-linked samples")
        geo_studies = _normalized_accessions(
            list(dict.fromkeys(geo_rows["study_accession"])),
            SDY_PATTERN,
            "GEO-linked study accessions",
        )

        links["study_accession"] = links["study_accession"].str.strip().str.upper()
        links["gse_accession"] = links["gse_accession"].str.strip().str.upper()
        invalid_studies = sorted(
            value
            for value in set(links["study_accession"])
            if not SDY_PATTERN.fullmatch(value)
        )
        invalid_gses = sorted(
            value for value in set(links["gse_accession"]) if not GSE_PATTERN.fullmatch(value)
        )
        if invalid_studies:
            raise PipelineToolError(
                f"ImmPort GEO series links contain invalid studies: {invalid_studies}"
            )
        if invalid_gses:
            raise PipelineToolError(
                f"ImmPort GEO series links contain invalid GSE accessions: {invalid_gses}"
            )

        linked_studies = set(links["study_accession"])
        missing_studies = sorted(set(geo_studies) - linked_studies)
        if missing_studies:
            raise PipelineToolError(
                "GEO-linked ImmPort studies have no GSE link: "
                + ", ".join(missing_studies)
            )
        selected = links.loc[
            links["study_accession"].isin(geo_studies),
            "gse_accession",
        ]
        accessions = list(dict.fromkeys(selected))

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

        self.geo_cache.mkdir(parents=True, exist_ok=True)
        provenance_path = self.geo_cache / GEO_FETCH_SELECTION_FILENAME
        provenance = {
            "status": _overall_status([item["status"] for item in results]),
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "tool": "mcp.data_normalizer.fetch_linked_geo",
            "tool_version": PIPELINE_TOOLS_VERSION,
            "inputs": {
                "sample_manifest": {
                    "path": str(self.immport_manifest),
                    "sha256": _sha256(self.immport_manifest),
                },
                "geo_series_links": {
                    "path": str(self.geo_series_links),
                    "sha256": _sha256(self.geo_series_links),
                },
            },
            "selected_gse_accessions": accessions,
            "outputs": {
                "cache_directory": str(self.geo_cache),
                "provenance_log": str(self.geo_cache / "provenance_log.jsonl"),
                "selection_provenance": str(provenance_path),
            },
            "results": results,
        }
        provenance_path.write_text(
            json.dumps(provenance, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return {
            "operation": "fetch_linked_geo",
            "status": provenance["status"],
            "linked_accessions": accessions,
            "cache_root": str(self.geo_cache),
            "provenance_path": str(provenance_path),
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

    def parse_linked_geo_matrices(self) -> dict[str, Any]:
        """Derive analysis units from parsed ImmPort links and parse cached matrices."""
        provenance = run_geo_matrix_parser(
            None,
            self.immport_manifest,
            self.geo_cache,
            self.geo_parsed,
            geo_series_links_path=self.geo_series_links,
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
            "operation": "parse_linked_geo_matrices",
            "status": provenance["status"],
            "selection_mode": provenance["selection_mode"],
            "selection_path": provenance["selection_path"],
            "counts": provenance["counts"],
            "units": units,
            "manifest_path": str(self.geo_parsed / GEO_PARSE_MANIFEST_FILENAME),
            "duration_sec": provenance["duration_sec"],
        }
