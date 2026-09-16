"""Read-only accessors for parsed GEO analysis units.

Inputs are the run manifest and unit files written by
``data/geo_matrix_parse_module.py``. Outputs are JSON-serializable dictionaries
intended for direct calls or MCP tools. No network access or data mutation occurs.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

RUN_MANIFEST_FILENAME = "geo_matrix_parse_manifest.json"
EXPRESSION_FILENAME = "expression.tsv.gz"
SAMPLES_FILENAME = "samples.tsv"
GEO_METADATA_FILENAME = "geo_sample_metadata_long.tsv"
DEFAULT_PARSED_ROOT = (
    Path(__file__).resolve().parents[1] / "data" / "geo_cache" / "parsed"
)
MAX_SAMPLE_PAGE_SIZE = 100
GSM_PATTERN = re.compile(r"^GSM\d+$", re.IGNORECASE)


class GeoMcpDataError(ValueError):
    """Raised when parsed GEO data cannot satisfy a read-only tool request."""


def analysis_unit_id(unit: dict[str, Any]) -> str:
    """Return the stable identifier for a parsed GEO analysis unit."""
    return "/".join(
        [
            str(unit["study_accession"]),
            str(unit["experiment_accession"]),
            f"{unit['gse_accession']}_{unit['gpl_accession']}",
        ]
    )


class GeoDataStore:
    """Load and query one GEO parser output root without modifying it."""

    def __init__(self, parsed_root: str | Path = DEFAULT_PARSED_ROOT) -> None:
        self.parsed_root = Path(parsed_root).expanduser().resolve()
        self.manifest_path = self.parsed_root / RUN_MANIFEST_FILENAME

    def _load_manifest(self) -> dict[str, Any]:
        if not self.manifest_path.is_file():
            raise GeoMcpDataError(
                f"GEO parse manifest does not exist: {self.manifest_path}"
            )
        try:
            payload = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise GeoMcpDataError(
                f"Could not read GEO parse manifest {self.manifest_path}: {exc}"
            ) from exc
        if not isinstance(payload, dict) or not isinstance(payload.get("units"), list):
            raise GeoMcpDataError(
                f"GEO parse manifest has no valid units list: {self.manifest_path}"
            )
        return payload

    def _unit_index(self) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
        manifest = self._load_manifest()
        index: dict[str, dict[str, Any]] = {}
        for unit in manifest["units"]:
            try:
                unit_id = analysis_unit_id(unit)
            except KeyError as exc:
                raise GeoMcpDataError(
                    f"GEO parse manifest unit is missing {exc.args[0]!r}"
                ) from exc
            if unit_id in index:
                raise GeoMcpDataError(f"Duplicate analysis unit ID: {unit_id}")
            index[unit_id] = unit
        return manifest, index

    def _resolve_unit(self, unit_id: str) -> tuple[dict[str, Any], Path]:
        _, index = self._unit_index()
        if unit_id not in index:
            available = ", ".join(sorted(index))
            raise GeoMcpDataError(
                f"Unknown analysis unit {unit_id!r}; available units: {available}"
            )
        unit = index[unit_id]
        unit_dir = (
            self.parsed_root
            / str(unit["study_accession"])
            / str(unit["experiment_accession"])
            / f"{unit['gse_accession']}_{unit['gpl_accession']}"
        ).resolve()
        if not unit_dir.is_relative_to(self.parsed_root):
            raise GeoMcpDataError(f"Analysis unit resolves outside parsed root: {unit_id}")
        if not unit_dir.is_dir():
            raise GeoMcpDataError(f"Analysis unit directory does not exist: {unit_dir}")
        return unit, unit_dir

    def _unit_file(self, unit_dir: Path, filename: str) -> Path:
        path = (unit_dir / filename).resolve()
        if not path.is_relative_to(self.parsed_root):
            raise GeoMcpDataError(f"Parsed GEO output resolves outside parsed root: {path}")
        if not path.is_file():
            raise GeoMcpDataError(f"Parsed GEO output does not exist: {path}")
        return path

    @staticmethod
    def _read_tsv(path: Path) -> pd.DataFrame:
        if not path.is_file():
            raise GeoMcpDataError(f"Parsed GEO output does not exist: {path}")
        try:
            return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            raise GeoMcpDataError(f"Could not read parsed GEO output {path}: {exc}") from exc

    @staticmethod
    def _output_record(unit: dict[str, Any], output_type: str, path: Path) -> dict[str, Any]:
        records = [item for item in unit.get("outputs", []) if item.get("type") == output_type]
        if len(records) != 1:
            raise GeoMcpDataError(
                f"Analysis unit must have one {output_type!r} provenance record"
            )
        if not path.is_file():
            raise GeoMcpDataError(f"Parsed GEO output does not exist: {path}")
        record = records[0]
        return {
            "type": output_type,
            "path": str(path),
            "size_bytes": record.get("size_bytes"),
            "sha256": record.get("sha256"),
        }

    def list_analysis_units(self) -> dict[str, Any]:
        """List parsed GEO units and their high-level counts."""
        manifest, index = self._unit_index()
        units = []
        for unit_id, unit in index.items():
            units.append(
                {
                    "analysis_unit_id": unit_id,
                    "study_accession": unit["study_accession"],
                    "experiment_accession": unit["experiment_accession"],
                    "gse_accession": unit["gse_accession"],
                    "gpl_accession": unit["gpl_accession"],
                    "status": unit.get("status"),
                    "counts": unit.get("counts", {}),
                }
            )
        return {
            "run_status": manifest.get("status"),
            "run_started_at_utc": manifest.get("run_started_at_utc"),
            "parser_version": manifest.get("parser_version"),
            "counts": manifest.get("counts", {}),
            "units": sorted(units, key=lambda item: item["analysis_unit_id"]),
        }

    def get_analysis_unit(self, unit_id: str) -> dict[str, Any]:
        """Return provenance and available metadata fields for one GEO unit."""
        unit, unit_dir = self._resolve_unit(unit_id)
        expression_path = self._unit_file(unit_dir, EXPRESSION_FILENAME)
        samples_path = self._unit_file(unit_dir, SAMPLES_FILENAME)
        metadata_path = self._unit_file(unit_dir, GEO_METADATA_FILENAME)
        samples = self._read_tsv(samples_path)
        metadata = self._read_tsv(metadata_path)
        required_metadata = {"attribute", "occurrence"}
        if not required_metadata.issubset(metadata.columns):
            missing = sorted(required_metadata - set(metadata.columns))
            raise GeoMcpDataError(
                f"GEO metadata is missing columns: {missing}"
            )
        occurrences = pd.to_numeric(metadata["occurrence"], errors="coerce")
        attributes = []
        for attribute in metadata["attribute"].drop_duplicates():
            mask = metadata["attribute"].eq(attribute)
            maximum = occurrences.loc[mask].max()
            attributes.append(
                {
                    "attribute": attribute,
                    "occurrences_per_sample": int(maximum) + 1 if pd.notna(maximum) else None,
                }
            )
        output_paths = {
            "expression": expression_path,
            "samples": samples_path,
            "geo_sample_metadata": metadata_path,
        }
        return {
            "analysis_unit_id": unit_id,
            "study_accession": unit["study_accession"],
            "experiment_accession": unit["experiment_accession"],
            "gse_accession": unit["gse_accession"],
            "gpl_accession": unit["gpl_accession"],
            "status": unit.get("status"),
            "run_started_at_utc": unit.get("run_started_at_utc"),
            "duration_sec": unit.get("duration_sec"),
            "counts": unit.get("counts", {}),
            "source_matrix": unit.get("source_matrix", {}),
            "sample_fields": list(samples.columns),
            "geo_metadata_attributes": attributes,
            "outputs": [
                self._output_record(unit, output_type, path)
                for output_type, path in output_paths.items()
            ],
        }

    def get_sample_metadata(
        self,
        unit_id: str,
        gsm_accessions: list[str] | None = None,
        attributes: list[str] | None = None,
        limit: int = 25,
        offset: int = 0,
    ) -> dict[str, Any]:
        """Return a bounded page of linked ImmPort and GEO sample metadata."""
        if limit < 1 or limit > MAX_SAMPLE_PAGE_SIZE:
            raise GeoMcpDataError(
                f"limit must be between 1 and {MAX_SAMPLE_PAGE_SIZE}, received {limit}"
            )
        if offset < 0:
            raise GeoMcpDataError(f"offset must be non-negative, received {offset}")

        _, unit_dir = self._resolve_unit(unit_id)
        samples = self._read_tsv(self._unit_file(unit_dir, SAMPLES_FILENAME))
        metadata = self._read_tsv(self._unit_file(unit_dir, GEO_METADATA_FILENAME))
        if "repository_accession" not in samples.columns:
            raise GeoMcpDataError("Linked sample table has no repository_accession column")
        required_metadata = {"gsm_accession", "attribute"}
        if not required_metadata.issubset(metadata.columns):
            missing = sorted(required_metadata - set(metadata.columns))
            raise GeoMcpDataError(
                f"GEO metadata is missing columns: {missing}"
            )

        samples["repository_accession"] = samples["repository_accession"].str.upper()
        available_gsms = list(samples["repository_accession"])
        if gsm_accessions is not None:
            requested = [str(value).strip().upper() for value in gsm_accessions]
            if not requested or any(not GSM_PATTERN.fullmatch(value) for value in requested):
                raise GeoMcpDataError("gsm_accessions must contain valid GSM accessions")
            if len(requested) != len(set(requested)):
                raise GeoMcpDataError("gsm_accessions contains duplicates")
            missing = sorted(set(requested) - set(available_gsms))
            if missing:
                raise GeoMcpDataError(f"GSM accessions are not in {unit_id}: {missing}")
            order = {gsm: position for position, gsm in enumerate(requested)}
            samples = samples.loc[samples["repository_accession"].isin(requested)].copy()
            samples["_order"] = samples["repository_accession"].map(order)
            samples = samples.sort_values("_order").drop(columns="_order")

        available_attributes = list(metadata["attribute"].drop_duplicates())
        if attributes is not None:
            requested_attributes = [str(value).strip() for value in attributes]
            if not requested_attributes or any(not value for value in requested_attributes):
                raise GeoMcpDataError("attributes must contain non-empty names")
            if len(requested_attributes) != len(set(requested_attributes)):
                raise GeoMcpDataError("attributes contains duplicates")
            missing_attributes = sorted(set(requested_attributes) - set(available_attributes))
            if missing_attributes:
                raise GeoMcpDataError(
                    f"Metadata attributes are not in {unit_id}: {missing_attributes}"
                )
            metadata = metadata.loc[metadata["attribute"].isin(requested_attributes)]

        total_matching = len(samples)
        page = samples.iloc[offset : offset + limit].copy()
        page_gsms = list(page["repository_accession"])
        metadata_page = metadata.loc[metadata["gsm_accession"].isin(page_gsms)].copy()
        sample_order = {gsm: position for position, gsm in enumerate(page_gsms)}
        metadata_page["_sample_order"] = metadata_page["gsm_accession"].map(sample_order)
        sort_columns = ["_sample_order"]
        sort_columns.extend(
            column for column in ["attribute_order", "occurrence"] if column in metadata_page
        )
        metadata_page = metadata_page.sort_values(sort_columns).drop(columns="_sample_order")
        return {
            "analysis_unit_id": unit_id,
            "total_unit_samples": len(available_gsms),
            "total_matching_samples": total_matching,
            "offset": offset,
            "limit": limit,
            "returned_samples": len(page),
            "next_offset": offset + len(page) if offset + len(page) < total_matching else None,
            "available_attributes": available_attributes,
            "linked_samples": page.to_dict(orient="records"),
            "geo_metadata": metadata_page.to_dict(orient="records"),
        }

    def get_expression_info(self, unit_id: str) -> dict[str, Any]:
        """Return matrix dimensions and sample columns without expression values."""
        unit, unit_dir = self._resolve_unit(unit_id)
        path = self._unit_file(unit_dir, EXPRESSION_FILENAME)
        output = self._output_record(unit, "expression", path)
        try:
            with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
                header = next(csv.reader([handle.readline()], delimiter="\t"))
        except (OSError, UnicodeError, csv.Error) as exc:
            raise GeoMcpDataError(f"Could not read expression header {path}: {exc}") from exc
        if not header or header[0] != "ID_REF":
            raise GeoMcpDataError(f"Expression matrix does not begin with ID_REF: {path}")
        counts = unit.get("counts", {})
        return {
            "analysis_unit_id": unit_id,
            "id_column": header[0],
            "sample_accessions": header[1:],
            "sample_count": counts.get("samples"),
            "probe_count": counts.get("probes"),
            "path": output["path"],
            "size_bytes": output["size_bytes"],
            "sha256": output["sha256"],
        }
