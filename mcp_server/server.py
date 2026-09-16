"""Expose the Hypothesis2Omics data pipeline and parsed GEO outputs through MCP.

The server reads ``HYPOTHESIS2OMICS_GEO_PARSED_ROOT`` when set and otherwise
uses ``data/geo_cache/parsed``. Executable tools use ``HYPOTHESIS2OMICS_DATA_ROOT``
when set and otherwise use ``data``. The server uses stdio by default.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

from fastmcp import FastMCP

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from mcp_server.geo_tools import GeoDataStore
from mcp_server.pipeline_tools import DEFAULT_DATA_ROOT, PipelineTools

PARSED_ROOT_ENV = "HYPOTHESIS2OMICS_GEO_PARSED_ROOT"
DATA_ROOT_ENV = "HYPOTHESIS2OMICS_DATA_ROOT"


def create_server(
    parsed_root: str | Path | None = None,
    data_root: str | Path | None = None,
) -> FastMCP:
    """Create a server bound to one parsed GEO output root."""
    configured_data_root = data_root or os.environ.get(DATA_ROOT_ENV, DEFAULT_DATA_ROOT)
    pipeline = PipelineTools(configured_data_root)
    configured_parsed_root = parsed_root or os.environ.get(
        PARSED_ROOT_ENV,
        pipeline.geo_parsed,
    )
    store = GeoDataStore(configured_parsed_root)
    server = FastMCP(
        name="hypothesis2omics",
        instructions=(
            "Execute provenance-backed ingestion stages and inspect parsed GEO units. "
            "Run stages separately and inspect each status before continuing."
        ),
    )

    @server.tool
    def list_analysis_units() -> dict[str, Any]:
        """List every parsed GEO analysis unit and its sample/probe counts."""
        return store.list_analysis_units()

    @server.tool
    def get_analysis_unit(unit_id: str) -> dict[str, Any]:
        """Get one unit's provenance, output hashes, and available metadata fields."""
        return store.get_analysis_unit(unit_id)

    @server.tool
    def get_sample_metadata(
        unit_id: str,
        gsm_accessions: list[str] | None = None,
        attributes: list[str] | None = None,
        limit: int = 25,
        offset: int = 0,
    ) -> dict[str, Any]:
        """Get a filtered, paginated page of linked ImmPort and GEO sample metadata."""
        return store.get_sample_metadata(
            unit_id,
            gsm_accessions=gsm_accessions,
            attributes=attributes,
            limit=limit,
            offset=offset,
        )

    @server.tool
    def get_expression_info(unit_id: str) -> dict[str, Any]:
        """Get expression dimensions, samples, path, and hash without matrix values."""
        return store.get_expression_info(unit_id)

    @server.tool
    def fetch_immport_studies(
        study_accessions: list[str],
        max_files_per_study: int | None = None,
    ) -> dict[str, Any]:
        """Fetch specified ImmPort studies using server-side credentials and caching."""
        return pipeline.fetch_immport_studies(
            study_accessions,
            max_files_per_study=max_files_per_study,
        )

    @server.tool
    def parse_immport_studies() -> dict[str, Any]:
        """Parse all cached ImmPort studies and combine their sample linkage."""
        return pipeline.parse_immport_studies()

    @server.tool
    def plan_geo_retrieval(batch_size: int = 100, timeout: int = 60) -> dict[str, Any]:
        """Resolve linked GSM accessions and write the cached GEO retrieval plan."""
        return pipeline.plan_geo_retrieval(batch_size=batch_size, timeout=timeout)

    @server.tool
    def fetch_planned_geo() -> dict[str, Any]:
        """Fetch only GSE accessions marked for download in the current GEO plan."""
        return pipeline.fetch_planned_geo()

    @server.tool
    def parse_geo_matrices() -> dict[str, Any]:
        """Parse downloaded GEO matrices into bounded analysis units."""
        return pipeline.parse_geo_matrices()

    return server


mcp = create_server()


if __name__ == "__main__":
    mcp.run(transport="stdio")
