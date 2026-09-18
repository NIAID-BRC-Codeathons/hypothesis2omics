"""Tests for fetching GEO series from parsed ImmPort study links."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
MCP_ROOT = REPOSITORY_ROOT / "mcp"
if str(MCP_ROOT) not in sys.path:
    sys.path.insert(0, str(MCP_ROOT))

from data_normalizer.pipeline_tools import PipelineToolError, PipelineTools  # noqa: E402


def _fetch_record(gse_accession: str) -> SimpleNamespace:
    return SimpleNamespace(
        gse_id=gse_accession,
        status="success",
        n_platforms=1,
        n_samples=10,
        duration_sec=0.1,
        error=None,
        artifacts=[],
    )


class PipelineLinkedGeoTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.pipeline = PipelineTools(Path(self.temporary_directory.name) / "data")
        self.pipeline.immport_parsed.mkdir(parents=True)
        self.pipeline.immport_manifest.write_text(
            "study_accession\trepository_name\trepository_accession\n"
            "SDY1\tGEO\tGSM1\n"
            "SDY2\tGEO\tGSM2\n",
            encoding="utf-8",
        )
        self.pipeline.geo_series_links.write_text(
            "study_accession\tgse_accession\n"
            "SDY1\tGSE10\n"
            "SDY2\tGSE20\n"
            "SDY2\tGSE20\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @patch("data_normalizer.pipeline_tools.fetch_geo_datasets")
    def test_fetches_deduplicated_linked_series_with_provenance(self, fetch_mock) -> None:
        fetch_mock.return_value = [_fetch_record("GSE10"), _fetch_record("GSE20")]

        result = self.pipeline.fetch_linked_geo()

        accessions = fetch_mock.call_args.args[0]
        self.assertEqual(accessions, ["GSE10", "GSE20"])
        self.assertEqual(result["status"], "success")
        provenance = json.loads(Path(result["provenance_path"]).read_text(encoding="utf-8"))
        self.assertEqual(provenance["selected_gse_accessions"], ["GSE10", "GSE20"])
        self.assertIn("sha256", provenance["inputs"]["geo_series_links"])
        self.assertIn("sha256", provenance["inputs"]["sample_manifest"])

    @patch("data_normalizer.pipeline_tools.fetch_geo_datasets")
    def test_rejects_geo_linked_study_without_gse(self, fetch_mock) -> None:
        self.pipeline.geo_series_links.write_text(
            "study_accession\tgse_accession\nSDY1\tGSE10\n",
            encoding="utf-8",
        )

        with self.assertRaisesRegex(PipelineToolError, "SDY2"):
            self.pipeline.fetch_linked_geo()

        fetch_mock.assert_not_called()

    @patch("data_normalizer.pipeline_tools.run_geo_matrix_parser")
    def test_parses_derived_linked_units(self, parser_mock) -> None:
        parser_mock.return_value = {
            "status": "success",
            "selection_mode": "derived_from_parsed_immport_links",
            "selection_path": "derived.tsv",
            "counts": {"download_units": 1},
            "units": [
                {
                    "study_accession": "SDY1",
                    "experiment_accession": "EXP1",
                    "gse_accession": "GSE10",
                    "gpl_accession": "GPL1",
                    "status": "success",
                    "counts": {"samples": 10},
                }
            ],
            "duration_sec": 0.1,
        }

        result = self.pipeline.parse_linked_geo_matrices()

        parser_mock.assert_called_once_with(
            None,
            self.pipeline.immport_manifest,
            self.pipeline.geo_cache,
            self.pipeline.geo_parsed,
            geo_series_links_path=self.pipeline.geo_series_links,
        )
        self.assertEqual(result["selection_mode"], "derived_from_parsed_immport_links")
        self.assertEqual(result["units"][0]["analysis_unit_id"], "SDY1/EXP1/GSE10_GPL1")


if __name__ == "__main__":
    unittest.main()
