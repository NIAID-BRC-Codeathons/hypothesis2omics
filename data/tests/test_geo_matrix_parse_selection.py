"""Tests for deriving GEO matrix parse units from parsed ImmPort links."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from data.geo_matrix_parse_module import (
    GeoMatrixParseError,
    MatrixHeader,
    _derive_download_plan,
)


def _header(gse: str, gsms: list[str], gpl: str) -> MatrixHeader:
    metadata = pd.DataFrame(
        [
            {
                "gsm_accession": gsm,
                "attribute": "platform_id",
                "occurrence": 0,
                "attribute_order": 0,
                "value": gpl,
            }
            for gsm in gsms
        ]
    )
    return MatrixHeader(gse, gsms, {gpl}, metadata)


class GeoMatrixParseSelectionTests(unittest.TestCase):
    def test_partitions_samples_across_linked_series(self) -> None:
        manifest = pd.DataFrame(
            [
                {
                    "study_accession": "SDY1",
                    "experiment_accession": "EXP1",
                    "repository_name": "GEO",
                    "repository_accession": gsm,
                }
                for gsm in ["GSM1", "GSM2", "GSM3"]
            ]
        )
        links = pd.DataFrame(
            [
                {"study_accession": "SDY1", "gse_accession": "GSE1"},
                {"study_accession": "SDY1", "gse_accession": "GSE2"},
            ]
        )
        headers = {
            "GSE1_series_matrix.txt.gz": _header("GSE1", ["GSM1"], "GPL1"),
            "GSE2_series_matrix.txt.gz": _header("GSE2", ["GSM2", "GSM3"], "GPL2"),
        }

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            for gse, filename in [
                ("GSE1", "GSE1_series_matrix.txt.gz"),
                ("GSE2", "GSE2_series_matrix.txt.gz"),
            ]:
                path = cache / gse / filename
                path.parent.mkdir()
                path.touch()
            with patch(
                "data.geo_matrix_parse_module._read_matrix_header",
                side_effect=lambda path: headers[path.name],
            ):
                plan = _derive_download_plan(links, manifest, cache)

        self.assertEqual(list(plan["gse_accession"]), ["GSE1", "GSE2"])
        self.assertEqual(list(plan["requested_sample_count"]), ["1", "2"])
        self.assertEqual(list(plan["gsm_accessions"]), ["GSM1", "GSM2;GSM3"])

    def test_rejects_sample_present_in_multiple_linked_matrices(self) -> None:
        manifest = pd.DataFrame(
            [
                {
                    "study_accession": "SDY1",
                    "experiment_accession": "EXP1",
                    "repository_name": "GEO",
                    "repository_accession": "GSM1",
                }
            ]
        )
        links = pd.DataFrame(
            [
                {"study_accession": "SDY1", "gse_accession": "GSE1"},
                {"study_accession": "SDY1", "gse_accession": "GSE2"},
            ]
        )

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            headers = {}
            for gse in ["GSE1", "GSE2"]:
                path = cache / gse / f"{gse}_series_matrix.txt.gz"
                path.parent.mkdir()
                path.touch()
                headers[path.name] = _header(gse, ["GSM1"], "GPL1")
            with patch(
                "data.geo_matrix_parse_module._read_matrix_header",
                side_effect=lambda path: headers[path.name],
            ):
                with self.assertRaisesRegex(GeoMatrixParseError, r"ambiguous=\['GSM1'\]"):
                    _derive_download_plan(links, manifest, cache)


if __name__ == "__main__":
    unittest.main()
