"""Tests for selecting GEO downloads from parsed ImmPort study links."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from data.geo_fetch_module import (
    GeoSelectionError,
    load_gse_ids_from_parsed_links,
    main,
)


class GeoFetchSelectionTests(unittest.TestCase):
    def test_loads_normalizes_and_deduplicates_parsed_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            links_path = Path(directory) / "geo_series_links.tsv"
            links_path.write_text(
                "study_accession\tgse_accession\n"
                "SDY1\tgse13485\n"
                "SDY2\tGSE13699\n"
                "SDY3\tGSE13485\n",
                encoding="utf-8",
            )

            accessions = load_gse_ids_from_parsed_links(links_path)

        self.assertEqual(accessions, ["GSE13485", "GSE13699"])

    def test_rejects_missing_accession_column(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            links_path = Path(directory) / "geo_series_links.tsv"
            links_path.write_text("study_accession\nSDY1\n", encoding="utf-8")

            with self.assertRaisesRegex(GeoSelectionError, "missing gse_accession"):
                load_gse_ids_from_parsed_links(links_path)

    @patch("data.geo_fetch_module.fetch_geo_datasets", return_value=[])
    def test_main_defaults_to_parsed_links_and_records_provenance(self, fetch_mock) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            links_path = root / "geo_series_links.tsv"
            cache_dir = root / "geo_cache"
            links_path.write_text(
                "study_accession\tgse_accession\n"
                "SDY1\tGSE13485\n"
                "SDY2\tGSE22768\n",
                encoding="utf-8",
            )
            argv = [
                "geo_fetch_module.py",
                "--links-file",
                str(links_path),
                "--destdir",
                str(cache_dir),
            ]

            with patch("sys.argv", argv):
                main()

            provenance = json.loads(
                (cache_dir / "geo_fetch_selection.provenance.json").read_text(
                    encoding="utf-8"
                )
            )

        self.assertEqual(fetch_mock.call_args.args[0], ["GSE13485", "GSE22768"])
        self.assertEqual(provenance["selection_mode"], "parsed_immport_geo_series_links")
        self.assertEqual(provenance["selected_gse_ids"], ["GSE13485", "GSE22768"])
        self.assertIn("sha256", provenance["input"])


if __name__ == "__main__":
    unittest.main()
