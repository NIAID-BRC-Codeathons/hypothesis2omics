"""Tests for deterministic downloaded-file inventory and parser classification."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
MCP_ROOT = REPOSITORY_ROOT / "mcp"
if str(MCP_ROOT) not in sys.path:
    sys.path.insert(0, str(MCP_ROOT))

from data.list_directory import build_file_inventory  # noqa: E402
from data_normalizer.geo_tools import DEFAULT_PARSED_ROOT  # noqa: E402
from data_normalizer.pipeline_tools import DEFAULT_DATA_ROOT, PipelineTools  # noqa: E402


class FileInventoryTests(unittest.TestCase):
    def test_data_normalizer_defaults_resolve_from_repository_root(self) -> None:
        self.assertEqual(DEFAULT_DATA_ROOT, REPOSITORY_ROOT / "data")
        self.assertEqual(
            DEFAULT_PARSED_ROOT,
            REPOSITORY_ROOT / "data" / "geo_cache" / "parsed",
        )

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.data_root = Path(self.temporary_directory.name) / "data"
        self.immport_study = self.data_root / "immport_cache" / "SDY123"
        self.geo_study = self.data_root / "geo_cache" / "GSE456"
        self.immport_study.mkdir(parents=True)
        self.geo_study.mkdir(parents=True)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def _rows(path: str | Path) -> list[dict[str, str]]:
        with Path(path).open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_classifies_candidates_and_excludes_credentials_and_generated_files(self) -> None:
        tab_archive = self.immport_study / "SDY123-DR1_Tab.zip"
        tab_archive.write_bytes(b"tab archive")
        (self.immport_study / "notes.tsv").write_text("a\tb\n", encoding="utf-8")
        key_path = self.data_root / "immport_cache" / "immport-key-secret.json"
        key_path.write_text("secret", encoding="utf-8")
        parsed = self.data_root / "geo_cache" / "parsed"
        parsed.mkdir()
        (parsed / "generated.tsv").write_text("generated", encoding="utf-8")
        matrix = self.geo_study / "GSE456_series_matrix.txt.gz"
        matrix.write_bytes(b"matrix")
        soft = self.geo_study / "GSE456_family.soft.gz"
        soft.write_bytes(b"soft")

        result = build_file_inventory(self.data_root)
        inventory = self._rows(result["outputs"][0]["path"])
        by_name = {row["file_name"]: row for row in inventory}

        self.assertEqual(by_name[tab_archive.name]["candidate_type"], "immport_tab_archive")
        self.assertEqual(by_name[matrix.name]["candidate_type"], "geo_series_matrix")
        self.assertEqual(by_name[soft.name]["candidate_type"], "geo_family_soft")
        self.assertEqual(by_name["notes.tsv"]["parser_support"], "potential")
        self.assertNotIn(key_path.name, by_name)
        self.assertNotIn("generated.tsv", by_name)

    def test_inventory_is_sorted_and_reconciles_portable_manifest_paths(self) -> None:
        later = self.immport_study / "z_Tab.zip"
        earlier = self.immport_study / "a_Tab.zip"
        later.write_bytes(b"z")
        earlier.write_bytes(b"a")
        manifest = [
            {
                "sdy_id": "SDY123",
                "artifacts": [
                    {
                        "artifact_type": "study_file",
                        "status": "success",
                        "source_path": "/SDY123/a_Tab.zip",
                        "local_file": "/different/machine/SDY123/a_Tab.zip",
                        "sha256": "abc123",
                    }
                ],
            }
        ]
        manifest_path = self.data_root / "immport_cache" / "manifest.json"
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

        result = build_file_inventory(self.data_root)
        inventory = self._rows(result["outputs"][0]["path"])
        study_rows = [row for row in inventory if row["dataset_accession"] == "SDY123"]

        self.assertEqual(
            [row["relative_path"] for row in study_rows],
            ["SDY123/a_Tab.zip", "SDY123/z_Tab.zip"],
        )
        self.assertEqual(study_rows[0]["manifest_recorded"], "true")
        self.assertEqual(study_rows[0]["sha256"], "abc123")

    def test_symlinks_are_not_followed_or_inventoried(self) -> None:
        outside = Path(self.temporary_directory.name) / "outside.tsv"
        outside.write_text("outside", encoding="utf-8")
        (self.immport_study / "linked.tsv").symlink_to(outside)

        result = build_file_inventory(self.data_root)
        inventory = self._rows(result["outputs"][0]["path"])

        self.assertFalse(any(row["file_name"] == "linked.tsv" for row in inventory))
        self.assertEqual(result["counts"]["skipped"]["symlink"], 1)

    def test_pipeline_wrapper_uses_configured_data_root(self) -> None:
        (self.geo_study / "GSE456_family.soft.gz").write_bytes(b"soft")

        result = PipelineTools(self.data_root).inventory_downloaded_files()

        self.assertEqual(result["status"], "success")
        self.assertEqual(len(result["supported_candidates"]), 1)
        resolved_root = self.data_root.resolve()
        self.assertTrue(Path(result["inventory_path"]).is_relative_to(resolved_root))
        self.assertTrue(Path(result["provenance_path"]).is_relative_to(resolved_root))


if __name__ == "__main__":
    unittest.main()
