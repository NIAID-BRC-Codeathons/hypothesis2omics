"""Focused tests for ImmPort file selection, progress, and CLI failure reporting."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from data.immport_batch_parse import discover_tab_sources
from data.immport_fetch_module import (
    ArtifactRecord,
    FetchRecord,
    ImmportSession,
    _download_progress,
    _format_bytes,
    _select_parser_required_entries,
    fetch_one,
    main,
)


def _artifact(status: str, path: str, error: str | None = None) -> ArtifactRecord:
    return ArtifactRecord(
        artifact_type="study_file",
        status=status,
        source_path=path,
        local_file=path,
        retrieval_tool="test",
        retrieval_tool_version="test",
        size_bytes=1024 if status != "failed" else None,
        error=error,
    )


class ImmportFetchProgressTests(unittest.TestCase):
    def test_formats_download_progress(self) -> None:
        self.assertEqual(_format_bytes(1024), "1.0 KiB")
        self.assertEqual(_download_progress(1024, 2048), "1.0 KiB/2.0 KiB")

    def test_selects_one_highest_release_direct_tab_archive(self) -> None:
        entries = [
            {"path": "/SDY123/SDY123-DR57_Tab.zip"},
            {"path": "/SDY123/SDY123-DR58_MySQL.zip"},
            {"path": "/SDY123/ResultFiles/SDY123-DR99_Tab.zip"},
            {"path": "/SDY123/SDY123-DR58_Tab.zip"},
        ]

        selected, release = _select_parser_required_entries("SDY123", entries)

        self.assertEqual(release, 58)
        self.assertEqual(selected, [entries[-1]])

    def test_rejects_ambiguous_highest_release_tab_archive(self) -> None:
        entries = [
            {"path": "/SDY123/SDY123-DR58_Tab.zip", "fileUUID": "first"},
            {"path": "/SDY123/SDY123-DR58_Tab.zip", "fileUUID": "second"},
        ]

        with self.assertRaisesRegex(ValueError, "Expected one .* but found 2"):
            _select_parser_required_entries("SDY123", entries)

    @patch("data.immport_fetch_module._download_file")
    @patch("data.immport_fetch_module._fetch_filepath_manifest")
    def test_fetch_one_logs_file_statuses(self, manifest_mock, download_mock) -> None:
        manifest_mock.return_value = [
            {
                "path": "/SDY123/StudyFiles/first.zip",
                "fileUUID": "first",
                "filesizeBytes": 1024,
            },
            {
                "path": "/SDY123/StudyFiles/second.zip",
                "fileUUID": "second",
                "filesizeBytes": 2048,
            },
        ]
        download_mock.side_effect = [
            _artifact("success", "first.zip"),
            _artifact("failed", "second.zip", "HTTP 500"),
        ]

        with tempfile.TemporaryDirectory() as directory:
            with self.assertLogs("immport_fetcher", level="INFO") as captured:
                result = fetch_one(
                    "SDY123",
                    Path(directory),
                    ImmportSession(api_key="test-key"),
                    all_files=True,
                )

        messages = "\n".join(captured.output)
        self.assertEqual(result.status, "partial")
        self.assertIn("Manifest lists 2 files; selected 2 (all_files); processing 2", messages)
        self.assertIn("first.zip: downloaded", messages)
        self.assertIn("second.zip: failed: HTTP 500", messages)

    @patch("data.immport_fetch_module._download_file")
    @patch("data.immport_fetch_module._fetch_filepath_manifest")
    def test_default_fetch_downloads_only_parser_tab_zip(
        self,
        manifest_mock,
        download_mock,
    ) -> None:
        manifest_mock.return_value = [
            {
                "path": "/SDY123/ResultFiles/result.txt",
                "fileUUID": "result",
                "filesizeBytes": 20,
            },
            {
                "path": "/SDY123/SDY123-DR58_Tab.zip",
                "fileUUID": "tab",
                "filesizeBytes": 1024,
            },
        ]
        download_mock.return_value = _artifact("success", "SDY123-DR58_Tab.zip")

        with tempfile.TemporaryDirectory() as directory:
            result = fetch_one(
                "SDY123",
                Path(directory),
                ImmportSession(api_key="test-key"),
            )
            tab_path = Path(directory) / "SDY123" / "SDY123-DR58_Tab.zip"
            tab_path.touch()
            discovered = discover_tab_sources(directory)

        self.assertEqual(download_mock.call_count, 1)
        self.assertEqual(download_mock.call_args.args[0], "/SDY123/SDY123-DR58_Tab.zip")
        self.assertEqual(result.status, "success")
        self.assertEqual(result.selection_mode, "parser_required")
        self.assertEqual(result.selected_data_release, 58)
        self.assertEqual(result.n_files_listed, 2)
        self.assertEqual(result.n_files_selected, 1)
        self.assertEqual(result.n_files_downloaded, 1)
        self.assertEqual(discovered["SDY123"][0], 58)

    @patch("data.immport_fetch_module._download_file")
    @patch("data.immport_fetch_module._fetch_filepath_manifest")
    def test_missing_tab_archive_fails_with_manifest_provenance(
        self,
        manifest_mock,
        download_mock,
    ) -> None:
        manifest_mock.return_value = [
            {
                "path": "/SDY123/ResultFiles/result.txt",
                "fileUUID": "result",
            }
        ]

        with tempfile.TemporaryDirectory() as directory:
            result = fetch_one(
                "SDY123",
                Path(directory),
                ImmportSession(api_key="test-key"),
            )
            manifest_path = Path(directory) / "SDY123" / "SDY123_filepath_manifest.json"
            self.assertTrue(manifest_path.is_file())

        download_mock.assert_not_called()
        self.assertEqual(result.status, "failed")
        self.assertEqual(result.n_files_selected, 0)
        self.assertIn("No direct SDY123-DR<n>_Tab.zip", result.error or "")
        self.assertEqual(result.artifacts[0].artifact_type, "filepath_manifest")

    @patch("data.immport_fetch_module.fetch_immport_datasets")
    @patch("data.immport_fetch_module._parse_args")
    def test_cli_writes_manifest_before_nonzero_exit(self, args_mock, fetch_mock) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args_mock.return_value = SimpleNamespace(
                sdy_ids=["SDY123"],
                destdir=directory,
                api_key_file="key.json",
                force=False,
                max_files_per_study=None,
                all_files=False,
            )
            fetch_mock.return_value = [
                FetchRecord(
                    sdy_id="SDY123",
                    status="partial",
                    run_started_at_utc="2026-01-01T00:00:00+00:00",
                    destdir=directory,
                    api_base_url="https://example.test",
                    filepath_manifest_url="https://example.test/manifest",
                    error="one file failed",
                )
            ]

            with self.assertRaisesRegex(SystemExit, "1"):
                main()

            payload = json.loads((Path(directory) / "manifest.json").read_text())
            self.assertEqual(payload[0]["status"], "partial")


if __name__ == "__main__":
    unittest.main()
