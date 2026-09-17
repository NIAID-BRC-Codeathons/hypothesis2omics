"""Focused tests for ImmPort download progress and CLI failure reporting."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from data.immport_fetch_module import (
    ArtifactRecord,
    FetchRecord,
    ImmportSession,
    _download_progress,
    _format_bytes,
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
                )

        messages = "\n".join(captured.output)
        self.assertEqual(result.status, "partial")
        self.assertIn("Manifest lists 2 files; processing 2", messages)
        self.assertIn("first.zip: downloaded", messages)
        self.assertIn("second.zip: failed: HTTP 500", messages)

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
