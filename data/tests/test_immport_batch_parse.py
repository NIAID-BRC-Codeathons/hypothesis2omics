"""Focused tests for ImmPort batch GEO-series link extraction."""

from __future__ import annotations

import json
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from data.immport_batch_parse import (
    GEO_SERIES_LINKS_FILENAME,
    run_batch,
)
from data.immport_parse_module import OUTPUT_COLUMNS, OUTPUT_FILENAME


def _write_tab_zip(root: Path, accession: str, study_link: str | None) -> None:
    study_dir = root / accession
    study_dir.mkdir()
    archive = study_dir / f"{accession}-DR3_Tab.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        if study_link is not None:
            handle.writestr(
                f"{accession}-DR3_Tab/Tab/study_link.txt",
                study_link,
            )


def _fake_sample_parser(source: Path, output_dir: Path) -> dict:
    accession = source.parent.name
    output_dir.mkdir(parents=True)
    row = {column: "" for column in OUTPUT_COLUMNS}
    row.update(
        {
            "study_accession": accession,
            "experiment_accession": "EXP1",
            "expsample_accession": "ES1",
            "repository_name": "GEO",
            "repository_accession": "GSM1",
        }
    )
    pd.DataFrame([row], columns=OUTPUT_COLUMNS).to_csv(
        output_dir / OUTPUT_FILENAME,
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    return {"counts": {"repository_linked_samples": 1}}


class ImmportBatchParseTests(unittest.TestCase):
    @patch("data.immport_batch_parse.run_parser", side_effect=_fake_sample_parser)
    def test_writes_structured_geo_series_links(self, _parser_mock) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _write_tab_zip(
                root,
                "SDY123",
                "STUDY_LINK_ID\tNAME\tSTUDY_ACCESSION\tTYPE\tVALUE\n"
                "1\tNCBI GEO\tSDY123\t\thttps://example.test/?acc=GSE456\n",
            )

            provenance = run_batch(root)
            links = pd.read_csv(
                root / "parsed" / GEO_SERIES_LINKS_FILENAME,
                sep="\t",
                dtype=str,
                keep_default_na=False,
            )
            batch_manifest = json.loads(
                (root / "parsed" / "batch_manifest.json").read_text(encoding="utf-8")
            )

        self.assertEqual(provenance["status"], "success")
        self.assertEqual(links.loc[0, "study_accession"], "SDY123")
        self.assertEqual(links.loc[0, "gse_accession"], "GSE456")
        self.assertEqual(links.loc[0, "data_release"], "3")
        self.assertEqual(batch_manifest["counts"]["geo_series_accessions"], 1)
        self.assertEqual(batch_manifest["studies"][0]["geo_series_link_status"], "resolved")

    @patch("data.immport_batch_parse.run_parser", side_effect=_fake_sample_parser)
    def test_marks_geo_linked_study_without_gse_as_partial(self, _parser_mock) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _write_tab_zip(
                root,
                "SDY123",
                "STUDY_LINK_ID\tNAME\tSTUDY_ACCESSION\tTYPE\tVALUE\n"
                "1\tPublication\tSDY123\t\thttps://example.test/article\n",
            )

            provenance = run_batch(root)

        self.assertEqual(provenance["status"], "partial")
        self.assertEqual(
            provenance["geo_linked_studies_without_series_link"],
            ["SDY123"],
        )


if __name__ == "__main__":
    unittest.main()
