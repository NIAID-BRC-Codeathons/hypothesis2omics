"""Focused tests for the Galaxy neutralizing-antibody result export."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

from data.immport_neut_ab_export import (
    OUTPUT_FILENAME,
    PROVENANCE_FILENAME,
    NeutralizingAntibodyExportError,
    run_export,
)


COLUMNS = [
    "RESULT_ID",
    "ARM_ACCESSION",
    "BIOSAMPLE_ACCESSION",
    "COMMENTS",
    "EXPERIMENT_ACCESSION",
    "EXPSAMPLE_ACCESSION",
    "REPOSITORY_ACCESSION",
    "REPOSITORY_NAME",
    "STUDY_ACCESSION",
    "STUDY_TIME_COLLECTED",
    "STUDY_TIME_COLLECTED_UNIT",
    "SUBJECT_ACCESSION",
    "UNIT_PREFERRED",
    "UNIT_REPORTED",
    "VALUE_PREFERRED",
    "VALUE_REPORTED",
    "VIRUS_STRAIN_PREFERRED",
    "VIRUS_STRAIN_REPORTED",
    "WORKSPACE_ID",
]


def _write_study(root: Path, accession: str, release: int, row: list[str]) -> None:
    study_dir = root / accession
    study_dir.mkdir()
    archive = study_dir / f"{accession}-DR{release}_Tab.zip"
    content = "\t".join(COLUMNS) + "\n" + "\t".join(row) + "\n"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(
            f"{accession}-DR{release}_Tab/Tab/neut_ab_titer_result.txt",
            content,
        )


class NeutralizingAntibodyExportTests(unittest.TestCase):
    def test_exports_raw_rows_and_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            output_dir = root / "galaxy"
            cache.mkdir()
            _write_study(
                cache,
                "SDY2",
                3,
                [
                    "2",
                    "ARM2",
                    "BS2",
                    "",
                    "EXP2",
                    "ES2",
                    "",
                    "",
                    "SDY2",
                    "28",
                    "Days",
                    "SUB2",
                    "Antibody titer",
                    "Antibody titer",
                    "",
                    "<10",
                    "Yellow fever virus 17D",
                    "YF17D",
                    "20",
                ],
            )
            _write_study(
                cache,
                "SDY1",
                4,
                [
                    "1",
                    "ARM1",
                    "BS1",
                    "",
                    "EXP1",
                    "ES1",
                    "",
                    "",
                    "SDY1",
                    "60",
                    "Days",
                    "SUB1",
                    "Antibody titer",
                    "Antibody titer",
                    "320",
                    "320",
                    "Yellow fever virus 17D",
                    "YF17D",
                    "10",
                ],
            )

            provenance = run_export(cache, output_dir)

            output = output_dir / OUTPUT_FILENAME
            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            sidecar = json.loads(
                (output_dir / PROVENANCE_FILENAME).read_text(encoding="utf-8")
            )

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["study_accession"], "SDY1")
        self.assertEqual(rows[0]["value_preferred"], "320")
        self.assertEqual(rows[1]["value_preferred"], "")
        self.assertEqual(rows[1]["value_reported"], "<10")
        self.assertEqual(provenance["counts"]["value_preferred_populated"], 1)
        self.assertEqual(provenance["counts"]["value_preferred_blank"], 1)
        self.assertEqual(sidecar["output"]["sha256"], provenance["output"]["sha256"])

    def test_rejects_missing_required_column(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cache = root / "cache"
            study = cache / "SDY1"
            study.mkdir(parents=True)
            archive = study / "SDY1-DR1_Tab.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr(
                    "SDY1-DR1_Tab/Tab/neut_ab_titer_result.txt",
                    "RESULT_ID\tSTUDY_ACCESSION\n1\tSDY1\n",
                )

            with self.assertRaisesRegex(
                NeutralizingAntibodyExportError,
                "missing required columns",
            ):
                run_export(cache, root / "output")


if __name__ == "__main__":
    unittest.main()
