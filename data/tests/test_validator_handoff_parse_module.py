"""Focused tests for configured validator handoff table generation."""

from __future__ import annotations

import csv
import gzip
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

from data.validator_handoff_parse_module import run_parser


class ValidatorHandoffParserTests(unittest.TestCase):
    def test_writes_configured_geo_and_immport_records(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            matrix = root / "matrix.txt.gz"
            with gzip.open(matrix, "wt", encoding="utf-8", newline="") as handle:
                handle.write('!Series_geo_accession\t"GSE13485"\n')
                handle.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\n')
                handle.write('!Sample_platform_id\t"GPL7567"\t"GPL7567"\n')
                handle.write("!series_matrix_table_begin\n")
                handle.write('"ID_REF"\t"GSM1"\t"GSM2"\n')
                handle.write('"Hs.412102_at"\t8.049598656\t7.5\n')
                handle.write("!series_matrix_table_end\n")

            archive = root / "study.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.writestr(
                    "study/Tab/lab_test.txt",
                    "\t".join(
                        [
                            "BIOSAMPLE_ACCESSION",
                            "NAME_PREFERRED",
                            "NAME_REPORTED",
                            "RESULT_UNIT_PREFERRED",
                            "RESULT_UNIT_REPORTED",
                            "RESULT_VALUE_PREFERRED",
                            "RESULT_VALUE_REPORTED",
                        ]
                    )
                    + "\nBS1\t\tAct CD8 T Cell Response\tpercentage\tpercentage\t7.68\t7.68\n",
                )
                handle.writestr(
                    "study/Tab/biosample.txt",
                    "\t".join(
                        [
                            "BIOSAMPLE_ACCESSION",
                            "STUDY_ACCESSION",
                            "STUDY_TIME_COLLECTED",
                            "STUDY_TIME_COLLECTED_UNIT",
                            "SUBJECT_ACCESSION",
                        ]
                    )
                    + "\nBS1\tSDY1264\t15\tDays\tSUB1\n",
                )

            manifest = root / "sample_manifest.tsv"
            manifest.write_text(
                "study_accession\tsubject_accession\nSDY1264\tSUB1\n",
                encoding="utf-8",
            )
            output_dir = root / "validator_input"
            geo_output = output_dir / "feature_expression.tsv"
            immport_output = output_dir / "quantitative_outcome.tsv"
            config = {
                "output_dir": str(output_dir),
                "manifest": {"source": str(manifest)},
                "geo": {
                    "study_accession": "SDY1264",
                    "gse_accession": "GSE13485",
                    "gpl_accession": "GPL7567",
                    "gene": "EIF2AK4",
                    "feature_id": "Hs.412102_at",
                    "source_matrix": str(matrix),
                },
                "immport": {
                    "study_accession": "SDY1264",
                    "outcome_names": ["Act CD8 T Cell Response"],
                    "timepoint_source_unit": "Days",
                    "field_precedence": "preferred_then_reported",
                    "source_tab_zip": str(archive),
                },
            }
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")

            result = run_parser(config_path)

            manifest_output = output_dir / "sample_manifest.tsv"
            with geo_output.open(encoding="utf-8", newline="") as handle:
                geo_rows = list(csv.DictReader(handle, delimiter="\t"))
            with immport_output.open(encoding="utf-8", newline="") as handle:
                immport_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(manifest_output.read_bytes(), manifest.read_bytes())
            self.assertEqual(result["manifest"]["selection"]["mode"], "verbatim_copy")
            self.assertEqual(result["geo"]["output"]["rows"], 2)
            self.assertEqual(geo_rows[0]["sample_accession"], "GSM1")
            self.assertEqual(geo_rows[0]["expression_value"], "8.049598656")
            self.assertEqual(result["immport"]["output"]["rows"], 1)
            self.assertEqual(immport_rows[0]["subject_accession"], "SUB1")
            self.assertEqual(immport_rows[0]["timepoint_day"], "15")
            self.assertEqual(immport_rows[0]["result_unit"], "percentage")
            self.assertTrue((output_dir / "sample_manifest.provenance.json").is_file())
            self.assertTrue((output_dir / "feature_expression.provenance.json").is_file())
            self.assertTrue((output_dir / "quantitative_outcome.provenance.json").is_file())
            bundle_path = output_dir / "validator_input_manifest.json"
            bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
            self.assertEqual(
                set(bundle["artifacts"]),
                {"sample_manifest", "feature_expression", "quantitative_outcome"},
            )


if __name__ == "__main__":
    unittest.main()
