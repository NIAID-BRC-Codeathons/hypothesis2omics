import json
import unittest
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from eligibility_engine import assess_dataset


class EligibilityEngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = json.loads((ROOT / "examples" / "yf17d_test_spec.json").read_text())

    def load(self, name):
        return json.loads((ROOT / "examples" / name).read_text())

    def test_known_positive_is_scientifically_eligible(self):
        result = assess_dataset(self.spec, self.load("SDY1264.json"))
        self.assertEqual(result["scientific_status"], "eligible")
        self.assertEqual(result["execution_status"], "needs_data_check")
        self.assertEqual(result["unresolved_criteria"], ["outcome_data_accessible"])

    def test_unresolved_science_is_not_promoted_to_eligible(self):
        result = assess_dataset(self.spec, self.load("SDY1529.json"))
        self.assertEqual(result["scientific_status"], "uncertain")
        self.assertEqual(result["execution_status"], "needs_scientific_resolution")
        self.assertIn("quantitative_cd8_response_available", result["unresolved_criteria"])

    def test_failed_required_criterion_excludes_dataset(self):
        result = assess_dataset(self.spec, self.load("synthetic_excluded.json"))
        self.assertEqual(result["scientific_status"], "excluded")
        self.assertEqual(result["execution_status"], "not_applicable")
        self.assertEqual(result["failed_scientific_criteria"], ["yf17d_vaccination"])

    def test_evidence_provenance_is_preserved(self):
        result = assess_dataset(self.spec, self.load("SDY1264.json"))
        eif = next(x for x in result["criterion_assessments"] if x["criterion"] == "eif2ak4_measurable")
        self.assertIn("PMID:19029902", eif["sources"])
        self.assertEqual(eif["confidence"], "high")


if __name__ == "__main__":
    unittest.main()
