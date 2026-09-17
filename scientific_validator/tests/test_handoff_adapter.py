"""Focused tests for scientific-validator handoff input path resolution."""

from pathlib import Path
import unittest

from scientific_validator.handoff_adapter import resolve_input_paths


class HandoffAdapterTests(unittest.TestCase):
    def test_resolves_canonical_bundle_paths(self) -> None:
        paths = resolve_input_paths("data/validator_input", None, None, None)

        expected = (
            str(Path("data/validator_input/sample_manifest.tsv")),
            str(Path("data/validator_input/feature_expression.tsv")),
            str(Path("data/validator_input/quantitative_outcome.tsv")),
        )

        self.assertEqual(paths, expected)

    def test_explicit_path_overrides_bundle_path(self) -> None:
        paths = resolve_input_paths(
            "data/validator_input",
            "custom/manifest.tsv",
            None,
            None,
        )

        self.assertEqual(paths[0], "custom/manifest.tsv")

    def test_preserves_existing_explicit_path_interface(self) -> None:
        paths = resolve_input_paths(
            None,
            "manifest.tsv",
            "features.tsv",
            "outcome.tsv",
        )

        self.assertEqual(paths, ("manifest.tsv", "features.tsv", "outcome.tsv"))


if __name__ == "__main__":
    unittest.main()