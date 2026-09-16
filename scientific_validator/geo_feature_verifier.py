#!/usr/bin/env python3
"""
Verify that a required gene/feature is represented in downloaded GEO data
and update an existing Hypothesis2Omics GEO evidence JSON.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import List


def read_gzip_lines(path: Path) -> List[str]:
    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="ignore",
    ) as handle:
        return list(handle)


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Verify gene/feature measurability from GEO "
            "platform annotation and expression matrix."
        )
    )

    parser.add_argument(
        "--evidence",
        required=True,
        type=Path,
        help="Existing GEO evidence JSON",
    )

    parser.add_argument(
        "--soft",
        required=True,
        type=Path,
        help="GEO family SOFT .gz file",
    )

    parser.add_argument(
        "--matrix",
        required=True,
        type=Path,
        help="GEO series matrix .gz file",
    )

    parser.add_argument(
        "--gene",
        required=True,
        help="Gene symbol, e.g. EIF2AK4",
    )

    parser.add_argument(
        "--feature",
        required=True,
        help="Platform feature/probe identifier",
    )

    parser.add_argument(
        "--annotation-id",
        default=None,
        help="Optional supporting accession such as GenBank ID",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Updated evidence JSON",
    )

    args = parser.parse_args()

    evidence = json.loads(
        args.evidence.read_text(
            encoding="utf-8"
        )
    )

    soft_lines = read_gzip_lines(
        args.soft
    )

    matrix_lines = read_gzip_lines(
        args.matrix
    )

    gene_pattern = re.compile(
        rf"(?<![A-Za-z0-9_])"
        rf"{re.escape(args.gene)}"
        rf"(?![A-Za-z0-9_])",
        re.IGNORECASE,
    )

    annotation_matches = []

    for line in soft_lines:

        feature_match = (
            args.feature.lower()
            in line.lower()
        )

        gene_match = bool(
            gene_pattern.search(line)
        )

        annotation_id_match = (
            args.annotation_id is not None
            and args.annotation_id.lower()
            in line.lower()
        )

        if (
            feature_match
            or gene_match
            or annotation_id_match
        ):
            annotation_matches.append(
                line.strip()
            )

    matrix_matches = [
        line.strip()
        for line in matrix_lines
        if args.feature.lower()
        in line.lower()
    ]

    feature_annotated = bool(
        annotation_matches
    )

    feature_in_matrix = bool(
        matrix_matches
    )

    measurable = (
        feature_annotated
        and feature_in_matrix
    )

    dataset_id = evidence.get(
        "dataset_id",
        "GEO_UNKNOWN",
    )

    source = (
        f"GEO:{dataset_id}:"
        f"{args.feature}"
    )

    if measurable:

        evidence["criteria"][
            "eif2ak4_measurable"
        ] = {
            "value": "true",
            "confidence": "high",
            "evidence": (
                f"{args.gene} is represented by feature "
                f"{args.feature} in the GEO platform annotation "
                "and the same feature has numeric expression "
                "values in the downloaded series matrix."
            ),
            "sources": [
                source,
            ],
        }

    else:

        evidence["criteria"][
            "eif2ak4_measurable"
        ] = {
            "value": "unknown",
            "confidence": "medium",
            "evidence": (
                f"Could not verify both platform annotation and "
                f"matrix-level expression for {args.gene} "
                f"({args.feature})."
            ),
            "sources": [
                source,
            ],
        }

    evidence[
        "feature_verification"
    ] = {
        "gene": args.gene,
        "feature_id": args.feature,
        "annotation_id": args.annotation_id,
        "platform_annotation_found":
            feature_annotated,
        "expression_matrix_feature_found":
            feature_in_matrix,
        "measurable": measurable,
        "annotation_match_count":
            len(annotation_matches),
        "matrix_match_count":
            len(matrix_matches),
    }

    text = json.dumps(
        evidence,
        indent=2,
    )

    print(text)

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.output.write_text(
        text + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()