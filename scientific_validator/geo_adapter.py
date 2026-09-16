#!/usr/bin/env python3
"""
GEO metadata adapter for Hypothesis2Omics.

Purpose
-------
Convert raw/normalized GEO study metadata into the evidence structure expected
by eligibility_engine.py.

Important:
- This adapter is intentionally conservative.
- Detecting an assay (e.g. flow cytometry) does NOT automatically prove that
  the exact hypothesis outcome is available.
- Exact predictor/outcome availability remains "unknown" unless the metadata
  contains sufficiently explicit evidence.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def flatten_text(obj: Any) -> str:
    """Recursively convert metadata into one searchable text string."""
    parts: List[str] = []

    if isinstance(obj, dict):
        for key, value in obj.items():
            parts.append(str(key))
            parts.append(flatten_text(value))

    elif isinstance(obj, list):
        for value in obj:
            parts.append(flatten_text(value))

    elif obj is not None:
        parts.append(str(obj))

    return " ".join(parts)


def contains_any(text: str, terms: List[str]) -> bool:
    text = text.lower()
    return any(term.lower() in text for term in terms)


def evidence_item(
    value: str,
    confidence: str,
    evidence: str,
    source: str
) -> Dict[str, Any]:
    return {
        "value": value,
        "confidence": confidence,
        "evidence": evidence,
        "sources": [source],
    }


def get_accession(metadata: Dict[str, Any]) -> str:
    """Try common GEO accession field names."""
    possible_keys = [
        "accession",
        "gse_accession",
        "series_accession",
        "dataset_id",
        "geo_accession",
    ]

    for key in possible_keys:
        if metadata.get(key):
            return str(metadata[key])

    # Search values for a GSE accession
    text = flatten_text(metadata)
    match = re.search(r"\bGSE\d+\b", text, flags=re.IGNORECASE)

    if match:
        return match.group(0).upper()

    return "GEO_UNKNOWN"


# ---------------------------------------------------------------------
# Assay detection
# ---------------------------------------------------------------------

ASSAY_TERMS = {
    "transcriptomics": [
        "rna-seq",
        "rna seq",
        "expression profiling",
        "gene expression",
        "microarray",
        "transcriptomic",
        "transcriptome",
        "mrna",
    ],

    "flow_cytometry": [
        "flow cytometry",
        "facs",
        "cytometry",
        "cd38",
        "hla-dr",
        "tetramer",
    ],

    "antibody_response": [
        "antibody",
        "antibodies",
        "igg",
        "igm",
        "serology",
        "humoral",
    ],

    "neutralization_assay": [
        "neutralizing antibody",
        "neutralisation",
        "neutralization",
        "prnt",
        "plaque reduction neutralization",
    ],

    "elisa": [
        "elisa",
        "enzyme-linked immunosorbent",
    ],

    "cytokine_assay": [
        "cytokine",
        "ifn-gamma",
        "ifng",
        "tnf",
        "il-2",
        "interleukin",
        "elispot",
    ],

    "proteomics": [
        "proteomics",
        "mass spectrometry",
        "protein profiling",
    ],

    "single_cell": [
        "single-cell",
        "single cell",
        "scrna-seq",
        "scrna seq",
    ],
}


def detect_assays(text: str) -> Dict[str, Any]:
    detected = {}

    for assay, terms in ASSAY_TERMS.items():
        matched = [
            term for term in terms
            if term.lower() in text.lower()
        ]

        detected[assay] = {
            "detected": bool(matched),
            "matched_terms": matched,
        }

    return detected


# ---------------------------------------------------------------------
# Evidence inference
# ---------------------------------------------------------------------

def infer_geo_evidence(metadata: Dict[str, Any]) -> Dict[str, Any]:

    accession = get_accession(metadata)
    source = f"GEO:{accession}"

    full_text = flatten_text(metadata)
    lower = full_text.lower()

    assays = detect_assays(full_text)

    # --------------------------------------------------------------
    # Human study
    # --------------------------------------------------------------

    human = contains_any(
        lower,
        [
            "homo sapiens",
            "human",
            "humans",
            "participants",
            "subjects",
            "volunteers",
        ],
    )

    human_evidence = evidence_item(
        "true" if human else "unknown",
        "high" if human else "low",
        (
            "Human study terminology detected in GEO metadata."
            if human
            else "Human study status could not be established from GEO metadata."
        ),
        source,
    )

    # --------------------------------------------------------------
    # YF-17D vaccination
    # --------------------------------------------------------------

    yf17d = contains_any(
        lower,
        [
            "yf-17d",
            "yf17d",
            "yellow fever vaccine",
            "yellow fever vaccination",
            "stamaril",
        ],
    )

    yf17d_evidence = evidence_item(
        "true" if yf17d else "unknown",
        "high" if yf17d else "low",
        (
            "YF-17D/yellow fever vaccination terminology detected."
            if yf17d
            else "YF-17D vaccination could not be established from GEO metadata."
        ),
        source,
    )

    # --------------------------------------------------------------
    # Transcriptomics
    # --------------------------------------------------------------

    transcriptomics = assays["transcriptomics"]["detected"]

    transcriptomics_evidence = evidence_item(
        "true" if transcriptomics else "unknown",
        "high" if transcriptomics else "low",
        (
            "Transcriptomic/gene-expression assay detected in GEO metadata."
            if transcriptomics
            else "Transcriptomic assay could not be established."
        ),
        source,
    )

    # --------------------------------------------------------------
    # EIF2AK4 measurable
    #
    # Conservative:
    # presence of transcriptomics alone does not prove feature coverage.
    # --------------------------------------------------------------

    eif2ak4_explicit = contains_any(
        lower,
        [
            "eif2ak4",
            "gcn2",
        ],
    )

    if eif2ak4_explicit:
        eif2ak4_value = "true"
        eif2ak4_confidence = "high"
        eif2ak4_message = (
            "EIF2AK4/GCN2 is explicitly referenced in the GEO metadata."
        )
    elif transcriptomics:
        eif2ak4_value = "unknown"
        eif2ak4_confidence = "medium"
        eif2ak4_message = (
            "Transcriptomics are available, but EIF2AK4 feature coverage "
            "has not yet been verified from platform/gene annotation."
        )
    else:
        eif2ak4_value = "unknown"
        eif2ak4_confidence = "low"
        eif2ak4_message = (
            "EIF2AK4 measurability cannot currently be established."
        )

    eif2ak4_evidence = evidence_item(
        eif2ak4_value,
        eif2ak4_confidence,
        eif2ak4_message,
        source,
    )

    # --------------------------------------------------------------
    # CD8 outcome candidate
    #
    # Detecting CD8/flow data makes this a candidate, but DOES NOT
    # automatically prove that a participant-level quantitative endpoint
    # is actually available.
    # --------------------------------------------------------------

    cd8_terms = contains_any(
        lower,
        [
            "cd8",
            "cd8+",
            "cd8 t cell",
            "cd8+ t cell",
            "cytotoxic t cell",
        ],
    )

    quantitative_terms = contains_any(
        lower,
        [
            "frequency",
            "percentage",
            "percent",
            "magnitude",
            "response",
            "activation",
            "cd38",
            "hla-dr",
            "tetramer",
        ],
    )

    flow_detected = assays["flow_cytometry"]["detected"]

    if cd8_terms and flow_detected and quantitative_terms:
        cd8_value = "unknown"
        cd8_confidence = "medium"
        cd8_message = (
            "GEO metadata contains CD8-related and flow-cytometry/quantitative "
            "response terminology. This suggests a candidate CD8 outcome, "
            "but participant-level quantitative values must still be verified."
        )

    elif cd8_terms:
        cd8_value = "unknown"
        cd8_confidence = "low"
        cd8_message = (
            "CD8 terminology is present, but the exact quantitative "
            "CD8 response endpoint has not been established."
        )

    else:
        cd8_value = "unknown"
        cd8_confidence = "low"
        cd8_message = (
            "A quantitative CD8+ T-cell response could not be established "
            "from GEO metadata."
        )

    cd8_evidence = evidence_item(
        cd8_value,
        cd8_confidence,
        cd8_message,
        source,
    )

    # --------------------------------------------------------------
    # Timepoints
    # --------------------------------------------------------------

    pre_terms = contains_any(
        lower,
        [
            "baseline",
            "pre-vaccination",
            "prevaccination",
            "day 0",
            "d0",
        ],
    )

    post_terms = bool(
        re.search(
            r"\b(day|d)\s*[-_ ]?(1|2|3|4|5|6|7|8|9|10|14|15|21|28|30|60|80|84)\b",
            lower,
        )
    )

    if pre_terms and post_terms:
        time_value = "true"
        time_confidence = "high"
        time_message = (
            "Both baseline/pre-vaccination and post-vaccination timepoints "
            "are represented in GEO metadata."
        )

    elif post_terms:
        time_value = "unknown"
        time_confidence = "medium"
        time_message = (
            "Post-vaccination timepoints are present, but required predictor "
            "and outcome timing has not been fully established."
        )

    else:
        time_value = "unknown"
        time_confidence = "low"
        time_message = (
            "Appropriate predictor/outcome timepoints could not be established."
        )

    time_evidence = evidence_item(
        time_value,
        time_confidence,
        time_message,
        source,
    )

    # --------------------------------------------------------------
    # Participant-level linkage
    #
    # GEO sample identifiers do not automatically prove linkage to an
    # immune phenotype stored elsewhere.
    # --------------------------------------------------------------

    participant_evidence = evidence_item(
        "unknown",
        "low",
        (
            "GEO samples may contain subject identifiers, but participant-level "
            "linkage between transcriptomic predictor and quantitative immune "
            "outcome must be verified explicitly."
        ),
        source,
    )

    # --------------------------------------------------------------
    # Data access
    # --------------------------------------------------------------

    predictor_access = evidence_item(
        "true" if transcriptomics else "unknown",
        "high" if transcriptomics else "low",
        (
            "GEO transcriptomic data are represented in the study metadata."
            if transcriptomics
            else "Predictor data accessibility could not be established."
        ),
        source,
    )

    outcome_access = evidence_item(
        "unknown",
        "low",
        (
            "Exact participant-level quantitative CD8 outcome values have not "
            "yet been verified in GEO. A linked source such as ImmPort may "
            "provide complementary immune phenotype data."
        ),
        source,
    )

    # --------------------------------------------------------------
    # Final normalized record
    # --------------------------------------------------------------

    return {
        "dataset_id": accession,
        "repository": "GEO",

        "criteria": {
            "human_study": human_evidence,
            "yf17d_vaccination": yf17d_evidence,
            "transcriptomics_available": transcriptomics_evidence,
            "eif2ak4_measurable": eif2ak4_evidence,
            "quantitative_cd8_response_available": cd8_evidence,
            "appropriate_timepoints": time_evidence,
            "participant_level_linkage": participant_evidence,
            "predictor_data_accessible": predictor_access,
            "outcome_data_accessible": outcome_access,
        },

        "assay_context": {
            "assays_detected": assays,

            "response_modalities": {
                "cellular_response": (
                    assays["flow_cytometry"]["detected"]
                    or assays["cytokine_assay"]["detected"]
                ),
                "antibody_response": assays["antibody_response"]["detected"],
                "neutralization_response": assays["neutralization_assay"]["detected"],
                "elisa": assays["elisa"]["detected"],
            },
        },

        "adapter_notes": [
            "Assay detection is metadata-based and does not prove availability "
            "of participant-level numeric outcome values.",

            "Unknown values should be resolved using linked repositories such "
            "as ImmPort, publications, supplementary files, or assay tables.",
        ],
    }


# ---------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------

def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Convert GEO metadata into Hypothesis2Omics eligibility evidence."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Raw or normalized GEO metadata JSON file",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional normalized evidence JSON output file",
    )

    args = parser.parse_args()

    metadata = json.loads(args.input.read_text(encoding="utf-8"))

    result = infer_geo_evidence(metadata)

    text = json.dumps(result, indent=2)

    print(text)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()