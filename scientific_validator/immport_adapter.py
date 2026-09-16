#!/usr/bin/env python3

"""
ImmPort metadata adapter for Hypothesis2Omics.

Converts ImmPort study metadata into the evidence format consumed by
eligibility_engine.py.

The adapter is intentionally conservative:
- detecting an assay does not automatically prove the exact required outcome;
- antibody/neutralization responses are kept distinct from CD8 cellular responses;
- missing evidence remains "unknown".
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List


def flatten_text(obj: Any) -> str:
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
    source: str,
) -> Dict[str, Any]:

    return {
        "value": value,
        "confidence": confidence,
        "evidence": evidence,
        "sources": [source],
    }


def get_accession(metadata: Dict[str, Any]) -> str:

    keys = [
        "study_accession",
        "STUDY_ACCESSION",
        "dataset_id",
        "accession",
    ]

    for key in keys:
        if metadata.get(key):
            return str(metadata[key])

    text = flatten_text(metadata)

    match = re.search(r"\bSDY\d+\b", text, flags=re.IGNORECASE)

    if match:
        return match.group(0).upper()

    return "IMM_UNKNOWN"


ASSAY_TERMS = {

    "transcriptomics": [
        "dna microarray",
        "microarray",
        "rna-seq",
        "rna seq",
        "gene expression",
        "transcriptomic",
        "transcriptome",
    ],

    "flow_cytometry": [
        "flow cytometry",
        "facs",
        "mass cytometry",
        "cytof",
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
        "humoral response",
    ],

    "neutralization_assay": [
        "neutralizing antibody",
        "neutralizing antibody titer",
        "neutralisation",
        "neutralization",
        "prnt",
        "plaque reduction",
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
    ],

    "elispot": [
        "elispot",
    ],

    "proteomics": [
        "proteomics",
        "mass spectrometry",
        "protein profiling",
    ],
}


def detect_assays(text: str) -> Dict[str, Any]:

    detected = {}

    for assay, terms in ASSAY_TERMS.items():

        matches = [
            term
            for term in terms
            if term.lower() in text.lower()
        ]

        detected[assay] = {
            "detected": bool(matches),
            "matched_terms": matches,
        }

    return detected


def extract_timepoints(text: str) -> List[int]:
    """
    Extract study day numbers from expressions such as:
    - day 7
    - day 80
    - days 3, 7, 14 and 84
    - day 3/day 7
    """

    text = text.lower()
    values = set()

    # Single expressions such as "day 7" or "d14"
    for match in re.findall(
        r"\b(?:day|days|d)\s*[-_:]?\s*(\d+)\b",
        text,
    ):
        values.add(int(match))

    # Lists such as "days 3, 7, 14 and 84"
    list_matches = re.findall(
        r"\bdays?\s+((?:\d+\s*(?:,|and|/)?\s*)+)",
        text,
    )

    for block in list_matches:
        for number in re.findall(r"\d+", block):
            values.add(int(number))

    return sorted(values)


def infer_immport_evidence(metadata: Dict[str, Any]) -> Dict[str, Any]:

    accession = get_accession(metadata)
    source = f"ImmPort:{accession}"

    full_text = flatten_text(metadata)
    lower = full_text.lower()

    assays = detect_assays(full_text)
    timepoints = extract_timepoints(full_text)

    # ----------------------------------------------------------
    # Human
    # ----------------------------------------------------------

    human = contains_any(
        lower,
        [
            "human",
            "humans",
            "participants",
            "subjects",
            "adults",
            "volunteers",
        ],
    )

    # ----------------------------------------------------------
    # YF-17D
    # ----------------------------------------------------------

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

    # ----------------------------------------------------------
    # Transcriptomics
    # ----------------------------------------------------------

    transcriptomics = assays["transcriptomics"]["detected"]

    # ----------------------------------------------------------
    # EIF2AK4 / GCN2
    # ----------------------------------------------------------

    eif2ak4 = contains_any(
        lower,
        [
            "eif2ak4",
            "gcn2",
        ],
    )

    if eif2ak4:
        eif_value = "true"
        eif_conf = "high"
        eif_message = (
            "EIF2AK4/GCN2 is explicitly referenced in ImmPort metadata."
        )

    elif transcriptomics:
        eif_value = "unknown"
        eif_conf = "medium"
        eif_message = (
            "Transcriptomic data are described, but EIF2AK4 feature "
            "coverage has not yet been verified."
        )

    else:
        eif_value = "unknown"
        eif_conf = "low"
        eif_message = (
            "EIF2AK4 measurability cannot be established."
        )

    # ----------------------------------------------------------
    # Quantitative CD8 response
    # ----------------------------------------------------------

    cd8 = contains_any(
        lower,
        [
            "cd8",
            "cd8+",
            "cd8 t cell",
            "cd8+ t-cell",
            "cytotoxic t cell",
        ],
    )

    quantitative = contains_any(
        lower,
        [
            "frequency",
            "percentage",
            "percent",
            "magnitude",
            "activation",
            "response",
            "cd38",
            "hla-dr",
            "tetramer",
        ],
    )

    flow = assays["flow_cytometry"]["detected"]

    if cd8 and flow and quantitative:

        cd8_value = "true"
        cd8_conf = "high"

        cd8_message = (
            "ImmPort metadata explicitly describes CD8-related "
            "flow-cytometry/quantitative response measurements."
        )

    elif cd8:

        cd8_value = "unknown"
        cd8_conf = "medium"

        cd8_message = (
            "CD8 response terminology is present, but an accessible "
            "participant-level quantitative CD8 endpoint has not yet "
            "been established."
        )

    else:

        cd8_value = "unknown"
        cd8_conf = "low"

        cd8_message = (
            "No explicit quantitative CD8+ T-cell endpoint was identified "
            "in the supplied ImmPort metadata."
        )

    # ----------------------------------------------------------
    # Appropriate timepoints
    #
    # For this hypothesis we need early predictor measurement and
    # a later CD8 outcome. Having many timepoints alone is not enough.
    # ----------------------------------------------------------

    early_time = any(day <= 7 for day in timepoints)
    later_time = any(day >= 14 for day in timepoints)

    if early_time and later_time and cd8_value == "true":

        time_value = "true"
        time_conf = "high"

        time_message = (
            "Early predictor and later quantitative CD8 response "
            "timepoints are represented."
        )

    elif early_time and later_time:

        time_value = "unknown"
        time_conf = "medium"

        time_message = (
            "Early and later study timepoints are present, but the later "
            "timepoint has not yet been verified as the required "
            "quantitative CD8 outcome."
        )

    else:

        time_value = "unknown"
        time_conf = "low"

        time_message = (
            "Required predictor/outcome timepoint relationship "
            "could not be established."
        )

    # ----------------------------------------------------------
    # Participant linkage
    # ----------------------------------------------------------

    linkage_terms = contains_any(
        lower,
        [
            "participant-level",
            "subject-level",
            "same participant",
            "same subject",
            "matched subjects",
        ],
    )

    participant_linkage = evidence_item(
        "true" if linkage_terms else "unknown",
        "high" if linkage_terms else "low",
        (
            "Participant-level linkage is explicitly described."
            if linkage_terms
            else
            "Participant-level linkage between transcriptomic predictor "
            "and immune outcome has not yet been explicitly verified."
        ),
        source,
    )

    # ----------------------------------------------------------
    # Access
    # ----------------------------------------------------------

    predictor_access = evidence_item(
        "true" if transcriptomics else "unknown",
        "medium" if transcriptomics else "low",
        (
            "Transcriptomic data are described in ImmPort study metadata."
            if transcriptomics
            else
            "Predictor data accessibility has not been established."
        ),
        source,
    )

    outcome_access = evidence_item(
        "unknown",
        "low",
        (
            "Presence of an immune-response assay does not by itself prove "
            "that participant-level numeric outcome values are accessible."
        ),
        source,
    )

    # ----------------------------------------------------------
    # Output
    # ----------------------------------------------------------

    return {

        "dataset_id": accession,
        "repository": "ImmPort",

        "criteria": {

            "human_study": evidence_item(
                "true" if human else "unknown",
                "high" if human else "low",
                (
                    "Human participants are described."
                    if human
                    else "Human study status was not established."
                ),
                source,
            ),

            "yf17d_vaccination": evidence_item(
                "true" if yf17d else "unknown",
                "high" if yf17d else "low",
                (
                    "YF-17D/yellow fever vaccination is described."
                    if yf17d
                    else "YF-17D vaccination was not established."
                ),
                source,
            ),

            "transcriptomics_available": evidence_item(
                "true" if transcriptomics else "unknown",
                "high" if transcriptomics else "low",
                (
                    "Transcriptomic assay detected in ImmPort metadata."
                    if transcriptomics
                    else "Transcriptomic assay was not established."
                ),
                source,
            ),

            "eif2ak4_measurable": evidence_item(
                eif_value,
                eif_conf,
                eif_message,
                source,
            ),

            "quantitative_cd8_response_available": evidence_item(
                cd8_value,
                cd8_conf,
                cd8_message,
                source,
            ),

            "appropriate_timepoints": evidence_item(
                time_value,
                time_conf,
                time_message,
                source,
            ),

            "participant_level_linkage": participant_linkage,

            "predictor_data_accessible": predictor_access,

            "outcome_data_accessible": outcome_access,
        },

        "assay_context": {

            "assays_detected": assays,

            "response_modalities": {

                "cellular_response": (
                    assays["flow_cytometry"]["detected"]
                    or assays["cytokine_assay"]["detected"]
                    or assays["elispot"]["detected"]
                ),

                "antibody_response":
                    assays["antibody_response"]["detected"],

                "neutralization_response":
                    assays["neutralization_assay"]["detected"],

                "elisa":
                    assays["elisa"]["detected"],

                "cytokine_response":
                    assays["cytokine_assay"]["detected"],

            },

            "timepoints_detected": timepoints,
        },

        "adapter_notes": [

            "ImmPort assay detection is metadata-based.",

            "Antibody and neutralization responses are kept separate "
            "from CD8 cellular-response evidence.",

            "Unknown fields should remain unresolved until participant-level "
            "data or assay tables are verified.",
        ],
    }


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Convert ImmPort metadata into Hypothesis2Omics "
            "scientific eligibility evidence."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="ImmPort metadata JSON",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional normalized evidence output",
    )

    args = parser.parse_args()

    metadata = json.loads(
        args.input.read_text(encoding="utf-8")
    )

    result = infer_immport_evidence(metadata)

    text = json.dumps(result, indent=2)

    print(text)

    if args.output:

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