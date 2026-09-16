#!/usr/bin/env python3
"""
Cross-repository evidence fusion for Hypothesis2Omics.

Combines normalized evidence from sources such as GEO and ImmPort before
sending the fused result to eligibility_engine.py.

Safety model
------------
- Cross-repository records are fused only when their relationship is present
  in a verified linkage registry.
- TRUE + UNKNOWN -> TRUE
- FALSE + UNKNOWN -> FALSE
- TRUE + FALSE -> UNKNOWN with a conflict flag
- UNKNOWN + UNKNOWN -> UNKNOWN
- Fusion does not by itself prove participant-level predictor/outcome linkage.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional


CONFIDENCE_RANK = {
    "low": 1,
    "medium": 2,
    "high": 3,
}

RANK_CONFIDENCE = {
    1: "low",
    2: "medium",
    3: "high",
}


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_value(value: Any) -> str:
    value = str(value).lower().strip()

    if value in {"true", "yes", "1"}:
        return "true"

    if value in {"false", "no", "0"}:
        return "false"

    return "unknown"


def confidence_rank(confidence: str) -> int:
    return CONFIDENCE_RANK.get(
        str(confidence).lower(),
        1,
    )


def highest_confidence(items: List[Dict[str, Any]]) -> str:
    if not items:
        return "low"

    rank = max(
        confidence_rank(item.get("confidence", "low"))
        for item in items
    )

    return RANK_CONFIDENCE[rank]


def unique_list(values: List[Any]) -> List[Any]:
    result = []

    for value in values:
        if value not in result:
            result.append(value)

    return result


def infer_repository(
    dataset_id: str,
    repository: Optional[str] = None,
) -> str:
    """
    Preserve an explicit repository when available.
    Otherwise infer from common accession prefixes.
    """
    if repository and str(repository).lower() not in {"", "unknown", "none"}:
        return str(repository)

    dataset_id = str(dataset_id).upper()

    if dataset_id.startswith("GSE") or dataset_id.startswith("GSM"):
        return "GEO"

    if dataset_id.startswith("SDY"):
        return "ImmPort"

    return "unknown"


# ---------------------------------------------------------------------
# Linkage registry
# ---------------------------------------------------------------------

def load_link_registry(path: Path) -> Dict[str, Any]:
    registry = load_json(path)

    if "links" not in registry or not isinstance(registry["links"], list):
        raise ValueError(
            "Linkage registry must contain a top-level 'links' list."
        )

    return registry


def linkage_is_allowed(
    documents: List[Dict[str, Any]],
    registry: Dict[str, Any],
    requested_group_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Verify that the supplied evidence documents match a linkage record.

    A GEO subseries is accepted when it is explicitly listed under
    `geo_subseries` for a registry record.
    """

    dataset_ids = {
        str(doc.get("dataset_id", "")).upper()
        for doc in documents
        if doc.get("dataset_id")
    }

    for link in registry.get("links", []):

        if requested_group_id:
            registry_group_id = str(
                link.get("group_id", "")
            )

            if registry_group_id != requested_group_id:
                continue

        immport = str(
            link.get("immport_accession", "")
        ).upper()

        geo = str(
            link.get("geo_accession", "")
        ).upper()

        geo_subseries = {
            str(value).upper()
            for value in link.get("geo_subseries", [])
        }

        immport_present = (
            bool(immport)
            and immport in dataset_ids
        )

        geo_present = (
            (bool(geo) and geo in dataset_ids)
            or bool(geo_subseries.intersection(dataset_ids))
        )

        if immport_present and geo_present:
            return {
                "verified": True,
                "link": link,
                "matched_dataset_ids": sorted(dataset_ids),
            }

    return {
        "verified": False,
        "link": None,
        "matched_dataset_ids": sorted(dataset_ids),
    }


# ---------------------------------------------------------------------
# Criterion fusion
# ---------------------------------------------------------------------

def fuse_criterion(
    criterion: str,
    evidence_items: List[Dict[str, Any]],
) -> Dict[str, Any]:

    true_items = []
    false_items = []
    unknown_items = []

    for item in evidence_items:

        value = normalize_value(
            item.get("value", "unknown")
        )

        if value == "true":
            true_items.append(item)

        elif value == "false":
            false_items.append(item)

        else:
            unknown_items.append(item)

    conflict = bool(true_items and false_items)

    if conflict:
        final_value = "unknown"
        final_confidence = "low"

        decision_note = (
            "Conflicting evidence was found across sources; "
            "manual resolution is required."
        )

    elif false_items:
        final_value = "false"
        final_confidence = highest_confidence(false_items)

        decision_note = (
            "At least one source explicitly contradicts this criterion "
            "and no source positively establishes it."
        )

    elif true_items:
        final_value = "true"
        final_confidence = highest_confidence(true_items)

        decision_note = (
            "At least one source positively establishes this criterion "
            "and no source contradicts it."
        )

    else:
        final_value = "unknown"
        final_confidence = highest_confidence(unknown_items)

        decision_note = (
            "Available sources do not yet establish this criterion."
        )

    evidence_strings = []
    sources = []
    source_values = []

    for item in evidence_items:

        item_sources = item.get("sources", [])

        if isinstance(item_sources, str):
            item_sources = [item_sources]

        sources.extend(item_sources)

        evidence = item.get("evidence")

        if evidence:
            evidence_strings.append(str(evidence))

        source_values.append(
            {
                "value": normalize_value(
                    item.get("value", "unknown")
                ),
                "confidence": item.get(
                    "confidence",
                    "low",
                ),
                "sources": item_sources,
            }
        )

    return {
        "value": final_value,
        "confidence": final_confidence,
        "evidence": " | ".join(
            unique_list(evidence_strings)
        ),
        "sources": unique_list(sources),
        "conflict": conflict,
        "fusion_note": decision_note,
        "source_values": source_values,
    }


# ---------------------------------------------------------------------
# Assay-context fusion
# ---------------------------------------------------------------------

def fuse_assay_context(
    documents: List[Dict[str, Any]]
) -> Dict[str, Any]:

    assays: Dict[str, Dict[str, Any]] = {}
    response_modalities: Dict[str, bool] = {}
    timepoints = []

    for document in documents:

        context = document.get(
            "assay_context",
            {},
        )

        source_assays = context.get(
            "assays_detected",
            {},
        )

        for assay_name, assay_info in source_assays.items():

            if assay_name not in assays:
                assays[assay_name] = {
                    "detected": False,
                    "matched_terms": [],
                }

            if assay_info.get("detected", False):
                assays[assay_name]["detected"] = True

            assays[assay_name]["matched_terms"].extend(
                assay_info.get(
                    "matched_terms",
                    [],
                )
            )

        modalities = context.get(
            "response_modalities",
            {},
        )

        for modality, detected in modalities.items():

            response_modalities[modality] = (
                response_modalities.get(
                    modality,
                    False,
                )
                or bool(detected)
            )

        timepoints.extend(
            context.get(
                "timepoints_detected",
                [],
            )
        )

    for assay_info in assays.values():

        assay_info["matched_terms"] = unique_list(
            assay_info["matched_terms"]
        )

    numeric_timepoints = []

    for value in timepoints:

        try:
            numeric_timepoints.append(int(value))
        except (TypeError, ValueError):
            pass

    return {
        "assays_detected": assays,
        "response_modalities": response_modalities,
        "timepoints_detected": sorted(
            set(numeric_timepoints)
        ),
    }


# ---------------------------------------------------------------------
# Document fusion
# ---------------------------------------------------------------------

def fuse_documents(
    documents: List[Dict[str, Any]],
    group_id: str,
    linkage_record: Dict[str, Any],
) -> Dict[str, Any]:

    all_criteria = set()
    linked_datasets = []

    for document in documents:

        criteria = document.get(
            "criteria",
            {},
        )

        all_criteria.update(criteria.keys())

        dataset_id = document.get(
            "dataset_id",
            "UNKNOWN",
        )

        linked_datasets.append(
            {
                "dataset_id": dataset_id,
                "repository": infer_repository(
                    dataset_id,
                    document.get("repository"),
                ),
            }
        )

    fused_criteria = {}

    for criterion in sorted(all_criteria):

        criterion_items = []

        for document in documents:

            item = (
                document
                .get("criteria", {})
                .get(criterion)
            )

            if item:
                criterion_items.append(item)

        fused_criteria[criterion] = fuse_criterion(
            criterion,
            criterion_items,
        )

    conflicts = [
        criterion
        for criterion, result
        in fused_criteria.items()
        if result.get("conflict")
    ]

    return {
        "dataset_id": group_id,
        "repository": "fused",
        "linked_datasets": linked_datasets,

        "cross_repository_linkage_verified": True,
        "linkage_record": linkage_record,

        "criteria": fused_criteria,

        "assay_context": fuse_assay_context(
            documents
        ),

        "fusion_summary": {
            "n_sources": len(documents),
            "n_criteria": len(fused_criteria),
            "conflicting_criteria": conflicts,
        },

        "fusion_notes": [
            (
                "Evidence was fused only after the supplied accessions "
                "matched a verified cross-repository linkage record."
            ),
            (
                "A GEO subseries may be used when that subseries is "
                "explicitly listed under the linked GEO SuperSeries."
            ),
            (
                "Positive evidence can resolve an unknown criterion, "
                "but explicit contradictory evidence creates a conflict "
                "that remains unresolved."
            ),
            (
                "Cross-repository evidence fusion does not itself prove "
                "participant-level predictor/outcome linkage."
            ),
        ],
    }


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Fuse normalized GEO/ImmPort evidence for "
            "Hypothesis2Omics using a verified linkage registry."
        )
    )

    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        type=Path,
        help=(
            "Two or more normalized evidence JSON files"
        ),
    )

    parser.add_argument(
        "--group-id",
        required=True,
        help=(
            "Group identifier. This must match a group_id "
            "in the linkage registry."
        ),
    )

    parser.add_argument(
        "--link-registry",
        required=True,
        type=Path,
        help=(
            "JSON registry of verified cross-repository study links"
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional fused evidence JSON output",
    )

    args = parser.parse_args()

    if len(args.inputs) < 2:
        parser.error(
            "Evidence fusion requires at least two inputs."
        )

    documents = [
        load_json(path)
        for path in args.inputs
    ]

    try:
        registry = load_link_registry(
            args.link_registry
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        parser.error(
            f"Could not load linkage registry: {exc}"
        )

    link_check = linkage_is_allowed(
        documents,
        registry,
        requested_group_id=args.group_id,
    )

    if not link_check["verified"]:

        supplied = ", ".join(
            link_check["matched_dataset_ids"]
        )

        parser.error(
            "Refusing to fuse records because the supplied accessions "
            f"({supplied}) do not match group_id '{args.group_id}' "
            "in the verified linkage registry."
        )

    result = fuse_documents(
        documents,
        args.group_id,
        link_check["link"],
    )

    text = json.dumps(
        result,
        indent=2,
    )

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