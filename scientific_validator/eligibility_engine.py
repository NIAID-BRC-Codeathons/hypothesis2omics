#!/usr/bin/env python3
"""Deterministic, config-driven dataset eligibility engine for Hypothesis2Omics.

The engine is deliberately domain-agnostic. A test specification supplies the
scientific and operational criteria; a dataset-evidence record supplies
true/false/unknown values plus provenance. The engine returns an auditable
eligibility decision suitable for orchestration by an LLM/MCP layer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable

ALLOWED_VALUES = {"true", "false", "unknown"}


def _criterion_ids(spec: Dict[str, Any], section: str) -> list[str]:
    items = spec.get("eligibility", {}).get(section, [])
    out: list[str] = []
    for item in items:
        if isinstance(item, str):
            out.append(item)
        elif isinstance(item, dict) and item.get("id"):
            out.append(str(item["id"]))
        else:
            raise ValueError(f"Invalid criterion in {section}: {item!r}")
    return out


def _value(evidence: Dict[str, Any], criterion_id: str) -> str:
    item = evidence.get(criterion_id, {})
    raw = item.get("value", "unknown") if isinstance(item, dict) else item
    value = str(raw).lower()
    if value not in ALLOWED_VALUES:
        raise ValueError(
            f"{criterion_id}: invalid value {value!r}; expected one of "
            f"{sorted(ALLOWED_VALUES)}"
        )
    return value


def _details(record: Dict[str, Any], criterion_ids: Iterable[str]) -> list[Dict[str, Any]]:
    evidence = record.get("criteria", {})
    details = []
    for cid in criterion_ids:
        item = evidence.get(cid, {})
        if not isinstance(item, dict):
            item = {"value": item}
        details.append(
            {
                "criterion": cid,
                "value": _value(evidence, cid),
                "confidence": item.get("confidence"),
                "evidence": item.get("evidence"),
                "sources": item.get("sources", []),
            }
        )
    return details


def assess_dataset(test_spec: Dict[str, Any], record: Dict[str, Any]) -> Dict[str, Any]:
    scientific_ids = _criterion_ids(test_spec, "scientific_criteria")
    operational_ids = _criterion_ids(test_spec, "operational_criteria")
    if not scientific_ids:
        raise ValueError("test specification has no scientific eligibility criteria")

    evidence = record.get("criteria", {})
    scientific = {cid: _value(evidence, cid) for cid in scientific_ids}
    operational = {cid: _value(evidence, cid) for cid in operational_ids}

    false_scientific = [cid for cid, value in scientific.items() if value == "false"]
    unknown_scientific = [cid for cid, value in scientific.items() if value == "unknown"]
    false_operational = [cid for cid, value in operational.items() if value == "false"]
    unknown_operational = [cid for cid, value in operational.items() if value == "unknown"]

    if false_scientific:
        scientific_status = "excluded"
    elif unknown_scientific:
        scientific_status = "uncertain"
    else:
        scientific_status = "eligible"

    if scientific_status == "excluded":
        execution_status = "not_applicable"
    elif scientific_status == "uncertain":
        execution_status = "needs_scientific_resolution"
    elif false_operational:
        execution_status = "blocked"
    elif unknown_operational:
        execution_status = "needs_data_check"
    else:
        execution_status = "ready"

    unresolved = unknown_scientific + unknown_operational
    failed = false_scientific

    if scientific_status == "excluded":
        next_action = "Exclude from this hypothesis test; failed scientific criteria: " + ", ".join(failed)
    elif scientific_status == "uncertain":
        next_action = "Resolve scientific evidence before analysis: " + ", ".join(unknown_scientific)
    elif execution_status == "blocked":
        next_action = "Scientific match is adequate, but required data are inaccessible: " + ", ".join(false_operational)
    elif execution_status == "needs_data_check":
        next_action = "Scientific match is adequate; verify/access required data: " + ", ".join(unknown_operational)
    else:
        next_action = "Proceed to harmonization and statistical workflow specification."

    return {
        "hypothesis_id": test_spec.get("hypothesis_id"),
        "dataset_id": record.get("dataset_id", "UNKNOWN"),
        "scientific_status": scientific_status,
        "execution_status": execution_status,
        "failed_scientific_criteria": failed,
        "unresolved_criteria": unresolved,
        "next_action": next_action,
        "criterion_assessments": _details(record, scientific_ids + operational_ids),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Assess whether a dataset can test a Hypothesis2Omics test specification.")
    parser.add_argument("--spec", required=True, type=Path, help="test specification JSON")
    parser.add_argument("--dataset", required=True, type=Path, help="dataset evidence JSON")
    parser.add_argument("--output", type=Path, default=None, help="optional output JSON")
    args = parser.parse_args()

    spec = json.loads(args.spec.read_text())
    record = json.loads(args.dataset.read_text())
    result = assess_dataset(spec, record)
    text = json.dumps(result, indent=2)
    print(text)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")


if __name__ == "__main__":
    main()
