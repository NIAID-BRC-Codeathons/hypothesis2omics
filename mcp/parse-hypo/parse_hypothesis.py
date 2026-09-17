#!/usr/bin/env python3
"""Step 1: a plain-English hypothesis -> its named components.

    plain text hypothesis
            |
            |  STEP 1  (LLM)  parse into components
            v
    01_parsed.yaml        hypothesis_text + intervention / predictor / cell_type /
                          outcome / relationship_type / directionality + warnings

That file is the whole interface to step 2, and the human gate sits on it: read it,
edit anything the parser got wrong, then run build_test_spec.py on it.

This module imports nothing from server.py and never touches ImmPort, so it can be
run and tested without the repository being reachable. Only Argo is required.

Install:
    pip install openai pyyaml

Run:
    OPENAI_API_KEY=... python3 parse_hypothesis.py \
        --hypothesis "GCN2/EIF2AK4 activity is associated with the magnitude of the
                      CD8+ T-cell response following YF-17D vaccination." \
        --outdir runs/HYP001

With the default gateway this requires a connection to the Argonne-auth network.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Any

import yaml

from h2o_common import (
    DEFAULT_MODEL,
    EmptyCompletion,
    _ask,
    _header,
    _provenance,
    _scalar,
    client as llm_client,
)

GENERATOR = "parse_hypothesis.py"

PARSED_FIELDS = (
    "intervention",
    "predictor",
    "cell_type",
    "outcome",
    "relationship_type",
    "directionality",
)

PARSE_SYSTEM = """You decompose a biological hypothesis into its components. You are the \
first step of an automated pipeline; a later step expands what you produce, so be \
faithful and literal rather than creative.

Return one JSON object, nothing else -- no prose, no code fences:

{
  "intervention": "the exposure, treatment or vaccine, or null if the hypothesis is observational",
  "predictor": "the biological factor proposed to explain variation in the outcome",
  "cell_type": "the cell type the outcome is measured in, or null",
  "outcome": "what varies -- the quantity being explained",
  "relationship_type": "association | causal | difference",
  "directionality": "positive | negative | non-directional",
  "warnings": ["short notes on anything you inferred rather than read"]
}

Rules:
* Use the hypothesis's own wording where you can. Do not substitute a broader concept.
* directionality is non-directional unless the hypothesis states a direction \
("higher X is associated with greater Y" is positive; "is associated with" alone is \
non-directional).
* relationship_type is causal only if the hypothesis claims causation, not correlation.
* Null is a valid answer. Do not invent an intervention for an observational hypothesis.

Warnings are read by the scientist who reviews your output before the pipeline \
continues, so they must point at real problems:
* Add one whenever you fill a field the hypothesis does not actually state, and name \
that field: "cell_type inferred; not stated in the hypothesis".
* Add one when the outcome names no measurable quantity ("matters for", "is important \
in"), when the predictor is too broad to measure, or when the hypothesis names no \
population.
* Keep each under twelve words. Return an empty list when the hypothesis states \
everything you used -- do not invent concerns to look thorough."""


# ------------------------------------------------------------------ validation


def validate_parsed(obj: dict[str, Any]) -> list[str]:
    problems = []
    for field in PARSED_FIELDS:
        if field not in obj:
            problems.append(f"missing field: {field}")
    if obj.get("relationship_type") not in (None, "association", "causal", "difference"):
        problems.append(f"relationship_type must be association|causal|difference, "
                        f"got {obj.get('relationship_type')!r}")
    if obj.get("directionality") not in (None, "positive", "negative", "non-directional"):
        problems.append(f"directionality must be positive|negative|non-directional, "
                        f"got {obj.get('directionality')!r}")
    if not obj.get("outcome"):
        problems.append("outcome is required; a hypothesis must have something that varies")
    if not obj.get("predictor"):
        problems.append("predictor is required")
    warnings = obj.get("warnings")
    if warnings is not None and (
        not isinstance(warnings, list) or any(not isinstance(w, str) for w in warnings)
    ):
        problems.append("warnings must be a list of strings")
    return problems


# ----------------------------------------------------------------------- YAML


def parsed_to_yaml(parsed: dict[str, Any], hypothesis: str,
                   provenance: dict | None = None) -> str:
    out = _header(provenance)
    # Carried so step 2 can work from the original wording instead of
    # reconstructing it by concatenating the fields below.
    out += [f"hypothesis_text: {_scalar(hypothesis.strip())}", ""]
    out += [f"{field}: {_scalar(parsed.get(field))}" for field in PARSED_FIELDS]

    # Advisory only: what the parser inferred rather than read. Nothing downstream
    # branches on these; they exist for the reviewer at the human gate.
    warnings = parsed.get("warnings") or []
    out += ["", "warnings:"] if warnings else ["", "warnings: []"]
    out += [f"  - {_scalar(w)}" for w in warnings]
    return "\n".join(out) + "\n"


def load_parsed(path: Path) -> tuple[dict[str, Any], str]:
    """Read a 01_parsed.yaml back, hand-edits and all.

    Returns the components and the hypothesis text. Raises ValueError listing
    every problem if the file is not a valid parse -- a hand edit that breaks the
    schema should stop the pipeline here, not surface downstream.
    """
    loaded = yaml.safe_load(path.read_text()) or {}
    problems = validate_parsed(loaded)
    if problems:
        raise ValueError(f"{path} is not a valid parsed hypothesis:\n"
                         + "\n".join(f"  - {p}" for p in problems))
    parsed = {f: loaded.get(f) for f in PARSED_FIELDS}
    if loaded.get("warnings"):
        parsed["warnings"] = loaded["warnings"]
    hypothesis = loaded.get("hypothesis_text") or " ".join(
        str(loaded[f]) for f in PARSED_FIELDS if loaded.get(f)
    )
    return parsed, hypothesis


# -------------------------------------------------------------------- pipeline


def parse_hypothesis(client, hypothesis: str, model: str, max_tokens: int) -> tuple[dict, dict]:
    """Step 1: plain text -> components."""
    parsed, meta = _ask(
        client, model, PARSE_SYSTEM, f"Hypothesis:\n{hypothesis.strip()}",
        validate_parsed, max_tokens, "step 1 (parse)",
    )
    prov = _provenance(GENERATOR, model, PARSE_SYSTEM, "1-parse", meta)
    prov["hypothesis_sha256"] = hashlib.sha256(hypothesis.strip().encode()).hexdigest()[:16]
    return parsed, prov


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse a plain-text biological hypothesis into its components."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--hypothesis", help="hypothesis text")
    source.add_argument("--hypothesis-file", type=Path, help="file containing the hypothesis")
    parser.add_argument("--outdir", type=Path, default=Path("."),
                        help="directory to write 01_parsed.yaml into")
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"model id on the configured gateway (default: {DEFAULT_MODEL})")
    parser.add_argument("--no-provenance", action="store_true",
                        help="omit the comment header")
    parser.add_argument("--stdout", action="store_true",
                        help="print the YAML instead of writing a file")
    parser.add_argument("--max-tokens", type=int, default=4000)
    args = parser.parse_args()

    hypothesis = args.hypothesis or args.hypothesis_file.read_text()

    try:
        parsed, prov = parse_hypothesis(llm_client(), hypothesis, args.model, args.max_tokens)
    except EmptyCompletion as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    text = parsed_to_yaml(parsed, hypothesis, None if args.no_provenance else prov)
    if args.stdout:
        print(text, end="")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    path = args.outdir / "01_parsed.yaml"
    path.write_text(text)
    print(f"wrote {path}", file=sys.stderr)
    print(f"review it, edit anything wrong, then: "
          f"python3 build_test_spec.py {path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
