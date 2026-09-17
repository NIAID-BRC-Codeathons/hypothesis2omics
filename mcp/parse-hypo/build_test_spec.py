#!/usr/bin/env python3
"""Steps 2-3: parsed components -> test specification -> ImmPort search spec.

    01_parsed.yaml        written by parse_hypothesis.py, reviewed and edited by hand
            |
            |  STEP 2  (LLM)  expand each component into synonyms, measurements
            v                 and checkable eligibility criteria
    02_test_spec.json     consumed by scientific_validator/eligibility_engine.py
            |
            |  STEP 3  (deterministic)  collect per-component search_terms
            v
    03_search_spec.yaml   consumed by server.py:search_spec

Running this script IS the human gate: step 1 stops at 01_parsed.yaml, and nothing
proceeds until someone points this at that file.

Search terms are never invented in step 3; they are lifted from the component that
owns them in step 2, so every term traces back to a named part of the specification.

Install:
    pip install openai httpx pyyaml mcp

Run:
    ARGO_USER=ac.yourname python3 build_test_spec.py runs/HYP001/01_parsed.yaml \
        --audit --search

Requires a connection to the Argonne-auth network.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

from h2o_common import (
    DEFAULT_MODEL,
    EmptyCompletion,
    _ask,
    _header,
    _provenance,
    _scalar,
    client as argo_client,
)
from parse_hypothesis import load_parsed

sys.path.insert(0, str(Path(__file__).resolve().parent))
from server import FACETS, list_facet_values, search_spec  # noqa: E402

GENERATOR = "build_test_spec.py"

# The four search groups of Example_input_output.md, in document order.
CANONICAL_GROUPS = ("intervention", "predictor", "outcome", "assay")

GROUNDED_FACETS = ("species", "assayMethod", "researchFocus")
MAX_VOCAB_VALUES = 60

SPEC_SCHEMA = """{
  "hypothesis_id": "HYP001",
  "hypothesis_text": "the hypothesis, verbatim",
  "analysis_intent": "association | causal | difference",
  "population": "who must have been studied, one phrase",

  "intervention": {
    "name": "canonical name",
    "search_terms": ["every synonym, abbreviation, expansion and spelling"]
  },
  "predictor": {
    "biological_factor": "canonical name",
    "possible_measurements": ["how it could be quantified in a dataset"],
    "search_terms": ["gene symbol, aliases, pathway names"]
  },
  "outcome": {
    "name": "canonical name",
    "possible_measurements": ["how it could be quantified in a dataset"],
    "search_terms": ["synonyms for the outcome concept"]
  },
  "assay": {
    "search_terms": ["measurement technologies that would produce the required data"]
  },

  "required_data": ["the data types a dataset must contain"],
  "preferred_predictor_timepoint": ["pre-vaccination", "early post-vaccination"],
  "preferred_outcome_timepoint": ["post-vaccination"],
  "preferred_sample_types": ["PBMC", "whole blood"],
  "required_features": ["short phrase", "short phrase"],

  "eligibility": {
    "scientific_criteria": [
      {"id": "snake_case_id",
       "description": "one sentence, checkable against study metadata",
       "exclusion": "the short phrase for failing it, e.g. 'non-human study'"}
    ],
    "operational_criteria": [
      {"id": "snake_case_id", "description": "whether the data can actually be retrieved"}
    ]
  },

  "primary_statistical_test": {"type": "association/regression"},
  "recommended_analysis": {
    "predictor": "the concrete measurement to use",
    "outcome": "the concrete measurement to use",
    "model": "a model formula, e.g. outcome ~ predictor + covariate"
  },

  "search": {
    "require": ["intervention"],
    "min_groups": 1,
    "filters": {"species": "Homo sapiens"},
    "limit": 100,
    "per_term": 200
  }
}"""

SPEC_SYSTEM = """You expand the parsed components of a biological hypothesis into a test \
specification: the synonyms needed to find candidate datasets in ImmPort, a human \
immunology study repository, and the criteria needed to decide whether a dataset can \
test the hypothesis.

You are a middle step of an automated pipeline. Downstream code queries ImmPort with \
your search_terms and applies your eligibility criteria mechanically. Your output is \
consumed by a program, never read by a person first.

EXPANDING search_terms -- these rules come from measurements against the live API:
* A query ANDs every word in it. "yellow fever CD8 transcriptomics" matches almost \
nothing. Keep each term to one concept, usually one to three words.
* There is no OR operator. The pipeline runs one query per term and unions the \
results, so terms within a component are effectively ORed. Be generous: every synonym, \
abbreviation, expansion and spelling a submitter might plausibly have used. Gene symbol \
and alias both. Hyphenated and unhyphenated both. Trade names and strain names.
* Never include a fragment of a compound identifier as its own term. Measured \
2026-09-16: for YF-17D the term "17D" returned 14 studies of which 9 were unrelated \
(influenza, HIV, kidney transplant), because short fragments collide with arm codes and \
protocol identifiers. The full forms "YF-17D", "YF17D", "YFV", "Stamaril" were clean. \
Spell out each synonym in full; never decompose one.
* ImmPort indexes study titles, conditions and research focus -- NOT molecular content. \
Measured 2026-09-16: "GCN2", "EIF2AK4", "ATF4" and every other gene symbol returned zero \
studies. Still populate predictor.search_terms: they cost nothing, they document intent, \
and they contribute to ranking. But NEVER put "predictor" in search.require -- that \
returns nothing at all. Whether a gene is measurable is decided later from the \
expression matrix, not from search.
* search.require names the components a study MUST match, and decides what comes back \
at all; every other component only affects ranking. A loose term in a required component \
costs precision directly; a loose term elsewhere is nearly free. Normally require only \
"intervention". Do not require "outcome" or "assay": a study can satisfy them without \
saying so, and requiring them silently drops good datasets.
* search.filters are ANDed with every query and MUST use the controlled vocabulary given \
below, verbatim. Prefer few filters: each one can only lose studies.
* required_features is a human-readable summary. Short noun phrases, at most six words \
each ("human subjects", "transcriptomics", "YF-17D vaccination"). It does not affect \
the query.

WRITING eligibility criteria:
* Each criterion is answered per dataset as only true, false, or unknown. Write criteria \
a reader of study metadata could answer that way. "Participants received YF-17D \
vaccination" is checkable; "the study is high quality" is not.
* scientific_criteria decide whether the dataset can test the hypothesis at all: species, \
exposure, required measurements, timepoints, and whether predictor and outcome can be \
linked within a participant. That last one is easy to forget and is usually the criterion \
that excludes a dataset.
* operational_criteria decide only whether the data can currently be retrieved.
* Any single unknown forces the dataset to "uncertain", so do not pad the list. Five to \
eight scientific criteria is typical.
* Every scientific criterion also needs an "exclusion": the short phrase naming the \
dataset that fails it, as it would appear in an exclusion list. "Study includes human \
participants" excludes a "non-human study". Six words at most. Write it; do not \
mechanically negate the description.

DO NOT emit "inclusion_criteria", "exclusion_criteria", or \
"primary_statistical_test.example_model". The pipeline derives all three from the fields \
above so the two renderings cannot drift apart.

Return one JSON object matching this schema, and nothing else -- no prose, no code \
fences:

""" + SPEC_SCHEMA


# ------------------------------------------------------------------ vocabulary


def _vocabulary() -> dict[str, list[str]]:
    """Controlled-vocabulary values for the facets we let the model filter on."""
    try:
        every = list_facet_values()
    except Exception as exc:  # network, ImmPort outage
        print(f"warning: could not fetch ImmPort facet vocabulary ({exc}); "
              f"proceeding without filter grounding", file=sys.stderr)
        return {}
    return {f: every[f][:MAX_VOCAB_VALUES] for f in GROUNDED_FACETS if every.get(f)}


def _vocabulary_block(vocab: dict[str, list[str]]) -> str:
    if not vocab:
        return "CONTROLLED VOCABULARY: unavailable. Omit `filters` entirely."
    lines = ["CONTROLLED VOCABULARY for search.filters -- use these strings verbatim:"]
    for facet, values in vocab.items():
        lines.append(f"  {facet}: {', '.join(values)}")
    lines.append(f"Valid facet names: {', '.join(sorted(FACETS))}")
    return "\n".join(lines)


# ------------------------------------------------------------------ validation


def validate_spec(obj: dict[str, Any], vocab: dict[str, list[str]]) -> list[str]:
    problems: list[str] = []

    # The four search groups sit at the top level, matching Step 2 of
    # Example_input_output.md. A group may be absent -- an observational
    # hypothesis has no intervention -- but one that is present must be usable.
    components = {g: obj[g] for g in CANONICAL_GROUPS if g in obj}
    if not components:
        return [f"at least one of {list(CANONICAL_GROUPS)} must be present"]

    for name, body in components.items():
        if not isinstance(body, dict):
            problems.append(f"{name} must be an object")
            continue
        terms = body.get("search_terms")
        if not isinstance(terms, list) or not terms:
            problems.append(f"{name}.search_terms must be a non-empty list")
            continue
        for term in terms:
            if not isinstance(term, str) or not term.strip():
                problems.append(f"{name}.search_terms has a bad term: {term!r}")

    search = obj.get("search") or {}
    if not isinstance(search, dict):
        problems.append("search must be an object")
        search = {}

    require = search.get("require") or []
    if not isinstance(require, list):
        problems.append("search.require must be a list")
        require = []
    for name in require:
        if name not in components:
            problems.append(f"search.require names unknown component {name!r}")
    if "predictor" in require:
        problems.append("search.require must not contain 'predictor': ImmPort does not "
                        "index gene symbols, so requiring it returns zero studies")
    if not require:
        problems.append("search.require must name at least one component, normally 'intervention'")

    filters = search.get("filters") or {}
    if not isinstance(filters, dict):
        problems.append("search.filters must be an object")
        filters = {}
    for facet, value in filters.items():
        if facet not in FACETS:
            problems.append(f"unknown facet {facet!r}; valid facets: {sorted(FACETS)}")
            continue
        allowed = vocab.get(facet)
        if not allowed:
            continue
        for v in value if isinstance(value, (list, tuple)) else [value]:
            if str(v) not in allowed:
                problems.append(f"search.filters.{facet}={v!r} is not in the controlled "
                                f"vocabulary; valid values include: {', '.join(allowed[:8])}")

    for field in ("limit", "per_term", "min_groups"):
        if field in search and not isinstance(search[field], int):
            problems.append(f"search.{field} must be an integer, got {search[field]!r}")

    for field in ("required_features", "required_data", "preferred_sample_types",
                  "preferred_predictor_timepoint", "preferred_outcome_timepoint"):
        value = obj.get(field)
        if value is not None and (
            not isinstance(value, list) or any(not isinstance(v, str) for v in value)
        ):
            problems.append(f"{field} must be a list of strings")

    eligibility = obj.get("eligibility")
    if not isinstance(eligibility, dict):
        problems.append("eligibility must be an object")
    else:
        scientific = eligibility.get("scientific_criteria")
        if not isinstance(scientific, list) or not scientific:
            problems.append("eligibility.scientific_criteria must be a non-empty list")
        seen: set[str] = set()
        for section in ("scientific_criteria", "operational_criteria"):
            for item in eligibility.get(section) or []:
                if not isinstance(item, dict) or not item.get("id"):
                    problems.append(f"{section}: every criterion needs an id: {item!r}")
                    continue
                cid = str(item["id"])
                if not re.fullmatch(r"[a-z0-9_]+", cid):
                    problems.append(f"criterion id must be snake_case: {cid!r}")
                if cid in seen:
                    problems.append(f"duplicate criterion id: {cid!r}")
                seen.add(cid)
                if section != "scientific_criteria":
                    continue
                # inclusion_criteria and exclusion_criteria are derived from these
                # two fields, so a missing one silently shortens a rendered list.
                if not str(item.get("description") or "").strip():
                    problems.append(f"scientific criterion {cid!r} needs a description")
                if not str(item.get("exclusion") or "").strip():
                    problems.append(f"scientific criterion {cid!r} needs an 'exclusion' phrase")

    return problems


# ------------------------------------------------------------------- derived


def add_derived_fields(spec: dict[str, Any]) -> dict[str, Any]:
    """Fill the Example_input_output.md keys that restate other fields.

    Step 2 of the document lists `inclusion_criteria`, `exclusion_criteria` and a
    statistical test; the eligibility engine reads `eligibility.*_criteria` and
    `recommended_analysis`. Both shapes have to be in the file, so the document's
    are computed from the engine's here rather than asked for twice -- a model
    writing the same content in two places is a drift waiting to happen.

    Mutates and returns `spec`.
    """
    scientific = (spec.get("eligibility") or {}).get("scientific_criteria") or []
    inclusion, exclusion = [], []
    for item in scientific:
        if not isinstance(item, dict):
            continue
        # validate_spec has already required both, so a missing one here would be
        # a bug in this function's caller rather than bad model output.
        inclusion.append(str(item.get("description", "")).strip())
        exclusion.append(str(item.get("exclusion", "")).strip())
    spec["inclusion_criteria"] = inclusion
    spec["exclusion_criteria"] = exclusion

    model = (spec.get("recommended_analysis") or {}).get("model")
    test = spec.get("primary_statistical_test")
    if not isinstance(test, dict):
        test = {} if test is None else {"type": str(test)}
    if model:
        test["example_model"] = model
    spec["primary_statistical_test"] = test
    return spec


# --------------------------------------------------------------------- step 3


def to_search_spec(spec: dict[str, Any]) -> dict[str, Any]:
    """Step 3, deterministic: collect each component's search_terms into groups."""
    search = spec.get("search") or {}
    groups = {
        name: list(dict.fromkeys(spec[name]["search_terms"]))
        for name in CANONICAL_GROUPS
        if isinstance(spec.get(name), dict) and spec[name].get("search_terms")
    }
    out: dict[str, Any] = {"repositories": ["ImmPort"], "search_terms": groups}
    if spec.get("required_features"):
        out["required_features"] = spec["required_features"]
    out["require"] = [g for g in search.get("require") or [] if g in groups]
    out["min_groups"] = int(search.get("min_groups", 1))
    if search.get("filters"):
        out["filters"] = search["filters"]
    out["limit"] = int(search.get("limit", 100))
    out["per_term"] = int(search.get("per_term", 200))
    return out


def search_spec_to_yaml(spec: dict[str, Any], poc_format: bool = False,
                        provenance: dict | None = None) -> str:
    """Render the Step 3 document of Example_input_output.md."""
    out = _header(provenance)
    out += ["repositories:", "  - ImmPort", "", "search_terms:"]

    groups = spec.get("search_terms") or {}
    for i, group in enumerate(groups):
        if i:
            out.append("")
        out.append(f"  {group}:")
        out += [f"    - {_scalar(t)}" for t in groups[group]]

    if spec.get("required_features"):
        out += ["", "required_features:"]
        out += [f"  - {_scalar(f)}" for f in spec["required_features"]]

    if poc_format:
        return "\n".join(out) + "\n"

    out += [
        "",
        "# Below this line is what actually constrains the query. `repositories` and",
        "# `required_features` above are documentation; search_spec ignores them.",
    ]
    if spec.get("require"):
        out.append("require:")
        out += [f"  - {_scalar(g)}" for g in spec["require"]]
    out.append(f"min_groups: {int(spec.get('min_groups', 1))}")
    if spec.get("filters"):
        out.append("filters:")
        for facet, value in spec["filters"].items():
            if isinstance(value, (list, tuple)):
                out.append(f"  {facet}:")
                out += [f"    - {_scalar(v)}" for v in value]
            else:
                out.append(f"  {facet}: {_scalar(value)}")
    out.append(f"limit: {int(spec.get('limit', 100))}")
    out.append(f"per_term: {int(spec.get('per_term', 200))}")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------- audit


def audit_terms(spec: dict[str, Any]) -> list[dict[str, Any]]:
    """Per-term hit counts, so a loose term cannot hide inside a union.

    `hits` is what the term returns on its own; `unique` is what only this term
    found. A term with unique=0 is redundant; a required-group term with high
    `hits` and high `unique` is where false positives enter.
    """
    import httpx

    from server import _filters, _hits

    filters = _filters(spec.get("filters") or {})
    per_term = max(1, min(int(spec.get("per_term", 200)), 1000))
    required = set(spec.get("require") or [])

    found: dict[tuple[str, str], set[str]] = {}
    with httpx.Client() as client:
        for group, terms in (spec.get("search_terms") or {}).items():
            for term in dict.fromkeys(terms or []):
                hits = _hits(client, str(term), filters, per_term, "study_accession")
                found[(group, str(term))] = {h["_id"] for h in hits}

    rows = []
    for (group, term), ids in found.items():
        others: set[str] = set()
        for key, other in found.items():
            if key != (group, term):
                others |= other
        rows.append({
            "group": group,
            "term": term,
            "required_group": group in required,
            "hits": len(ids),
            "unique": len(ids - others),
        })
    rows.sort(key=lambda r: (not r["required_group"], -r["unique"], -r["hits"]))
    return rows


def _print_audit(rows: list[dict[str, Any]]) -> None:
    print(f"\n{'group':<14} {'term':<32} {'req':<4} {'hits':>5} {'uniq':>5}", file=sys.stderr)
    for r in rows:
        flag = ""
        if r["hits"] == 0:
            flag = "  <- no hits"
        elif r["required_group"] and r["unique"] >= 3:
            flag = "  <- check: widens the required set"
        print(f"{r['group']:<14} {r['term'][:32]:<32} {'yes' if r['required_group'] else '':<4} "
              f"{r['hits']:>5} {r['unique']:>5}{flag}", file=sys.stderr)


# -------------------------------------------------------------------- pipeline


def build_test_spec(client, hypothesis: str, parsed: dict, model: str,
                    max_tokens: int) -> tuple[dict, dict]:
    """Step 2: components -> synonyms, measurements and eligibility criteria."""
    vocab = _vocabulary()
    system = SPEC_SYSTEM + "\n\n" + _vocabulary_block(vocab)
    user = (
        f"Hypothesis:\n{hypothesis.strip()}\n\n"
        f"Parsed components:\n{json.dumps(parsed, indent=2)}"
    )
    spec, meta = _ask(
        client, model, system, user, lambda o: validate_spec(o, vocab),
        max_tokens, "step 2 (test spec)",
    )
    spec.setdefault("hypothesis_text", hypothesis.strip())
    add_derived_fields(spec)
    prov = _provenance(GENERATOR, model, system, "2-test-spec", meta)
    prov["grounded_facets"] = ",".join(sorted(vocab)) or "none"
    return spec, prov


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Expand a parsed hypothesis into a test spec and an ImmPort search spec."
    )
    parser.add_argument("parsed", type=Path,
                        help="01_parsed.yaml from parse_hypothesis.py, hand-edits and all")
    parser.add_argument("--outdir", type=Path,
                        help="where to write 02/03 (default: alongside the parsed file)")
    parser.add_argument("--model", default=DEFAULT_MODEL,
                        help=f"Argo model (default: {DEFAULT_MODEL})")
    parser.add_argument("--hypothesis-id", help="override the generated hypothesis_id")
    parser.add_argument("--stop-after-spec", action="store_true",
                        help="stop after 02_test_spec.json")
    parser.add_argument("--poc-format", action="store_true",
                        help="emit only the Example_input_output.md keys, omitting require/filters")
    parser.add_argument("--no-provenance", action="store_true", help="omit the comment headers")
    parser.add_argument("--audit", action="store_true",
                        help="report per-term hit and unique counts")
    parser.add_argument("--search", action="store_true", help="also run the search against ImmPort")
    parser.add_argument("--max-tokens", type=int, default=4000)
    args = parser.parse_args()

    outdir = args.outdir or args.parsed.parent
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        parsed, hypothesis = load_parsed(args.parsed)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    # ---- step 2
    try:
        spec, prov = build_test_spec(argo_client(), hypothesis, parsed,
                                     args.model, args.max_tokens)
    except EmptyCompletion as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.hypothesis_id:
        spec["hypothesis_id"] = args.hypothesis_id
    spec["provenance"] = prov
    path = outdir / "02_test_spec.json"
    path.write_text(json.dumps(spec, indent=2) + "\n")
    print(f"wrote {path}", file=sys.stderr)

    if args.stop_after_spec:
        return 0

    # ---- step 3, deterministic
    search = to_search_spec(spec)
    path = outdir / "03_search_spec.yaml"
    header = None if args.no_provenance else {
        "step": "3-search-spec",
        "derived_from": "02_test_spec.json",
        "method": "deterministic collection of components[*].search_terms",
    }
    path.write_text(search_spec_to_yaml(search, args.poc_format, header))
    print(f"wrote {path}", file=sys.stderr)

    if args.audit:
        _print_audit(audit_terms(search))

    if args.search:
        accessions = search_spec(search)
        print(f"\n{len(accessions)} studies: {accessions}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
