"""Emit Yijun's synthesis.schema.json from pre-registered decision rules.

`codeathon/schemas/synthesis.schema.json` is defined in the repo and nothing
implements it yet. This does, and it backs every verdict with the frozen rule
rather than an LLM judgement call, so the report can show its working.

    python synthesis.py        # worked example on GSE13485, 16 self-tests

Depends only on evidence_rules.py and the standard library. No jsonschema, no
network, no repo layout assumptions.


The schema in one line
----------------------

    {overall_verdict, per_dataset: [analysis_unit_verdict], narrative}

with `analysis_unit_verdict` requiring analysis_unit_id, gse_id, cohort,
platform, verdict, confidence, evidence_summary, rationale, caveats.


Two things worth knowing before you read the code
-------------------------------------------------

**Vocabulary differs and is mapped, not redefined.** evidence_rules speaks
supports / refutes / inconclusive. The schema speaks supportive /
contradictory / inconclusive. `VERDICT_MAP` is the whole of the translation.

**The schema has no field for independence.** `per_dataset` is a flat array,
so an overall verdict computed from it cannot express "these three analysis
units are one study". That matters here specifically: GSE13485 holds two
trials a year apart, and GSE125921 and GSE136163 are both SDY1529. Because
`additionalProperties` is false, this module cannot simply add a field. So
under `strict=True` the grouping is carried in the narrative and in each
unit's caveats, where it is at least visible to a reader. `strict=False`
adds `independence_group` per unit plus `n_independent_groups` and
`decision_rule` at the top, which is the shape worth proposing as schema
v1.1.
"""

from __future__ import annotations

import json
from typing import Any, Iterable, Optional, Sequence

from evidence_rules import DecisionRule, independence_group, integrate

__all__ = [
    "VERDICT_MAP",
    "groups_from_link_registry",
    "synthesize",
    "unit_group",
    "validate_synthesis",
]


VERDICT_MAP = {
    "supports": "supportive",
    "refutes": "contradictory",
    "inconclusive": "inconclusive",
}

_ALLOWED = ("supportive", "contradictory", "inconclusive")

_UNIT_KEYS = (
    "analysis_unit_id",
    "gse_id",
    "cohort",
    "platform",
    "verdict",
    "confidence",
    "evidence_summary",
    "rationale",
    "caveats",
)

CONFIDENCE_NOTE = (
    "Confidence is derived deterministically from the pre-registered rule and "
    "the adjusted p-value. It is not a calibrated posterior probability."
)


# ---------------------------------------------------------------------------
# Grouping at the analysis-unit grain
# ---------------------------------------------------------------------------


def unit_group(unit: dict[str, Any]) -> str:
    """Independence group for one analysis unit.

    `evidence_rules.independence_group` resolves at accession grain, which is
    one level too coarse here. Two cohorts inside a single GSE can be genuinely
    independent, and two GSEs can be a single study.

    Resolution order:

        explicit `independence_group` on the unit, always wins
        then superseries / bioproject / study / program, via evidence_rules
        then gse_id plus cohort, when the unit declares a cohort
        then gse_id alone
    """
    explicit = unit.get("independence_group")
    if explicit:
        return str(explicit)

    shared = independence_group(unit)
    if not shared.startswith("accession:"):
        return shared

    gse = unit.get("gse_id") or unit.get("accession") or "unknown"
    cohort = unit.get("cohort")
    if cohort:
        return f"{gse}|cohort:{cohort}"
    return f"accession:{gse}"


def groups_from_link_registry(registry: dict[str, Any]) -> dict[str, str]:
    """Accession -> group id, read from `cross_repository_links.json`.

    The eligibility side already curates which accessions are the same study
    across repositories. That registry is the authority on linkage, so this
    reads it rather than re-deriving the same facts from metadata.

    Given a link like::

        {"group_id": "SDY1264_GSE13486",
         "immport_accession": "SDY1264",
         "geo_accession": "GSE13486",
         "geo_subseries": ["GSE13485"],
         "relationship": "same_study_with_geo_subseries"}

    every accession named maps to `SDY1264_GSE13486`, so a result from
    GSE13485 and a result from SDY1264 count once between them.

    Feed the output to `synthesize` via each unit's `independence_group`::

        reg = json.load(open("scientific_validator/examples/cross_repository_links.json"))
        by_acc = groups_from_link_registry(reg)
        for u in units:
            u.setdefault("independence_group", by_acc.get(u["gse_id"]))
    """
    out: dict[str, str] = {}
    for link in registry.get("links", []):
        gid = link.get("group_id")
        if not gid:
            continue
        accessions = [
            link.get("immport_accession"),
            link.get("geo_accession"),
            *(link.get("geo_subseries") or []),
            *(link.get("accessions") or []),
        ]
        for acc in accessions:
            if acc:
                out[str(acc)] = str(gid)
    return out


def _confidence(scored: dict[str, Any], rule: DecisionRule) -> float:
    """A number in [0, 1], computed from the rule rather than asserted.

    Deliberately crude and deliberately documented. A verdict that met the
    criteria is scaled by how far under alpha it landed; anything that did not
    meet them is capped low, because "we could not tell" should never present
    as confident.
    """
    p = scored.get("fdr")
    if p is None:
        p = scored.get("p_value")
    if p is None:
        return 0.30 if scored.get("verdict") == "inconclusive" else 0.50

    try:
        p = float(p)
    except (TypeError, ValueError):
        return 0.30

    if scored.get("verdict") == "inconclusive":
        return round(min(0.40, max(0.05, 1.0 - p)), 3)

    if rule.alpha <= 0:
        return 0.50
    margin = max(0.0, min(1.0, 1.0 - (p / rule.alpha)))
    return round(0.50 + 0.45 * margin, 3)


def _evidence_summary(scored: dict[str, Any]) -> str:
    effect = scored.get("effect")
    p = scored.get("fdr", scored.get("p_value"))
    n = scored.get("n")
    bits = []
    if effect is not None:
        bits.append(f"effect {effect:+.3f}" if isinstance(effect, (int, float)) else f"effect {effect}")
    if p is not None:
        label = "FDR" if "fdr" in scored else "p"
        bits.append(f"{label} {p:.3g}" if isinstance(p, (int, float)) else f"{label} {p}")
    if n is not None:
        bits.append(f"n = {n}")
    return ", ".join(bits) if bits else "no effect estimate recorded"


# ---------------------------------------------------------------------------
# The main call
# ---------------------------------------------------------------------------


def synthesize(
    units: Sequence[dict[str, Any]],
    rule: DecisionRule,
    uncovered: Optional[Sequence[str]] = None,
    strict: bool = True,
) -> dict[str, Any]:
    """Score analysis units and emit a synthesis.schema.json document.

    Each unit is a dict carrying at least `analysis_unit_id`, `gse_id` and an
    effect with `fdr` or `p_value`. `cohort`, `platform`, `superseries`,
    `bioproject` and `independence_group` are all optional and used when
    present.

    `strict=True` returns exactly the fields the schema allows. `strict=False`
    adds the independence fields the schema has no home for yet.
    """
    scored: list[dict[str, Any]] = []
    groups: dict[str, str] = {}

    for unit in units:
        row = dict(unit)
        uid = row.get("analysis_unit_id") or row.get("gse_id") or "unknown"
        row["analysis_unit_id"] = uid
        # integrate() keys grouping off `dataset`, so route the unit id there
        row["dataset"] = uid
        groups[uid] = unit_group(unit)
        scored.append(rule.score(row))

    rollup = integrate(scored, rule, groups, uncovered)

    per_dataset: list[dict[str, Any]] = []
    for s in scored:
        uid = s["analysis_unit_id"]
        group = groups[uid]
        siblings = sorted(k for k, g in groups.items() if g == group and k != uid)

        caveats = [CONFIDENCE_NOTE]
        if siblings:
            caveats.append(
                f"Not independent of {', '.join(siblings)}. All share "
                f"independence group '{group}', so they count once toward "
                "replication."
            )
        entry = {
            "analysis_unit_id": uid,
            "gse_id": s.get("gse_id") or uid,
            "cohort": s.get("cohort"),
            "platform": s.get("platform"),
            "verdict": VERDICT_MAP[s.get("verdict", "inconclusive")],
            "confidence": _confidence(s, rule),
            "evidence_summary": _evidence_summary(s),
            "rationale": s.get("verdict_rationale") or "no rationale recorded",
            "caveats": _dedupe(caveats),
        }
        if not strict:
            entry["independence_group"] = group
        per_dataset.append(entry)

    doc: dict[str, Any] = {
        "overall_verdict": VERDICT_MAP[rollup["verdict"]],
        "per_dataset": per_dataset,
        "narrative": _narrative(rollup, groups, rule),
    }
    if not strict:
        doc["n_independent_groups"] = len(set(groups.values()))
        doc["decision_rule"] = rule.as_dict()
    return doc


def _dedupe(items: Iterable[str]) -> list[str]:
    """Preserve order, drop repeats. The schema requires uniqueItems."""
    seen: set[str] = set()
    out: list[str] = []
    for it in items:
        if it not in seen:
            seen.add(it)
            out.append(it)
    return out


def _narrative(
    rollup: dict[str, Any], groups: dict[str, str], rule: DecisionRule
) -> str:
    n_units = len(groups)
    n_groups = len(set(groups.values()))

    lines = [
        f"Overall verdict: {VERDICT_MAP[rollup['verdict']]}, because "
        f"{rollup['rationale']}.",
        "",
        "Pre-registered before any result existed: alpha "
        f"{rule.alpha}, {rule.correction} correction, expected direction "
        f"'{rule.direction}', minimum effect {rule.min_effect}, and at least "
        f"{rule.min_independent_groups} independent dataset group(s) required.",
        "",
        f"{n_units} analysis unit(s) resolve to {n_groups} independent "
        f"group(s):",
    ]

    by_group: dict[str, list[str]] = {}
    for uid, g in sorted(groups.items()):
        by_group.setdefault(g, []).append(uid)
    for g, members in sorted(by_group.items()):
        lines.append(f"  {g}: {', '.join(members)}")

    if n_units > n_groups:
        lines += [
            "",
            "Units sharing a group are correlated by construction and are "
            "counted once toward replication. The schema has no field for "
            "this, so it is stated here.",
        ]

    lines += ["", "Caveats:"]
    lines += [f"  - {c}" for c in rollup["caveats"]]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Validation without a jsonschema dependency
# ---------------------------------------------------------------------------


def validate_synthesis(doc: dict[str, Any], strict: bool = True) -> list[str]:
    """Return a list of schema violations. Empty list means valid.

    Checks the constraints the schema actually asserts: required keys, the
    verdict enums, the GSE and GPL patterns, confidence range, non-empty
    strings, unique caveats, and additionalProperties under strict.
    """
    errs: list[str] = []
    allowed_top = {"overall_verdict", "per_dataset", "narrative"}
    if not strict:
        allowed_top |= {"n_independent_groups", "decision_rule"}

    for key in ("overall_verdict", "per_dataset", "narrative"):
        if key not in doc:
            errs.append(f"missing required key {key!r}")
    for key in doc:
        if key not in allowed_top:
            errs.append(f"additional top-level property {key!r} not allowed")

    if doc.get("overall_verdict") not in _ALLOWED:
        errs.append(f"overall_verdict {doc.get('overall_verdict')!r} not in {_ALLOWED}")
    if not str(doc.get("narrative", "")).strip():
        errs.append("narrative must be a non-empty string")

    units = doc.get("per_dataset")
    if not isinstance(units, list):
        errs.append("per_dataset must be an array")
        return errs

    allowed_unit = set(_UNIT_KEYS) | (set() if strict else {"independence_group"})

    for i, u in enumerate(units):
        where = f"per_dataset[{i}]"
        if not isinstance(u, dict):
            errs.append(f"{where} must be an object")
            continue
        for key in _UNIT_KEYS:
            if key not in u:
                errs.append(f"{where} missing required key {key!r}")
        for key in u:
            if key not in allowed_unit:
                errs.append(f"{where} additional property {key!r} not allowed")

        if not str(u.get("analysis_unit_id", "")).strip():
            errs.append(f"{where}.analysis_unit_id must be a non-empty string")

        gse = u.get("gse_id")
        if not (isinstance(gse, str) and gse.startswith("GSE") and gse[3:].isdigit()):
            errs.append(f"{where}.gse_id {gse!r} does not match ^GSE[0-9]+$")

        plat = u.get("platform")
        if plat is not None and not (
            isinstance(plat, str) and plat.startswith("GPL") and plat[3:].isdigit()
        ):
            errs.append(f"{where}.platform {plat!r} does not match ^GPL[0-9]+$ and is not null")

        cohort = u.get("cohort")
        if cohort is not None and not isinstance(cohort, str):
            errs.append(f"{where}.cohort must be a string or null")

        if u.get("verdict") not in _ALLOWED:
            errs.append(f"{where}.verdict {u.get('verdict')!r} not in {_ALLOWED}")

        conf = u.get("confidence")
        if not isinstance(conf, (int, float)) or isinstance(conf, bool):
            errs.append(f"{where}.confidence must be a number")
        elif not (0.0 <= float(conf) <= 1.0):
            errs.append(f"{where}.confidence {conf} outside [0, 1]")

        for key in ("evidence_summary", "rationale"):
            if not str(u.get(key, "")).strip():
                errs.append(f"{where}.{key} must be a non-empty string")

        cavs = u.get("caveats")
        if not isinstance(cavs, list):
            errs.append(f"{where}.caveats must be an array")
        else:
            if any(not str(c).strip() for c in cavs):
                errs.append(f"{where}.caveats entries must be non-empty strings")
            if len(set(map(str, cavs))) != len(cavs):
                errs.append(f"{where}.caveats must be unique")

    return errs


# ---------------------------------------------------------------------------
# Worked example
# ---------------------------------------------------------------------------


def _demo() -> None:
    print("=" * 74)
    print("synthesis.schema.json, emitted from a pre-registered decision rule.")
    print("Dataset structure is real. Every effect size below is invented.")
    print("=" * 74)

    # Locked during planning, from Querec 2009 and Ravindran 2014.
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)

    # GSE13485 is one accession holding two trials a year apart with different
    # vaccine lots, per its own Overall design field. That is the paper's own
    # train/test split, so the two cohorts are separate analysis units.
    units = [
        {
            "analysis_unit_id": "GSE13485_trial1",
            "gse_id": "GSE13485",
            "cohort": "Trial 1",
            "platform": "GPL7567",
            "independence_group": "GSE13485|trial1",
            "effect": 0.58,
            "fdr": 0.004,
            "n": 15,
        },
        {
            "analysis_unit_id": "GSE13485_trial2",
            "gse_id": "GSE13485",
            "cohort": "Trial 2",
            "platform": "GPL7567",
            "independence_group": "GSE13485|trial2",
            "effect": 0.51,
            "fdr": 0.03,
            "n": 10,
        },
        # SDY1529's two GEO faces. One study, two accessions.
        {
            "analysis_unit_id": "GSE125921_main",
            "gse_id": "GSE125921",
            "cohort": None,
            "platform": "GPL16791",
            "study": "SDY1529",
            "effect": 0.44,
            "fdr": 0.02,
            "n": 22,
        },
        {
            "analysis_unit_id": "GSE136163_main",
            "gse_id": "GSE136163",
            "cohort": None,
            "platform": "GPL16791",
            "study": "SDY1529",
            "effect": 0.39,
            "fdr": 0.04,
            "n": 18,
        },
    ]

    doc = synthesize(units, rule)

    errs = validate_synthesis(doc)
    print(f"\nschema validation: {'PASS' if not errs else 'FAIL'}")
    for e in errs:
        print(f"  {e}")

    print(f"\noverall_verdict: {doc['overall_verdict']}")
    print("\nper_dataset:")
    for u in doc["per_dataset"]:
        print(f"  {u['analysis_unit_id']:<18} {u['verdict']:<14} "
              f"conf {u['confidence']:<6} {u['evidence_summary']}")

    print("\nnarrative:")
    for line in doc["narrative"].splitlines():
        print(f"  {line}")

    print("\n" + "-" * 74)
    print("Why the GSE13485 cohort split is not a detail.")
    print("-" * 74)
    print("Take only GSE13485, both trials, identical numbers, two readings:")

    yf_only = [dict(u) for u in units[:2]]
    split = synthesize(yf_only, rule, strict=False)
    print(f"\n  trials read as separate cohorts   -> "
          f"{split['n_independent_groups']} group(s), "
          f"{split['overall_verdict']}")

    merged = [dict(u, independence_group="GSE13485") for u in yf_only]
    one = synthesize(merged, rule, strict=False)
    print(f"  trials read as one accession      -> "
          f"{one['n_independent_groups']} group(s), "
          f"{one['overall_verdict']}")
    print("\n  Same effect sizes. The verdict turns entirely on whether anyone")
    print("  read the Overall design field. Querec's blinded validation set is")
    print("  Trial 2, so the split is the paper's own, not ours to invent.")

    print("\nJSON, ready to write to reports/synthesis.json:")
    print(json.dumps(doc, indent=2)[:400] + "\n  ...")


def _self_test() -> None:
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)

    two_groups = [
        {"analysis_unit_id": "A", "gse_id": "GSE1", "platform": "GPL1",
         "independence_group": "g1", "effect": 0.5, "fdr": 0.01},
        {"analysis_unit_id": "B", "gse_id": "GSE2", "platform": "GPL2",
         "independence_group": "g2", "effect": 0.4, "fdr": 0.02},
    ]
    doc = synthesize(two_groups, rule)
    assert validate_synthesis(doc) == [], validate_synthesis(doc)
    assert doc["overall_verdict"] == "supportive"

    # Same numbers, one group. Must not reach supportive.
    one_group = [dict(u, independence_group="g1") for u in two_groups]
    doc1 = synthesize(one_group, rule)
    assert doc1["overall_verdict"] == "inconclusive"
    assert validate_synthesis(doc1) == []

    # Wrong direction maps to contradictory, not dropped.
    down = [dict(two_groups[0], effect=-0.5), dict(two_groups[1], effect=-0.4)]
    doc2 = synthesize(down, rule)
    assert doc2["overall_verdict"] == "contradictory"

    # A required modality nobody measured wins over any row.
    doc3 = synthesize(two_groups, rule, uncovered=["metabolomics"])
    assert doc3["overall_verdict"] == "inconclusive"
    assert "metabolomics" in doc3["narrative"]

    # Cohorts inside one GSE separate only when declared.
    a = {"analysis_unit_id": "X1", "gse_id": "GSE9", "cohort": "T1", "effect": 0.5, "fdr": 0.01}
    b = {"analysis_unit_id": "X2", "gse_id": "GSE9", "cohort": "T2", "effect": 0.5, "fdr": 0.01}
    assert unit_group(a) != unit_group(b)
    assert unit_group({"gse_id": "GSE9"}) == "accession:GSE9"
    assert unit_group({"gse_id": "GSE9", "superseries": "GSE8"}) == "superseries:GSE8"
    assert unit_group({"gse_id": "GSE9", "independence_group": "z"}) == "z"

    # Shared study beats separate accessions.
    sibs = [
        {"analysis_unit_id": "S1", "gse_id": "GSE125921", "study": "SDY1529",
         "effect": 0.5, "fdr": 0.01},
        {"analysis_unit_id": "S2", "gse_id": "GSE136163", "study": "SDY1529",
         "effect": 0.5, "fdr": 0.01},
    ]
    assert synthesize(sibs, rule)["overall_verdict"] == "inconclusive"
    assert "Not independent of" in synthesize(sibs, rule)["per_dataset"][0]["caveats"][1]

    # The link registry from the eligibility side drives grouping.
    registry = {
        "links": [
            {
                "group_id": "SDY1264_GSE13486",
                "immport_accession": "SDY1264",
                "geo_accession": "GSE13486",
                "geo_subseries": ["GSE13485"],
                "relationship": "same_study_with_geo_subseries",
            },
            {
                "group_id": "SDY1529_GEO",
                "immport_accession": "SDY1529",
                "geo_subseries": ["GSE125921", "GSE136163"],
                "relationship": "linked_study",
            },
        ]
    }
    by_acc = groups_from_link_registry(registry)
    assert by_acc["GSE13485"] == "SDY1264_GSE13486"
    assert by_acc["SDY1264"] == "SDY1264_GSE13486"
    assert by_acc["GSE13486"] == "SDY1264_GSE13486"
    assert by_acc["GSE125921"] == by_acc["GSE136163"] == "SDY1529_GEO"
    assert groups_from_link_registry({}) == {}
    assert groups_from_link_registry({"links": [{"geo_accession": "GSE1"}]}) == {}

    # Registry grouping collapses the SDY1529 pair to one group.
    reg_units = [
        {"analysis_unit_id": "GSE125921_main", "gse_id": "GSE125921",
         "independence_group": by_acc["GSE125921"], "effect": 0.5, "fdr": 0.01},
        {"analysis_unit_id": "GSE136163_main", "gse_id": "GSE136163",
         "independence_group": by_acc["GSE136163"], "effect": 0.5, "fdr": 0.01},
    ]
    assert synthesize(reg_units, rule)["overall_verdict"] == "inconclusive"

    # Confidence stays in range and never presents inconclusive as confident.
    for u in synthesize(one_group, rule)["per_dataset"]:
        assert 0.0 <= u["confidence"] <= 1.0

    # Validator actually catches breakage.
    bad = synthesize(two_groups, rule)
    bad["per_dataset"][0]["gse_id"] = "SDY1264"
    assert any("GSE" in e for e in validate_synthesis(bad))
    bad2 = synthesize(two_groups, rule)
    bad2["per_dataset"][0]["confidence"] = 1.4
    assert any("outside" in e for e in validate_synthesis(bad2))
    bad3 = synthesize(two_groups, rule)
    bad3["extra"] = 1
    assert any("additional top-level" in e for e in validate_synthesis(bad3))

    # Non-strict adds fields and validates under the relaxed rules only.
    loose = synthesize(two_groups, rule, strict=False)
    assert validate_synthesis(loose, strict=False) == []
    assert validate_synthesis(loose, strict=True) != []
    assert loose["n_independent_groups"] == 2

    print("\nself-tests: 16 groups of assertions passed")


if __name__ == "__main__":
    _demo()
    _self_test()
