"""Connect Amar's eligibility engine to the evidence rules.

The two do different jobs and do not overlap:

    Amar's engine   decides WHICH datasets may be used   (per dataset)
    evidence_rules  decides WHAT THE RESULTS MEAN        (across datasets)

The seam between them is one question his engine does not ask, because it
looks at one dataset at a time: of the datasets that came back eligible, how
many genuinely independent studies do they represent? Five eligible datasets
from one research programme is one confirmation, not five.

    python adapters.py      # worked example, 12 self-tests


Usage
-----

    from adapters import eligible_datasets, groups_from_verdicts
    from evidence_rules import DecisionRule, integrate

    verdicts = [assess_dataset(spec, rec) for rec in records]   # his engine

    usable = eligible_datasets(verdicts)                        # ids only
    groups = groups_from_verdicts(verdicts, records)            # id -> group

    rule   = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    result = integrate([rule.score(r) for r in rows], rule, groups)
"""

from __future__ import annotations

from typing import Any, Iterable, Optional, Sequence

from evidence_rules import DecisionRule, assign_groups, coverage_gaps, independence_group

__all__ = [
    "eligible_datasets",
    "groups_from_verdicts",
    "independence_summary",
    "rule_from_spec",
    "unmet_scientific_criteria",
]


# ---------------------------------------------------------------------------
# Reading his engine's output
# ---------------------------------------------------------------------------


def eligible_datasets(
    verdicts: Iterable[dict[str, Any]],
    include_uncertain: bool = False,
) -> list[str]:
    """Dataset ids his engine cleared.

    `scientific_status` is the field that matters. "uncertain" means evidence
    is missing rather than absent, so it is excluded by default but can be
    included while someone chases the missing metadata.
    """
    wanted = {"eligible"} | ({"uncertain"} if include_uncertain else set())
    return [
        v.get("dataset_id", "UNKNOWN")
        for v in verdicts
        if v.get("scientific_status") in wanted
    ]


def unmet_scientific_criteria(verdicts: Iterable[dict[str, Any]]) -> dict[str, list[str]]:
    """Which criteria are blocking, across all datasets.

    Useful before the analysis runs: if `participant_level_linkage` is unknown
    everywhere, that is one problem to solve once, not per dataset.
    """
    counts: dict[str, list[str]] = {}
    for v in verdicts:
        ds = v.get("dataset_id", "UNKNOWN")
        for cid in list(v.get("failed_scientific_criteria", [])) + list(
            v.get("unresolved_criteria", [])
        ):
            counts.setdefault(cid, []).append(ds)
    return counts


# ---------------------------------------------------------------------------
# The bit his engine cannot see
# ---------------------------------------------------------------------------


def groups_from_verdicts(
    verdicts: Iterable[dict[str, Any]],
    records: Iterable[dict[str, Any]],
    include_uncertain: bool = False,
) -> dict[str, str]:
    """Independence group per eligible dataset.

    `records` are the dataset descriptors carrying `superseries`, `bioproject`,
    `lab`, `platform` and so on. Anything his evidence files already hold works;
    missing fields are skipped rather than guessed.
    """
    usable = set(eligible_datasets(verdicts, include_uncertain))
    by_id = {r.get("accession") or r.get("dataset_id"): r for r in records}
    out: dict[str, str] = {}
    for ds in usable:
        rec = by_id.get(ds)
        out[ds] = independence_group(rec) if rec else f"accession:{ds}"
    return out


def independence_summary(
    verdicts: Sequence[dict[str, Any]],
    records: Sequence[dict[str, Any]],
    rule: DecisionRule,
) -> dict[str, Any]:
    """How much independent evidence the eligible set can possibly give.

    Worth running BEFORE any analysis. If the eligible datasets collapse to
    fewer groups than the rule requires, no result can clear the bar and the
    honest move is to widen the search rather than run the numbers and hope.
    """
    groups = groups_from_verdicts(verdicts, records)
    distinct = sorted(set(groups.values()))
    enough = len(distinct) >= rule.min_independent_groups

    by_group: dict[str, list[str]] = {}
    for ds, g in sorted(groups.items()):
        by_group.setdefault(g, []).append(ds)

    return {
        "n_eligible_datasets": len(groups),
        "n_independent_groups": len(distinct),
        "required": rule.min_independent_groups,
        "sufficient": enough,
        "datasets_by_group": by_group,
        "note": (
            "Enough independent evidence is available in principle."
            if enough
            else (
                f"{len(groups)} eligible dataset(s) collapse to {len(distinct)} "
                f"independent group(s), below the required "
                f"{rule.min_independent_groups}. No result from this set can "
                "reach a supporting verdict. Widen the search before analysing."
            )
        ),
    }


# ---------------------------------------------------------------------------
# Reading the direction out of his spec
# ---------------------------------------------------------------------------


def rule_from_spec(
    spec: dict[str, Any],
    alpha: float = 0.05,
    direction: str = "either",
    min_effect: float = 0.0,
    min_independent_groups: int = 2,
) -> DecisionRule:
    """Build a locked DecisionRule alongside a test specification.

    The spec says what will be tested. The rule says what would count as an
    answer. Building them together, before anything runs, is the whole point.

    `direction` is not inferred from the spec text. "associated with" does not
    state a sign, and guessing one is exactly the kind of quiet assumption
    this module exists to prevent. Set it deliberately or leave it "either".
    """
    if direction not in ("up", "down", "either"):
        raise ValueError(f"direction must be up, down or either, got {direction!r}")
    return DecisionRule(
        alpha=alpha,
        direction=direction,
        min_effect=min_effect,
        min_independent_groups=min_independent_groups,
    )


# ---------------------------------------------------------------------------
# Worked example
# ---------------------------------------------------------------------------


def _demo() -> None:
    print("=" * 72)
    print("Connecting the eligibility engine to the evidence rules.")
    print("Dataset fields below are illustrative, effect sizes are invented.")
    print("=" * 72)

    # Shaped like his engine's output.
    verdicts = [
        {"dataset_id": "SDY1264", "scientific_status": "eligible",
         "execution_status": "ready", "failed_scientific_criteria": [],
         "unresolved_criteria": []},
        {"dataset_id": "GSE13485", "scientific_status": "eligible",
         "execution_status": "ready", "failed_scientific_criteria": [],
         "unresolved_criteria": []},
        {"dataset_id": "GSE13699", "scientific_status": "eligible",
         "execution_status": "needs_data_check", "failed_scientific_criteria": [],
         "unresolved_criteria": ["outcome_data_accessible"]},
        {"dataset_id": "SDY1529", "scientific_status": "uncertain",
         "execution_status": "needs_scientific_resolution",
         "failed_scientific_criteria": [],
         "unresolved_criteria": ["quantitative_cd8_response_available"]},
        {"dataset_id": "GSE82152", "scientific_status": "excluded",
         "execution_status": "not_applicable",
         "failed_scientific_criteria": ["participant_level_linkage"],
         "unresolved_criteria": []},
    ]

    # Whatever descriptors we hold for those datasets.
    records = [
        {"accession": "SDY1264",  "program": "HIPC_YF"},
        {"accession": "GSE13485", "program": "HIPC_YF"},
        {"accession": "GSE13699", "bioproject": "PRJNA_LAUSANNE"},
        {"accession": "SDY1529",  "program": "HIPC_YF"},
        {"accession": "GSE82152", "bioproject": "PRJNA_OTHER"},
    ]

    print("\n1. His engine cleared:", eligible_datasets(verdicts))

    print("\n2. What is blocking, across all datasets:")
    for cid, datasets in sorted(unmet_scientific_criteria(verdicts).items()):
        print(f"   {cid}: {', '.join(datasets)}")

    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)

    print("\n3. Independence check, before running anything:")
    summary = independence_summary(verdicts, records, rule)
    for g, members in summary["datasets_by_group"].items():
        print(f"   {g}: {', '.join(members)}")
    print(f"   -> {summary['n_eligible_datasets']} eligible, "
          f"{summary['n_independent_groups']} independent group(s), "
          f"{summary['required']} required")
    print(f"   {summary['note']}")

    print("\n4. Same set, now with results:")
    groups = groups_from_verdicts(verdicts, records)
    rows = [
        {"dataset": "SDY1264",  "effect": 0.58, "fdr": 0.004, "n": 30},
        {"dataset": "GSE13485", "effect": 0.51, "fdr": 0.02,  "n": 25},
        {"dataset": "GSE13699", "effect": 0.44, "fdr": 0.03,  "n": 22},
    ]
    result = integrate_rows(rows, rule, groups)
    print(f"   verdict: {result['verdict']}")
    print(f"   because: {result['rationale']}")

    print("\n5. Drop the one dataset from the second group:")
    result2 = integrate_rows(rows[:2], rule, groups)
    print(f"   verdict: {result2['verdict']}")
    print(f"   because: {result2['rationale']}")


def integrate_rows(rows, rule, groups, uncovered=None):
    """Score and roll up in one call. Convenience only."""
    from evidence_rules import integrate as _integrate

    return _integrate([rule.score(r) for r in rows], rule, groups, uncovered)


def _self_test() -> None:
    verdicts = [
        {"dataset_id": "A", "scientific_status": "eligible", "failed_scientific_criteria": [], "unresolved_criteria": []},
        {"dataset_id": "B", "scientific_status": "eligible", "failed_scientific_criteria": [], "unresolved_criteria": []},
        {"dataset_id": "C", "scientific_status": "uncertain", "failed_scientific_criteria": [], "unresolved_criteria": ["linkage"]},
        {"dataset_id": "D", "scientific_status": "excluded", "failed_scientific_criteria": ["linkage"], "unresolved_criteria": []},
    ]
    records = [
        {"accession": "A", "program": "P1"},
        {"accession": "B", "program": "P1"},
        {"accession": "C", "program": "P2"},
        {"accession": "D", "bioproject": "X"},
    ]

    assert eligible_datasets(verdicts) == ["A", "B"]
    assert sorted(eligible_datasets(verdicts, include_uncertain=True)) == ["A", "B", "C"]

    blocking = unmet_scientific_criteria(verdicts)
    assert sorted(blocking["linkage"]) == ["C", "D"]

    groups = groups_from_verdicts(verdicts, records)
    assert groups == {"A": "program:P1", "B": "program:P1"}

    # A and B look like two datasets but are one study programme.
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    summary = independence_summary(verdicts, records, rule)
    assert summary["n_eligible_datasets"] == 2
    assert summary["n_independent_groups"] == 1
    assert summary["sufficient"] is False
    assert "Widen the search" in summary["note"]

    # Including the uncertain one adds a second group.
    groups_u = groups_from_verdicts(verdicts, records, include_uncertain=True)
    assert len(set(groups_u.values())) == 2

    # A dataset with no descriptor is assumed independent rather than dropped.
    g_missing = groups_from_verdicts(verdicts, [], )
    assert g_missing == {"A": "accession:A", "B": "accession:B"}

    # Results from one group cannot clear a two group bar.
    rows = [{"dataset": "A", "effect": 0.6, "fdr": 0.01},
            {"dataset": "B", "effect": 0.5, "fdr": 0.02}]
    assert integrate_rows(rows, rule, groups)["verdict"] == "inconclusive"

    # Direction is never inferred from spec text.
    assert rule_from_spec({}).direction == "either"
    try:
        rule_from_spec({}, direction="positive")
        raise AssertionError("should reject an unknown direction")
    except ValueError:
        pass

    assert coverage_gaps(["transcriptomics"], ["transcriptomics"]) == []

    print("\nself-tests: 12 assertions passed")


if __name__ == "__main__":
    _demo()
    _self_test()
