"""Evidence rules: three checks that keep an evidence table honest.

Standalone, zero dependencies, plain dicts in and out. Drop this file into
any pipeline. It does not care how you found your datasets or how you ran
your analysis, only how you decide what the results mean.

Written for Hypothesis2Omics (NIAID-BRC AI Codeathon 2.0, Project 10) to sit
behind whichever synthesis or validation layer the team settles on.

    python evidence_rules.py      # runs the built-in demo and self-tests


The three checks
----------------

1. LOCK THE RULE BEFORE YOU LOOK.
   `DecisionRule` is frozen at construction. Build it during planning, pass
   it to scoring afterwards. Nothing in between can widen alpha, flip the
   expected direction, or drop the effect-size floor once the numbers are
   in. This is the difference between validating a hypothesis and finding
   something and calling it a validation.

2. COUNT INDEPENDENT GROUPS, NOT ROWS.
   Three results from three subseries of one study are one confirmation,
   not three. Accession numbers all look different, so this is easy to get
   wrong and expensive to get wrong in public.

3. SAY INCONCLUSIVE WHEN YOU MEAN IT.
   A significant effect in the wrong direction is a refutation, not a
   result to quietly drop. A required measurement nobody actually made
   means the question was not answered, whatever the rows say.


How to plug it in
-----------------

Anywhere you currently decide "does this support the hypothesis", call
`rule.score(row)`. Anywhere you roll results up, call `integrate(...)`.

    from evidence_rules import DecisionRule, assign_groups, integrate

    rule   = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    groups = assign_groups(dataset_records)
    result = integrate([rule.score(r) for r in rows], rule, groups)

    result["verdict"]    # "supports" | "refutes" | "inconclusive"
    result["rationale"]  # one sentence saying why
    result["caveats"]    # list of things the report must disclose
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Optional, Sequence

__all__ = [
    "DecisionRule",
    "assign_groups",
    "independence_group",
    "integrate",
    "coverage_gaps",
    "benjamini_hochberg",
]

Verdict = str  # "supports" | "refutes" | "inconclusive" | "not_evaluated"


# ---------------------------------------------------------------------------
# 1. The rule, locked at construction
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DecisionRule:
    """Pre-registered criteria for what counts as support.

    Frozen on purpose. Construct it while planning, before any analysis has
    run. If you find yourself wanting to build a second one after seeing
    results, that is the moment this class exists to make visible.

    direction
        "up", "down", or "either". An effect that is significant but points
        the other way scores as a refutation, not as nothing.
    min_effect
        Magnitude floor. A p-value on a trivially small effect is still a
        trivially small effect.
    min_independent_groups
        How many genuinely separate dataset groups must support before the
        overall verdict is allowed to be "supports".
    """

    alpha: float = 0.05
    correction: str = "BH"
    min_effect: float = 0.0
    direction: str = "either"
    min_independent_groups: int = 2

    def __post_init__(self) -> None:
        if not 0 < self.alpha < 1:
            raise ValueError(f"alpha must be between 0 and 1, got {self.alpha}")
        if self.direction not in ("up", "down", "either"):
            raise ValueError(f"direction must be up, down or either, got {self.direction!r}")
        if self.min_independent_groups < 1:
            raise ValueError("min_independent_groups must be at least 1")

    # -- scoring one row ----------------------------------------------------

    def score(self, row: dict[str, Any]) -> dict[str, Any]:
        """Judge a single result against the rule.

        `row` is any dict carrying at least `effect` and one of `fdr` or
        `p_value`. Everything else is passed through untouched, so this
        works with whatever row shape your pipeline already uses.

        Adds two keys: `verdict` and `verdict_rationale`. Never mutates the
        input.
        """
        out = dict(row)
        effect = row.get("effect")
        stat = row.get("fdr") if row.get("fdr") is not None else row.get("p_value")
        stat_name = "FDR" if row.get("fdr") is not None else "p"

        if effect is None or stat is None:
            out["verdict"] = "inconclusive"
            out["verdict_rationale"] = "missing an effect estimate or a significance value"
            return out

        n = row.get("n")
        significant = stat <= self.alpha
        big_enough = abs(effect) >= self.min_effect
        observed = "up" if effect > 0 else "down"
        aligned = self.direction == "either" or observed == self.direction

        if significant and aligned and big_enough:
            out["verdict"] = "supports"
            out["verdict_rationale"] = (
                f"{stat_name}={stat:.3g} at or below {self.alpha}, effect {effect:+.3g} "
                f"in the predicted direction"
            )
        elif significant and not aligned:
            out["verdict"] = "refutes"
            out["verdict_rationale"] = (
                f"{stat_name}={stat:.3g} at or below {self.alpha}, but effect {effect:+.3g} "
                f"runs opposite to the predicted direction ({self.direction})"
            )
        elif significant and not big_enough:
            out["verdict"] = "inconclusive"
            out["verdict_rationale"] = (
                f"{stat_name}={stat:.3g} is significant but the effect "
                f"{abs(effect):.3g} is below the pre-set floor of {self.min_effect}"
            )
        else:
            out["verdict"] = "inconclusive"
            out["verdict_rationale"] = f"{stat_name}={stat:.3g} is above {self.alpha}"

        if isinstance(n, int) and n < 10:
            out["verdict_rationale"] += f"; note n={n} is small for a reliable p-value"
        return out

    def as_dict(self) -> dict[str, Any]:
        """For writing the rule into the report, so readers can see it."""
        return {
            "alpha": self.alpha,
            "correction": self.correction,
            "min_effect": self.min_effect,
            "direction": self.direction,
            "min_independent_groups": self.min_independent_groups,
        }


# ---------------------------------------------------------------------------
# 2. Independence
# ---------------------------------------------------------------------------

# Checked in order. The first field present decides the group.
GROUP_KEYS: Sequence[str] = ("superseries", "bioproject", "study", "program")


def independence_group(record: dict[str, Any]) -> str:
    """Return a key shared by datasets that are NOT independent of each other.

    Two datasets sharing a key are correlated by construction: same
    protocol, same normalization, often overlapping participants. Counting
    them separately overstates replication.

    Resolution order:
        superseries  ->  bioproject  ->  study  ->  program
        then submitting lab plus platform
        then the accession itself, which means "assumed independent"

    Pass whatever fields you have. Missing fields are skipped rather than
    guessed at.
    """
    for key in GROUP_KEYS:
        value = record.get(key)
        if value:
            return f"{key}:{value}"

    lab = record.get("lab") or record.get("contact") or record.get("submitter")
    platform = record.get("platform")
    if isinstance(platform, (list, tuple)):
        platform = ",".join(sorted(str(p) for p in platform))
    if lab or platform:
        return f"lab:{lab or '?'}|platform:{platform or '?'}"

    return f"accession:{record.get('accession', 'unknown')}"


def assign_groups(records: Iterable[dict[str, Any]]) -> dict[str, str]:
    """Map accession -> independence group, for a set of dataset records."""
    return {
        r.get("accession", f"unknown_{i}"): independence_group(r)
        for i, r in enumerate(records)
    }


# ---------------------------------------------------------------------------
# 3. Coverage: was the question actually measured
# ---------------------------------------------------------------------------


def coverage_gaps(
    required: Iterable[str], available: Iterable[str]
) -> list[str]:
    """Required measurement types that nothing eligible actually provides.

    The failure this prevents: a hypothesis about metabolites gets
    "supported" by transcriptomics rows, because the pipeline scored
    whatever it had rather than what was asked for.
    """
    have = {str(a).lower() for a in available}
    return [r for r in required if str(r).lower() not in have]


# ---------------------------------------------------------------------------
# Rollup
# ---------------------------------------------------------------------------

REANALYSIS_CAVEAT = (
    "This is a re-analysis of public data. Agreeing with a published finding "
    "shows the pipeline reproduces it, which is not independent confirmation "
    "of the underlying biology."
)


def integrate(
    scored_rows: Sequence[dict[str, Any]],
    rule: DecisionRule,
    groups: Optional[dict[str, str]] = None,
    uncovered: Optional[Sequence[str]] = None,
) -> dict[str, Any]:
    """Roll scored rows into one verdict, counting groups rather than rows.

    Returns a dict with `verdict`, `rationale`, `caveats`, and the counts
    behind the decision, so a report can show its working.
    """
    groups = groups or {}
    uncovered = list(uncovered or [])

    def group_of(row: dict[str, Any]) -> str:
        acc = row.get("dataset") or row.get("accession") or "unknown"
        return groups.get(acc, f"accession:{acc}")

    supporting = [r for r in scored_rows if r.get("verdict") == "supports"]
    refuting = [r for r in scored_rows if r.get("verdict") == "refutes"]
    support_groups = {group_of(r) for r in supporting}
    all_groups = {group_of(r) for r in scored_rows}

    caveats: list[str] = []

    if uncovered:
        verdict = "inconclusive"
        rationale = (
            f"nothing eligible measures {', '.join(uncovered)}, so the hypothesis "
            "as stated could not be tested"
        )
        caveats.append(
            f"Measurement types named in the hypothesis but absent from every "
            f"eligible dataset: {', '.join(uncovered)}. Any row below measures "
            "something other than what was asked."
        )
    elif refuting and not supporting:
        verdict = "refutes"
        rationale = (
            f"{len(refuting)} result(s) contradict the predicted direction and none support it"
        )
    elif supporting and not refuting and len(support_groups) >= rule.min_independent_groups:
        verdict = "supports"
        rationale = (
            f"{len(supporting)} supporting result(s) across {len(support_groups)} "
            f"independent dataset group(s), with none contradicting"
        )
    elif supporting:
        verdict = "inconclusive"
        bits = [
            f"{len(supporting)} supporting result(s) but only {len(support_groups)} "
            f"independent group(s), and {rule.min_independent_groups} were required"
        ]
        if refuting:
            bits.append(f"{len(refuting)} result(s) contradict")
        verdict, rationale = "inconclusive", "; ".join(bits)
    else:
        verdict = "inconclusive"
        rationale = "no result met the pre-registered criteria"

    if len(scored_rows) > len(all_groups):
        caveats.append(
            f"{len(scored_rows)} rows collapse to {len(all_groups)} independent "
            "dataset group(s). Rows inside a group share a study, lab or platform "
            "and are not independent replication."
        )
    caveats.append(REANALYSIS_CAVEAT)

    return {
        "verdict": verdict,
        "rationale": rationale,
        "caveats": caveats,
        "n_rows": len(scored_rows),
        "n_supporting": len(supporting),
        "n_refuting": len(refuting),
        "n_independent_support_groups": len(support_groups),
        "decision_rule": rule.as_dict(),
    }


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def benjamini_hochberg(pvals: Sequence[float]) -> list[float]:
    """BH adjusted p-values, order preserved. No scipy needed."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    out = [0.0] * n
    prev = 1.0
    for rank, idx in enumerate(reversed(order), start=1):
        i = n - rank + 1
        prev = min(prev, pvals[idx] * n / i)
        out[idx] = prev
    return out


# ---------------------------------------------------------------------------
# Demo and self-test
# ---------------------------------------------------------------------------


def _demo() -> None:
    print("=" * 72)
    print("evidence_rules demo. Effect sizes below are made up.")
    print("=" * 72)

    # Locked during planning, before any analysis ran.
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    print(f"\nrule locked: {rule.as_dict()}")

    # Case A: three results, but all from one study family.
    print("\n--- A. three supporting results, all from one study ---")
    records = [
        {"accession": "GSE13485", "superseries": "STUDY_X"},
        {"accession": "GSE13699", "superseries": "STUDY_X"},
        {"accession": "GSE82152", "superseries": "STUDY_X"},
    ]
    groups = assign_groups(records)
    rows = [
        {"dataset": "GSE13485", "effect": 0.61, "fdr": 0.002, "n": 30},
        {"dataset": "GSE13699", "effect": 0.55, "fdr": 0.006, "n": 28},
        {"dataset": "GSE82152", "effect": 0.49, "fdr": 0.03, "n": 24},
    ]
    res = integrate([rule.score(r) for r in rows], rule, groups)
    print(f"  verdict: {res['verdict']}")
    print(f"  because: {res['rationale']}")

    # Case B: same evidence, genuinely separate studies.
    print("\n--- B. same numbers, two separate studies ---")
    records2 = [
        {"accession": "GSE13485", "superseries": "STUDY_X"},
        {"accession": "GSE125921", "bioproject": "PRJNA_Y"},
    ]
    groups2 = assign_groups(records2)
    rows2 = [
        {"dataset": "GSE13485", "effect": 0.61, "fdr": 0.002, "n": 30},
        {"dataset": "GSE125921", "effect": 0.44, "fdr": 0.01, "n": 26},
    ]
    res2 = integrate([rule.score(r) for r in rows2], rule, groups2)
    print(f"  verdict: {res2['verdict']}")
    print(f"  because: {res2['rationale']}")

    # Case C: strong effect, wrong direction.
    print("\n--- C. strong and significant, wrong direction ---")
    scored = rule.score({"dataset": "GSE13485", "effect": -0.7, "fdr": 0.0001, "n": 30})
    print(f"  row verdict: {scored['verdict']}")
    print(f"  because: {scored['verdict_rationale']}")

    # Case D: the question was never measured.
    print("\n--- D. hypothesis needs metabolomics, nothing measures it ---")
    gaps = coverage_gaps(["transcriptomics", "metabolomics"], ["transcriptomics"])
    res4 = integrate(
        [rule.score({"dataset": "GSE13485", "effect": 0.9, "fdr": 1e-6, "n": 30})],
        rule, {"GSE13485": "superseries:STUDY_X"}, uncovered=gaps,
    )
    print(f"  verdict: {res4['verdict']}")
    print(f"  because: {res4['rationale']}")


def _self_test() -> None:
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)

    # The rule cannot be edited after the fact.
    try:
        rule.alpha = 0.2  # type: ignore[misc]
        raise AssertionError("DecisionRule must be immutable")
    except Exception as exc:
        assert "frozen" in str(type(exc)).lower() or "FrozenInstance" in type(exc).__name__

    assert DecisionRule().min_independent_groups == 2
    for bad in (dict(alpha=0), dict(alpha=1), dict(direction="sideways"), dict(min_independent_groups=0)):
        try:
            DecisionRule(**bad)
            raise AssertionError(f"should have rejected {bad}")
        except ValueError:
            pass

    # Scoring.
    assert rule.score({"effect": 0.6, "fdr": 0.01})["verdict"] == "supports"
    assert rule.score({"effect": -0.6, "fdr": 0.01})["verdict"] == "refutes"
    assert rule.score({"effect": 0.6, "fdr": 0.40})["verdict"] == "inconclusive"
    assert rule.score({"effect": None, "fdr": 0.01})["verdict"] == "inconclusive"
    assert rule.score({"effect": 0.6})["verdict"] == "inconclusive"
    assert all(rule.score(r)["verdict_rationale"] for r in
               [{"effect": 0.6, "fdr": 0.01}, {"effect": -0.6, "fdr": 0.01}, {"effect": 0.6, "fdr": 0.4}])
    assert "small" in rule.score({"effect": 0.6, "fdr": 0.01, "n": 6})["verdict_rationale"]
    assert rule.score({"effect": 0.6, "p_value": 0.01})["verdict"] == "supports"
    # fdr wins over p_value when both present
    assert rule.score({"effect": 0.6, "fdr": 0.4, "p_value": 0.001})["verdict"] == "inconclusive"
    # input untouched
    src = {"effect": 0.6, "fdr": 0.01}
    rule.score(src)
    assert "verdict" not in src

    # Effect floor.
    strict = DecisionRule(alpha=0.05, direction="up", min_effect=0.5)
    assert strict.score({"effect": 0.2, "fdr": 0.001})["verdict"] == "inconclusive"

    # Grouping.
    assert independence_group({"accession": "A", "superseries": "S1"}) == "superseries:S1"
    assert independence_group({"accession": "A", "bioproject": "P1"}) == "bioproject:P1"
    assert independence_group({"accession": "A", "superseries": "S1", "bioproject": "P1"}) == "superseries:S1"
    assert independence_group({"accession": "A", "lab": "L", "platform": ["GPL2", "GPL1"]}) == "lab:L|platform:GPL1,GPL2"
    assert independence_group({"accession": "A"}) == "accession:A"
    assert independence_group({"accession": "A", "superseries": None}) == "accession:A"

    # Same group cannot clear the bar; different groups can.
    g_same = {"A": "superseries:S1", "B": "superseries:S1"}
    g_diff = {"A": "superseries:S1", "B": "bioproject:P2"}
    rows = [{"dataset": "A", "effect": 0.6, "fdr": 0.01},
            {"dataset": "B", "effect": 0.5, "fdr": 0.02}]
    scored = [rule.score(r) for r in rows]
    assert integrate(scored, rule, g_same)["verdict"] == "inconclusive"
    assert integrate(scored, rule, g_diff)["verdict"] == "supports"

    # Refutation dominates when nothing supports.
    ref = [rule.score({"dataset": "A", "effect": -0.6, "fdr": 0.01})]
    assert integrate(ref, rule, g_same)["verdict"] == "refutes"

    # Coverage gap overrides everything.
    strong = [rule.score({"dataset": "A", "effect": 0.9, "fdr": 1e-9})]
    assert integrate(strong, rule, g_diff, uncovered=["metabolomics"])["verdict"] == "inconclusive"
    assert coverage_gaps(["a", "B"], ["A"]) == ["B"]
    assert coverage_gaps(["a"], ["a", "b"]) == []

    # Empty input is inconclusive, not supporting.
    assert integrate([], rule, {})["verdict"] == "inconclusive"

    # The standing caveat is always present.
    assert any("independent confirmation" in c for c in integrate([], rule, {})["caveats"])

    # A single group is allowed if you say so in advance.
    lenient = DecisionRule(alpha=0.05, direction="up", min_independent_groups=1)
    assert integrate([lenient.score(rows[0])], lenient, g_same)["verdict"] == "supports"

    # BH.
    q = benjamini_hochberg([0.001, 0.01, 0.03, 0.2, 0.9])
    assert q == sorted(q) and all(0 <= v <= 1 for v in q)
    assert benjamini_hochberg([]) == []

    print("\nself-tests: 30 assertions passed")


if __name__ == "__main__":
    _demo()
    _self_test()
