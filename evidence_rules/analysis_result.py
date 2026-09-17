"""One shape for a result, whichever branch produced it.

    python analysis_result.py        # worked example, self-tests

Standard library only.


The problem
-----------

The pipeline forks after eligibility. A dataset with a usable continuous
outcome goes to a regression, one row per subject, and comes back with a
correlation. A dataset whose outcome has been split into groups goes to
Galaxy limma, a probe-by-sample matrix plus a factors file, and comes back
with a logFC per probe.

Both are legitimate tests of the same hypothesis. Neither is interchangeable
with the other. An `effect` of 0.6 means a correlation in one and a log2 fold
change in the other, and nothing in a bare row says which.

`AnalysisResult` is the shape both branches produce. It carries what the
number is, not only what it equals.


The two fields that matter most
-------------------------------

**`effect_type`** says what scale the number is on. Magnitudes are comparable
only within a type. `comparability()` refuses to pool across types rather than
averaging a correlation with a fold change.

**`contrast`** says, for a group comparison, which group the effect is
measured against and whether that group is the high-outcome one. Without it a
positive logFC is ambiguous: it could mean the predictor is higher in strong
responders or in weak ones, depending on how someone wrote the factors file.
`oriented_effect()` uses it to put every result on one convention, positive
meaning more predictor goes with more outcome, so directions can be compared
across branches even when magnitudes cannot.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from typing import Any, Optional, Sequence

__all__ = [
    "ANALYSIS_TYPES",
    "EFFECT_TYPES",
    "AnalysisResult",
    "comparability",
    "from_group_contrast",
    "from_regression",
    "oriented_effect",
    "to_synthesis_unit",
    "validate_result",
]


ANALYSIS_TYPES = ("group_contrast", "regression", "correlation")

# effect_type -> (needs a declared contrast, human label)
EFFECT_TYPES = {
    "logFC":        (True,  "log2 fold change between groups"),
    "cohens_d":     (True,  "standardised mean difference between groups"),
    "pearson_r":    (False, "Pearson correlation"),
    "spearman_rho": (False, "Spearman rank correlation"),
    "beta":         (False, "regression coefficient, per unit predictor"),
    "partial_r":    (False, "partial correlation, covariates held constant"),
}


@dataclass(frozen=True)
class AnalysisResult:
    """One estimate, from one analysis unit, on a declared scale.

    Frozen on purpose. A result is a record of what was run, not a working
    value to be adjusted afterwards.
    """

    analysis_unit_id: str
    gse_id: str
    analysis_type: str
    effect: float
    effect_type: str
    n: int

    predictor: dict[str, Any]
    outcome: dict[str, Any]

    p_value: Optional[float] = None
    fdr: Optional[float] = None

    cohort: Optional[str] = None
    platform: Optional[str] = None
    independence_group: Optional[str] = None

    contrast: Optional[dict[str, Any]] = None
    covariates: tuple[str, ...] = ()
    software: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, **kw: Any) -> str:
        return json.dumps(self.as_dict(), **kw)


# ---------------------------------------------------------------------------
# Orientation
# ---------------------------------------------------------------------------


def oriented_effect(result: AnalysisResult) -> float:
    """The effect signed so that positive means more predictor, more outcome.

    Correlations and regression coefficients already carry that meaning. A
    group contrast does not: its sign depends on which group the factors file
    happened to put first. `contrast["test_is_higher_outcome"]` resolves it,
    and a group contrast without a contrast block has no defined direction,
    so this raises rather than guessing.
    """
    needs_contrast, _ = EFFECT_TYPES[result.effect_type]
    if not needs_contrast:
        return float(result.effect)

    if not result.contrast or "test_is_higher_outcome" not in result.contrast:
        raise ValueError(
            f"{result.analysis_unit_id}: {result.effect_type} needs a contrast "
            "declaring which group is the higher-outcome one. Without it the "
            "sign of the effect is not interpretable."
        )
    higher = result.contrast["test_is_higher_outcome"]
    return float(result.effect) if higher else -float(result.effect)


# ---------------------------------------------------------------------------
# Constructors, one per branch
# ---------------------------------------------------------------------------


def from_group_contrast(
    analysis_unit_id: str,
    gse_id: str,
    effect: float,
    n: int,
    predictor: dict[str, Any],
    outcome: dict[str, Any],
    test_group: str,
    reference_group: str,
    test_is_higher_outcome: bool,
    effect_type: str = "logFC",
    **kw: Any,
) -> AnalysisResult:
    """The Galaxy limma branch: an effect between two groups.

    `test_is_higher_outcome` is not optional and has no default, because a
    wrong value silently flips every conclusion drawn from this result.
    """
    outcome = dict(outcome)
    outcome.setdefault("dichotomised", True)
    return AnalysisResult(
        analysis_unit_id=analysis_unit_id,
        gse_id=gse_id,
        analysis_type="group_contrast",
        effect=float(effect),
        effect_type=effect_type,
        n=int(n),
        predictor=dict(predictor),
        outcome=outcome,
        contrast={
            "test_group": test_group,
            "reference_group": reference_group,
            "test_is_higher_outcome": bool(test_is_higher_outcome),
        },
        **kw,
    )


def from_regression(
    analysis_unit_id: str,
    gse_id: str,
    effect: float,
    n: int,
    predictor: dict[str, Any],
    outcome: dict[str, Any],
    effect_type: str = "pearson_r",
    covariates: Sequence[str] = (),
    **kw: Any,
) -> AnalysisResult:
    """The per-subject branch: an association across subjects.

    The outcome keeps its scale, so `dichotomised` defaults to False.
    """
    outcome = dict(outcome)
    outcome.setdefault("dichotomised", False)
    return AnalysisResult(
        analysis_unit_id=analysis_unit_id,
        gse_id=gse_id,
        analysis_type="correlation" if effect_type in ("pearson_r", "spearman_rho",
                                                       "partial_r") else "regression",
        effect=float(effect),
        effect_type=effect_type,
        n=int(n),
        predictor=dict(predictor),
        outcome=outcome,
        covariates=tuple(covariates),
        **kw,
    )


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_result(r: AnalysisResult) -> list[str]:
    """Return a list of problems. Empty means the record is well formed."""
    errs: list[str] = []

    if r.analysis_type not in ANALYSIS_TYPES:
        errs.append(f"analysis_type {r.analysis_type!r} not in {ANALYSIS_TYPES}")
    if r.effect_type not in EFFECT_TYPES:
        errs.append(f"effect_type {r.effect_type!r} not in {tuple(EFFECT_TYPES)}")
    if not str(r.analysis_unit_id).strip():
        errs.append("analysis_unit_id is empty")
    if not (isinstance(r.n, int) and r.n > 0):
        errs.append(f"n must be a positive integer, got {r.n!r}")
    if r.p_value is None and r.fdr is None:
        errs.append("neither p_value nor fdr is present")
    for name, v in (("p_value", r.p_value), ("fdr", r.fdr)):
        if v is not None and not (0.0 <= float(v) <= 1.0):
            errs.append(f"{name} {v} outside [0, 1]")

    for key in ("name",):
        if key not in r.predictor:
            errs.append(f"predictor is missing {key!r}")
        if key not in r.outcome:
            errs.append(f"outcome is missing {key!r}")

    if r.effect_type in ("pearson_r", "spearman_rho", "partial_r"):
        if not (-1.0 <= float(r.effect) <= 1.0):
            errs.append(f"{r.effect_type} of {r.effect} is outside [-1, 1]")

    needs_contrast, _ = EFFECT_TYPES.get(r.effect_type, (False, ""))
    if needs_contrast:
        if not r.contrast:
            errs.append(f"{r.effect_type} requires a contrast block")
        else:
            for key in ("test_group", "reference_group", "test_is_higher_outcome"):
                if key not in r.contrast:
                    errs.append(f"contrast is missing {key!r}")
    elif r.contrast:
        errs.append(f"{r.effect_type} is not a group comparison but carries a contrast")

    if r.effect_type == "partial_r" and not r.covariates:
        errs.append("partial_r without any covariates listed")

    if r.outcome.get("dichotomised") and r.analysis_type != "group_contrast":
        errs.append("outcome marked dichotomised but analysis_type is not group_contrast")

    return errs


# ---------------------------------------------------------------------------
# Comparability
# ---------------------------------------------------------------------------


POOLING_NOTE = (
    "Effect magnitudes are comparable only within an effect_type. Directions "
    "are comparable across types once oriented, so these results can support "
    "or contradict each other, but their effect sizes must not be averaged."
)


def comparability(results: Sequence[AnalysisResult]) -> dict[str, Any]:
    """Can these results be pooled, and if not, what can still be done.

    Returns what it found plus a note fit to print in a report. It never
    silently harmonises anything.
    """
    types = sorted({r.effect_type for r in results})
    dich = sorted({bool(r.outcome.get("dichotomised")) for r in results})
    predictors = sorted({r.predictor.get("name", "?") for r in results})
    outcomes = sorted({r.outcome.get("name", "?") for r in results})

    can_pool = len(types) == 1
    notes: list[str] = []

    if not can_pool:
        notes.append(
            f"{len(types)} different effect types present ({', '.join(types)}). "
            + POOLING_NOTE
        )
    if len(dich) > 1:
        notes.append(
            "Some results dichotomised the outcome and others kept it continuous. "
            "The dichotomised ones discard within-group ordering and are typically "
            "less powered, so a null there is weaker evidence than a null from the "
            "continuous branch."
        )
    if len(predictors) > 1:
        notes.append(f"Predictors differ across results: {', '.join(predictors)}.")
    if len(outcomes) > 1:
        notes.append(f"Outcomes differ across results: {', '.join(outcomes)}.")

    return {
        "n_results": len(results),
        "effect_types": types,
        "can_pool_magnitudes": can_pool,
        "directions_comparable": True,
        "predictors": predictors,
        "outcomes": outcomes,
        "notes": notes or ["All results share one effect type and one predictor "
                           "and outcome definition."],
    }


# ---------------------------------------------------------------------------
# Handing off to synthesis
# ---------------------------------------------------------------------------


def to_synthesis_unit(result: AnalysisResult) -> dict[str, Any]:
    """Shape one result for `synthesis.synthesize`.

    The effect handed over is the oriented one, so a group contrast and a
    correlation agree on what a positive sign means.
    """
    errs = validate_result(result)
    if errs:
        raise ValueError(f"{result.analysis_unit_id}: " + "; ".join(errs))

    unit: dict[str, Any] = {
        "analysis_unit_id": result.analysis_unit_id,
        "gse_id": result.gse_id,
        "cohort": result.cohort,
        "platform": result.platform,
        "effect": oriented_effect(result),
        "n": result.n,
        "effect_type": result.effect_type,
    }
    if result.independence_group:
        unit["independence_group"] = result.independence_group
    if result.fdr is not None:
        unit["fdr"] = result.fdr
    else:
        unit["p_value"] = result.p_value
    return unit


# ---------------------------------------------------------------------------
# Worked example
# ---------------------------------------------------------------------------


def _demo() -> None:
    print("=" * 74)
    print("One hypothesis, two branches, two shapes of result.")
    print("Effect sizes below are from real runs; see the PR for provenance.")
    print("=" * 74)

    # The regression branch: SDY1264, continuous CD8 outcome, per subject.
    reg = from_regression(
        analysis_unit_id="GSE13485_trial1",
        gse_id="GSE13485",
        effect=0.545,
        p_value=0.044,
        n=15,
        effect_type="partial_r",
        covariates=("EIF2AK4_day0",),
        predictor={"name": "EIF2AK4", "feature_id": "Hs.412102_at",
                   "timepoint": "day 7", "transform": "none"},
        outcome={"name": "Act CD8 T Cell Response", "timepoint": "day 15",
                 "scale": "percentage"},
        cohort="Trial 1",
        platform="GPL7567",
        independence_group="SDY1264|ARM4368",
        software={"tool": "scipy.stats", "version": "1.18.1"},
    )

    # The Galaxy branch: GSE125921, outcome split at the median.
    grp = from_group_contrast(
        analysis_unit_id="GSE125921_all",
        gse_id="GSE125921",
        effect=0.229,
        fdr=4.70e-03,
        n=36,
        effect_type="logFC",
        test_group="high",
        reference_group="low",
        test_is_higher_outcome=True,
        predictor={"name": "UTY", "feature_id": "ILMN_1739587",
                   "timepoint": "baseline", "transform": "log2 quantile"},
        outcome={"name": "day-84 neutralising antibody", "timepoint": "day 84",
                 "scale": "median split at 1280"},
        platform="GPL10558",
        independence_group="SDY1529_GEO",
        software={"tool": "limma", "version": "3.58.1+galaxy0"},
    )

    for r in (reg, grp):
        errs = validate_result(r)
        print(f"\n{r.analysis_unit_id}")
        print(f"  {r.analysis_type:<15} {r.effect_type:<13} effect {r.effect:+.3f}  n={r.n}")
        print(f"  oriented effect: {oriented_effect(r):+.3f}")
        print(f"  dichotomised outcome: {bool(r.outcome.get('dichotomised'))}")
        print(f"  valid: {'yes' if not errs else errs}")

    print("\nComparability across the two:")
    c = comparability([reg, grp])
    print(f"  effect types        {c['effect_types']}")
    print(f"  pool magnitudes     {c['can_pool_magnitudes']}")
    print(f"  compare directions  {c['directions_comparable']}")
    for note in c["notes"]:
        print(f"  - {note}")

    print("\n" + "-" * 74)
    print("The failure this prevents: a contrast written the other way round.")
    print("-" * 74)
    flipped = from_group_contrast(
        analysis_unit_id="GSE125921_flipped",
        gse_id="GSE125921",
        effect=0.229, fdr=4.70e-03, n=36,
        test_group="low", reference_group="high",
        test_is_higher_outcome=False,
        predictor=grp.predictor, outcome=grp.outcome,
    )
    print(f"  same raw effect  {flipped.effect:+.3f}")
    print(f"  oriented effect  {oriented_effect(flipped):+.3f}   <- sign flips, correctly")
    print("  Identical number, opposite meaning. Only the contrast block says which.")


def _self_test() -> None:
    base = dict(
        analysis_unit_id="u1", gse_id="GSE1", n=10,
        predictor={"name": "GENE"}, outcome={"name": "OUT"},
    )

    r = from_regression(effect=0.5, p_value=0.01, **base)
    assert validate_result(r) == []
    assert r.analysis_type == "correlation"
    assert r.outcome["dichotomised"] is False
    assert oriented_effect(r) == 0.5

    b = from_regression(effect=2.4, p_value=0.01, effect_type="beta", **base)
    assert b.analysis_type == "regression"
    assert validate_result(b) == []          # beta is not bounded by [-1, 1]

    g = from_group_contrast(effect=1.2, fdr=0.01, test_group="high",
                            reference_group="low", test_is_higher_outcome=True, **base)
    assert validate_result(g) == []
    assert g.analysis_type == "group_contrast"
    assert g.outcome["dichotomised"] is True
    assert oriented_effect(g) == 1.2

    # The whole point: same number, reversed contrast, opposite meaning.
    gf = from_group_contrast(effect=1.2, fdr=0.01, test_group="low",
                             reference_group="high", test_is_higher_outcome=False, **base)
    assert oriented_effect(gf) == -1.2

    # A group contrast with no contrast block has no defined direction.
    bad = AnalysisResult(analysis_unit_id="u", gse_id="GSE1",
                         analysis_type="group_contrast", effect=1.0,
                         effect_type="logFC", n=5,
                         predictor={"name": "G"}, outcome={"name": "O"}, fdr=0.01)
    assert any("requires a contrast" in e for e in validate_result(bad))
    try:
        oriented_effect(bad)
        raise AssertionError("should refuse to orient without a contrast")
    except ValueError as e:
        assert "not interpretable" in str(e)

    # A correlation must not carry one.
    stray = AnalysisResult(analysis_unit_id="u", gse_id="GSE1",
                           analysis_type="correlation", effect=0.3,
                           effect_type="pearson_r", n=5,
                           predictor={"name": "G"}, outcome={"name": "O"},
                           p_value=0.2, contrast={"test_group": "a"})
    assert any("carries a contrast" in e for e in validate_result(stray))

    # Range and completeness checks.
    assert any("outside [-1, 1]" in e for e in
               validate_result(from_regression(effect=1.4, p_value=0.1, **base)))
    assert any("neither p_value nor fdr" in e for e in
               validate_result(from_regression(effect=0.1, **base)))
    assert any("positive integer" in e for e in validate_result(
        from_regression(effect=0.1, p_value=0.1,
                        **{**base, "n": 0})))
    assert any("partial_r without any covariates" in e for e in validate_result(
        from_regression(effect=0.1, p_value=0.1, effect_type="partial_r", **base)))

    # Comparability refuses to pool across types but keeps directions usable.
    c = comparability([r, g])
    assert c["can_pool_magnitudes"] is False
    assert c["directions_comparable"] is True
    assert any("different effect types" in n for n in c["notes"])
    assert any("dichotomised" in n for n in c["notes"])
    assert comparability([r, from_regression(effect=0.4, p_value=0.02,
                                             **{**base, "analysis_unit_id": "u2"})]
                         )["can_pool_magnitudes"] is True

    # Handoff carries the oriented effect, not the raw one.
    u = to_synthesis_unit(gf)
    assert u["effect"] == -1.2
    assert u["fdr"] == 0.01 and "p_value" not in u
    assert to_synthesis_unit(r)["p_value"] == 0.01
    try:
        to_synthesis_unit(bad)
        raise AssertionError("should refuse an invalid result")
    except ValueError:
        pass

    # Frozen.
    try:
        r.effect = 0.9          # type: ignore[misc]
        raise AssertionError("AnalysisResult should be immutable")
    except Exception:
        pass

    assert json.loads(r.to_json())["effect_type"] == "pearson_r"

    print("\nself-tests: 14 groups of assertions passed")


if __name__ == "__main__":
    _demo()
    _self_test()
