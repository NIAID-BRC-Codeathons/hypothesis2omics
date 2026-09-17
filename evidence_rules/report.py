"""Evidence report: the last stage, and the only one a human reads.

    python report.py                      # worked example to stdout
    python report.py --self-test          # assertions only

`synthesis.py` emits a schema-valid JSON document. Nobody reads JSON aloud in
a review. This turns it into a Markdown report that a domain expert can grade,
which is what the project's own evaluation criteria ask for.

The generator is deliberately dull. It has no opinions and no adjectives. It
prints what the rule promised, what the data was, what the numbers came out
at, and what was not looked at. Every sentence it writes is either a field
from the inputs or a mechanical consequence of the decision rule. There is no
place in this module where prose is generated that asserts more than the
numbers do, and that is the point: the reader should be able to disagree with
the verdict using only the report.

Three things this refuses to do
-------------------------------

**It will not headline a verdict without a declared primary specification.**
The YF-17D benchmark is the reason. Three predictor specifications on the same
25 subjects give three different overall verdicts. A report that shows one of
them is not wrong, it is a choice, and the choice has to have been made
before the numbers. So `render_report` requires exactly one `Specification`
with `role="primary"` and raises otherwise. Naming the primary spec after
seeing the results is still possible, but you have to do it in the open.

**It will not hide the specifications that lost.** Every non-primary spec is
printed in a sensitivity section with its own verdict. If the primary says
inconclusive and a sensitivity spec says supportive, both appear, adjacent.

**It will not report only what was analysed.** A report listing four analysed
datasets and nothing else implies four datasets were all there was.
`ReportContext.not_analysed` and `uncovered` exist so the denominator is
visible: what was found, what was dropped, at which stage, and why.

Controls
--------

A result carrying `role="positive_control"` or `"negative_control"` never
appears in the evidence table. It gets its own section, labelled as machinery
evidence. `analysis_result.to_synthesis_unit` already refuses controls at the
boundary; this module prints them, separately, because a reviewer does want to
see that the sex-marker check came out right.

Stdlib only.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Optional, Sequence

try:  # package import
    from .evidence_rules import DecisionRule, integrate
    from .synthesis import CONFIDENCE_NOTE, synthesize, unit_group, validate_synthesis
except ImportError:  # run as a script
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from evidence_rules import DecisionRule, integrate  # type: ignore
    from synthesis import (  # type: ignore
        CONFIDENCE_NOTE, synthesize, unit_group, validate_synthesis,
    )

__all__ = [
    "ReportContext",
    "Specification",
    "render_report",
    "report_json",
    "write_report",
    "yf17d_context",
]


# Printed once, in the report, so a reader knows what confidence is not.
CONFIDENCE_DISCLAIMER = (
    "Confidence is a deterministic function of the pre-registered alpha and "
    "the adjusted p-value. It is not a calibrated posterior probability and "
    "should not be read as one."
)

REPRODUCE_NOTE = (
    "This report is generated, not written. Re-running the command above on "
    "the same inputs reproduces it byte for byte apart from the date line."
)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReportContext:
    """Everything the report needs that is not a number.

    hypothesis
        The hypothesis as a human stated it, verbatim. Not the parsed form.
        A reviewer needs to see the sentence that was actually asked about.
    literature
        Dicts with `citation`, optional `pmid`, and `claim`. This is where the
        pre-registered direction comes from, so the report can show that the
        direction was taken from published work rather than guessed from the
        hypothesis text.
    not_analysed
        Dicts with `accession`, `stage`, `reason`. Datasets that entered the
        pipeline and left it before producing a result. Without this the
        report has no denominator.
    uncovered
        Measurement types named in the hypothesis that no eligible dataset
        provides. Non-empty means the overall verdict is inconclusive by
        construction, and the report says so in the verdict section.
    """

    hypothesis: str
    hypothesis_id: Optional[str] = None
    prepared_by: Optional[str] = None
    prepared_on: Optional[str] = None

    repositories: tuple[str, ...] = ()
    search_terms: tuple[str, ...] = ()
    n_candidates: Optional[int] = None

    literature: tuple[dict[str, Any], ...] = ()
    not_analysed: tuple[dict[str, Any], ...] = ()
    uncovered: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()

    command: Optional[str] = None
    commit: Optional[str] = None

    def __post_init__(self) -> None:
        if not str(self.hypothesis).strip():
            raise ValueError("a report needs the hypothesis it is about")
        for i, entry in enumerate(self.not_analysed):
            for key in ("accession", "stage", "reason"):
                if not str(entry.get(key, "")).strip():
                    raise ValueError(
                        f"not_analysed[{i}] needs a non-empty {key!r}. A dataset "
                        "dropped without a recorded stage and reason is exactly "
                        "what this section exists to prevent."
                    )


@dataclass(frozen=True)
class Specification:
    """One predictor specification, and the analysis units it produced.

    role
        "primary" for the pre-registered specification, "sensitivity" for the
        rest. Exactly one primary is required across the set. This is not
        bookkeeping: it is the difference between a result and a choice of
        result.
    prereg_basis
        Why this is the primary, in one line, ideally citing the paper whose
        method it matches. Required on the primary, so the justification is
        recorded next to the verdict it produces.
    """

    spec_id: str
    label: str
    units: tuple[dict[str, Any], ...]
    role: str = "sensitivity"
    prereg_basis: Optional[str] = None
    note: Optional[str] = None

    def __post_init__(self) -> None:
        if self.role not in ("primary", "sensitivity"):
            raise ValueError(
                f"role must be 'primary' or 'sensitivity', got {self.role!r}"
            )
        if not self.units:
            raise ValueError(f"{self.spec_id}: a specification with no units")
        if self.role == "primary" and not str(self.prereg_basis or "").strip():
            raise ValueError(
                f"{self.spec_id} is the primary specification and needs a "
                "prereg_basis. If there is no stated reason for preferring it, "
                "it was chosen after the fact and the report should not "
                "present it as pre-registered."
            )


def _as_specs(specs: Sequence[Specification]) -> tuple[Specification, list[Specification]]:
    """Split into the one primary and the rest. Raises if that is not possible."""
    if not specs:
        raise ValueError("a report needs at least one specification")
    primary = [s for s in specs if s.role == "primary"]
    if len(primary) != 1:
        raise ValueError(
            f"exactly one specification must carry role='primary', found "
            f"{len(primary)}. The overall verdict depends on which predictor "
            "specification counts, and that has to be declared rather than "
            "inferred from which one came out best."
        )
    others = [s for s in specs if s.role != "primary"]
    return primary[0], others


# ---------------------------------------------------------------------------
# Small formatting helpers
# ---------------------------------------------------------------------------


def _num(x: Any, fmt: str = "{:+.3f}") -> str:
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return "n/a"
    return fmt.format(x)


def _p(x: Any) -> str:
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return "n/a"
    return "<0.001" if 0 <= x < 0.001 else f"{x:.3g}"


def _cell(x: Any) -> str:
    """A Markdown table cell. Pipes would break the row."""
    if x is None:
        return "-"
    return str(x).replace("|", "\\|").replace("\n", " ")


def _table(header: Sequence[str], rows: Sequence[Sequence[Any]],
           align: Optional[Sequence[str]] = None) -> list[str]:
    if not rows:
        return []
    align = align or ["---"] * len(header)
    out = ["| " + " | ".join(_cell(h) for h in header) + " |",
           "|" + "|".join(align) + "|"]
    for r in rows:
        out.append("| " + " | ".join(_cell(c) for c in r) + " |")
    return out


def _stat_of(unit: dict[str, Any]) -> tuple[str, Any]:
    """Which significance value this unit carries, and its label."""
    if unit.get("fdr") is not None:
        return "FDR", unit["fdr"]
    return "p", unit.get("p_value")


# ---------------------------------------------------------------------------
# Scoring a specification
# ---------------------------------------------------------------------------


def _resolve(spec: Specification, rule: DecisionRule,
             uncovered: Sequence[str]) -> dict[str, Any]:
    """Synthesize one specification and keep the scored rows alongside.

    Returns the non-strict document, so the independence fields the schema has
    no home for are available to print, plus the strict document, which is
    what gets validated.
    """
    units = list(spec.units)
    doc = synthesize(units, rule, uncovered=uncovered, strict=False)
    strict_doc = synthesize(units, rule, uncovered=uncovered, strict=True)
    errs = validate_synthesis(strict_doc, strict=True)
    if errs:
        raise ValueError(f"{spec.spec_id}: synthesis output is not schema-valid: {errs}")

    groups = {
        (u.get("analysis_unit_id") or u.get("gse_id") or "unknown"): unit_group(u)
        for u in units
    }
    scored = []
    for u in units:
        row = dict(u)
        uid = u.get("analysis_unit_id") or u.get("gse_id") or "unknown"
        row["analysis_unit_id"] = uid
        row["dataset"] = uid  # integrate() keys grouping off `dataset`
        scored.append(rule.score(row))

    # integrate()'s own caveats never reach `per_dataset`; they go into the
    # narrative string. The most important one there is that this is a
    # re-analysis, which a reader needs in the caveats section rather than
    # buried in prose, so the rollup is kept and read from directly.
    rollup = integrate(scored, rule, groups, uncovered)

    return {
        "spec": spec,
        "doc": doc,
        "strict_doc": strict_doc,
        "rollup": rollup,
        "units": units,
        "scored": scored,
        "groups": groups,
        "n_groups": len(set(groups.values())),
        "support_groups": sorted(
            {groups[s["analysis_unit_id"]] for s in scored if s["verdict"] == "supports"}
        ),
    }


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


def _sec_header(ctx: ReportContext, primary: Specification) -> list[str]:
    out = ["# Evidence report", ""]
    out.append(f"**Hypothesis.** {ctx.hypothesis.strip()}")
    out.append("")
    meta = []
    if ctx.hypothesis_id:
        meta.append(("Hypothesis ID", ctx.hypothesis_id))
    meta.append(("Primary specification", f"`{primary.spec_id}`: {primary.label}"))
    if ctx.prepared_by:
        meta.append(("Prepared by", ctx.prepared_by))
    meta.append(("Prepared on", ctx.prepared_on or date.today().isoformat()))
    if ctx.commit:
        meta.append(("Pipeline commit", f"`{ctx.commit}`"))
    out += _table(["", ""], [[k, v] for k, v in meta])
    return out


def _sec_verdict(res: dict[str, Any], rule: DecisionRule,
                 ctx: ReportContext) -> list[str]:
    doc = res["doc"]
    verdict = doc["overall_verdict"]
    out = ["", "## Verdict", "", f"**{verdict.upper()}**", ""]

    n_sup = sum(1 for s in res["scored"] if s["verdict"] == "supports")
    n_ref = sum(1 for s in res["scored"] if s["verdict"] == "refutes")
    out.append(
        f"{len(res['units'])} analysis unit(s) in {res['n_groups']} independent "
        f"group(s). {n_sup} support, {n_ref} contradict, "
        f"{len(res['units']) - n_sup - n_ref} inconclusive. "
        f"Support came from {len(res['support_groups'])} independent group(s); "
        f"{rule.min_independent_groups} were required."
    )

    if ctx.uncovered:
        out += [
            "",
            "The verdict is inconclusive by construction, not by result. "
            f"No eligible dataset measures {', '.join(ctx.uncovered)}, which the "
            "hypothesis names. Every number below measures something other than "
            "what was asked.",
        ]
    elif verdict == "inconclusive":
        out += [
            "",
            "Inconclusive is a verdict, not a failure to produce one. It says "
            "the pre-registered criteria were not met, and it is not the same "
            "claim as a refutation.",
        ]
    return out


def _sec_prereg(rule: DecisionRule, primary: Specification,
                ctx: ReportContext) -> list[str]:
    out = ["", "## What was pre-registered", "",
           "Set during planning, before any of the numbers below existed.", ""]
    d = rule.as_dict()
    out += _table(
        ["Criterion", "Value", "What it means"],
        [
            ["Significance", d["alpha"],
             "a unit must reach this to count either way"],
            ["Multiplicity", d["correction"], "correction applied across features"],
            ["Expected direction", d["direction"],
             "a significant effect the other way scores as a refutation, not as nothing"],
            ["Minimum effect", d["min_effect"],
             "magnitude floor; a p-value on a trivial effect is still trivial"],
            ["Independent groups required", d["min_independent_groups"],
             "how many separate dataset groups must support before the overall "
             "verdict is allowed to be supportive"],
        ],
    )
    out += ["", f"**Primary specification.** `{primary.spec_id}`: {primary.label}. "
            f"{primary.prereg_basis}"]

    if ctx.literature:
        out += ["", "**Where the expected direction comes from.**", ""]
        rows = []
        for lit in ctx.literature:
            ref = lit.get("citation", "?")
            if lit.get("pmid"):
                ref += f" (PMID {lit['pmid']})"
            rows.append([ref, lit.get("claim", "")])
        out += _table(["Source", "Claim relied on"], rows)
        out += ["", "The direction is taken from these, not inferred from the "
                "wording of the hypothesis."]
    return out


def _sec_coverage(ctx: ReportContext) -> list[str]:
    if not (ctx.repositories or ctx.not_analysed or ctx.uncovered
            or ctx.n_candidates is not None):
        return []
    out = ["", "## What was searched, and what was not analysed", ""]
    if ctx.repositories:
        out.append(f"Repositories: {', '.join(ctx.repositories)}.")
    if ctx.search_terms:
        out.append(f"Search terms: {', '.join(f'`{t}`' for t in ctx.search_terms)}.")
    if ctx.n_candidates is not None:
        out.append(f"Candidate datasets considered: {ctx.n_candidates}.")

    if ctx.not_analysed:
        out += ["", "Datasets that entered the pipeline and produced no result:", ""]
        out += _table(
            ["Accession", "Left at stage", "Reason"],
            [[e.get("accession"), e.get("stage"), e.get("reason")]
             for e in ctx.not_analysed],
        )
        out += ["", "These are listed so the analysed set is not mistaken for "
                "everything that was available."]

    if ctx.uncovered:
        out += ["", "**Measurement types named in the hypothesis and absent from "
                f"every eligible dataset:** {', '.join(ctx.uncovered)}."]
    return out


def _sec_evidence(res: dict[str, Any]) -> list[str]:
    doc, groups = res["doc"], res["groups"]
    out = ["", "## Evidence", "",
           f"Specification `{res['spec'].spec_id}`: {res['spec'].label}.", ""]

    rows = []
    for entry, unit in zip(doc["per_dataset"], res["units"]):
        label, stat = _stat_of(unit)
        rows.append([
            entry["analysis_unit_id"],
            entry.get("gse_id"),
            entry.get("cohort"),
            entry.get("platform"),
            unit.get("n"),
            _num(unit.get("effect")),
            f"{label} {_p(stat)}",
            entry["verdict"],
            f"{entry['confidence']:.2f}",
            groups.get(entry["analysis_unit_id"], "-"),
        ])
    out += _table(
        ["Analysis unit", "GSE", "Cohort", "Platform", "n", "Effect",
         "Significance", "Verdict", "Conf.", "Independence group"],
        rows,
        align=["---", "---", "---", "---", "--:", "--:", "--:", "---", "--:", "---"],
    )
    out += ["", CONFIDENCE_DISCLAIMER, "", "Why each unit scored as it did:", ""]
    for entry in doc["per_dataset"]:
        out.append(f"- **{entry['analysis_unit_id']}**: {entry['rationale']}.")
    return out


def _sec_independence(res: dict[str, Any], rule: DecisionRule) -> list[str]:
    groups = res["groups"]
    by_group: dict[str, list[str]] = {}
    for uid, g in sorted(groups.items()):
        by_group.setdefault(g, []).append(uid)

    out = ["", "## Independence accounting", "",
           f"{len(groups)} analysis unit(s) resolve to {res['n_groups']} "
           "independent group(s). Replication is counted in groups, not rows.", ""]
    out += _table(
        ["Independence group", "Units", "Counts as"],
        [[g, ", ".join(m), "1 replication"] for g, m in sorted(by_group.items())],
    )
    if len(groups) > res["n_groups"]:
        out += ["", "Units sharing a group share a study, a lab or a platform. "
                "They are correlated by construction, so they contribute once "
                "toward the replication requirement. Counting them separately "
                "is the easiest way to overstate a result."]
    out += ["", f"Support was required from at least "
            f"{rule.min_independent_groups} group(s) and came from "
            f"{len(res['support_groups'])}"
            + (f": {', '.join(res['support_groups'])}." if res["support_groups"] else ".")]
    return out


def _sec_sensitivity(primary_res: dict[str, Any],
                     other_res: Sequence[dict[str, Any]]) -> list[str]:
    if not other_res:
        return []
    out = ["", "## Sensitivity to the predictor specification", "",
           "The same subjects, analysed under other defensible definitions of "
           "the predictor. The primary row is the pre-registered one.", ""]

    rows = []
    for res in [primary_res, *other_res]:
        spec = res["spec"]
        per_unit = "; ".join(
            f"{e.get('cohort') or e['analysis_unit_id']} "
            f"{_num(u.get('effect'))} {_stat_of(u)[0]}={_p(_stat_of(u)[1])}"
            for e, u in zip(res["doc"]["per_dataset"], res["units"])
        )
        rows.append([
            f"`{spec.spec_id}`" + (" **(primary)**" if spec.role == "primary" else ""),
            spec.label,
            per_unit,
            res["doc"]["overall_verdict"],
        ])
    out += _table(["Specification", "Predictor", "Per unit", "Overall"], rows)

    verdicts = {r["doc"]["overall_verdict"] for r in [primary_res, *other_res]}
    if len(verdicts) > 1:
        out += [
            "",
            f"**The verdict is specification-sensitive.** {len(verdicts)} different "
            f"overall verdicts arise on the same data ({', '.join(sorted(verdicts))}). "
            "The primary specification is the one reported above, because it was "
            "declared first. Nothing here licenses reporting whichever came out "
            "best.",
        ]
    else:
        out += ["", "All specifications agree on the overall verdict, which is "
                "the more reassuring outcome of the two."]

    notes = [f"- `{r['spec'].spec_id}`: {r['spec'].note}"
             for r in other_res if r["spec"].note]
    if notes:
        out += ["", *notes]
    return out


def _sec_comparability(comparability_report: Optional[dict[str, Any]]) -> list[str]:
    if not comparability_report:
        return []
    c = comparability_report
    out = ["", "## Comparability of the results pooled above", ""]
    out += _table(
        ["", ""],
        [
            ["Results supplied", c.get("n_results")],
            ["Counted as evidence", c.get("n_scoring")],
            ["Held out as controls", c.get("n_controls")],
            ["Estimands present", ", ".join(c.get("estimands") or []) or "-"],
            ["Effect types present", ", ".join(c.get("effect_types") or []) or "-"],
            ["Magnitudes poolable", "yes" if c.get("can_pool_magnitudes") else "no"],
            ["Directions comparable", "yes" if c.get("directions_comparable") else "no"],
        ],
    )
    notes = c.get("notes") or []
    if notes:
        out += [""]
        out += [f"- {n}" for n in notes]
    if not c.get("can_pool_magnitudes"):
        out += ["", "Because magnitudes are not poolable, no averaged effect size "
                "appears anywhere in this report. The units can still agree or "
                "disagree in direction, and that is what the verdict uses."]
    return out


def _sec_controls(results: Sequence[Any]) -> list[str]:
    """Controls, from AnalysisResult objects. Never part of the evidence table."""
    if not results:
        return []
    try:  # only needed when results are supplied
        try:
            from .analysis_result import EVIDENCE_ROLES, oriented_effect
        except ImportError:
            from analysis_result import EVIDENCE_ROLES, oriented_effect  # type: ignore
    except ImportError:
        return []

    controls = [r for r in results if not EVIDENCE_ROLES.get(r.role, True)]
    if not controls:
        return []

    rows = []
    for r in controls:
        label, stat = ("FDR", r.fdr) if r.fdr is not None else ("p", r.p_value)
        rows.append([
            r.role, r.analysis_unit_id,
            r.predictor.get("name", "?"), r.outcome.get("name", "?"),
            r.n, _num(oriented_effect(r)), f"{label} {_p(stat)}",
        ])
    return ["", "## Controls", "",
            "These test whether the pipeline works, not whether the hypothesis "
            "holds. They are excluded from the evidence table and from the "
            "verdict.", "",
            *_table(["Role", "Analysis unit", "Predictor", "Outcome", "n",
                     "Effect", "Significance"], rows,
                    align=["---", "---", "---", "---", "--:", "--:", "--:"]),
            "",
            "A positive control with a known answer that comes out right is "
            "evidence the machinery is sound. A negative control is a question "
            "the data cannot answer, and it is expected to return inconclusive; "
            "anything else means the pipeline concludes when it should refuse."]


def _sec_what_would_change(res: dict[str, Any], rule: DecisionRule,
                           ctx: ReportContext) -> list[str]:
    """Mechanical, not editorial. Every line is read off the rule."""
    verdict = res["doc"]["overall_verdict"]
    out = ["", "## What would have changed the verdict", ""]

    if ctx.uncovered:
        out.append(
            f"A dataset measuring {', '.join(ctx.uncovered)} at the timepoint the "
            "hypothesis names. Until one exists, no amount of analysis on the "
            "present data addresses the question as stated."
        )
        return out

    lines: list[str] = []
    shortfall = rule.min_independent_groups - len(res["support_groups"])
    if verdict != "supportive" and res["support_groups"] and shortfall > 0:
        lines.append(
            f"Support from {shortfall} more independent dataset group(s). "
            f"{len(res['support_groups'])} supported and "
            f"{rule.min_independent_groups} were required, so the shortfall is "
            "in replication, not in the size of the effect."
        )

    for s, unit in zip(res["scored"], res["units"]):
        if s["verdict"] != "inconclusive":
            continue
        label, stat = _stat_of(unit)
        if isinstance(stat, (int, float)) and not isinstance(stat, bool) and stat > rule.alpha:
            lines.append(
                f"`{s['analysis_unit_id']}` reaching {label} <= {rule.alpha}. It "
                f"came out at {_p(stat)} on n={unit.get('n')}."
            )
        effect = unit.get("effect")
        if (isinstance(effect, (int, float)) and not isinstance(effect, bool)
                and abs(effect) < rule.min_effect):
            lines.append(
                f"`{s['analysis_unit_id']}` reaching the effect floor of "
                f"{rule.min_effect}. Its effect was {_num(effect)}."
            )

    for s, unit in zip(res["scored"], res["units"]):
        if s["verdict"] == "refutes":
            lines.append(
                f"`{s['analysis_unit_id']}` is significant in the opposite "
                f"direction ({_num(unit.get('effect'))} against an expected "
                f"'{rule.direction}'). That is a contradiction of the "
                "hypothesis, and no change to the other units removes it."
            )

    if verdict == "supportive":
        lines.append(
            "Nothing in the present data. A contradicting result from one more "
            "independent group would move this back to inconclusive, which is "
            "the asymmetry the rule was written with."
        )

    out += [f"- {l}" for l in (lines or ["Nothing determinable from the rule alone."])]
    out += ["", "These are read off the pre-registered rule, not chosen. They are "
            "not a plan to reach a particular verdict."]
    return out


def _sec_caveats(res: dict[str, Any], ctx: ReportContext) -> list[str]:
    """Rollup caveats, then per-unit ones, then the caller's own limitations.

    The confidence note is dropped here because the evidence section already
    states it verbatim next to the column it applies to.
    """
    seen: set[str] = {CONFIDENCE_NOTE}
    cavs: list[str] = []
    per_unit = [c for entry in res["doc"]["per_dataset"]
                for c in entry.get("caveats", [])]
    for c in [*res["rollup"].get("caveats", []), *per_unit, *ctx.limitations]:
        if c not in seen:
            seen.add(c)
            cavs.append(c)
    if not cavs:
        return []
    return ["", "## Caveats", "", *[f"- {c}" for c in cavs]]


def _sec_provenance(ctx: ReportContext, res: dict[str, Any]) -> list[str]:
    out = ["", "## Reproducing this", ""]
    if ctx.command:
        out += ["```", ctx.command, "```", ""]
    out.append(REPRODUCE_NOTE)
    out += ["", "The machine-readable twin of this report is the "
            "`synthesis.schema.json` document below, which is what downstream "
            "tooling should consume. The Markdown above is a rendering of it "
            "and adds no numbers of its own.", "",
            "<details><summary>synthesis.schema.json document</summary>", "",
            "```json", json.dumps(res["strict_doc"], indent=2, sort_keys=True), "```",
            "", "</details>"]
    return out


# ---------------------------------------------------------------------------
# The main calls
# ---------------------------------------------------------------------------


def render_report(
    specs: Sequence[Specification],
    rule: DecisionRule,
    context: ReportContext,
    comparability_report: Optional[dict[str, Any]] = None,
    results: Sequence[Any] = (),
) -> str:
    """Render the evidence report as Markdown.

    `specs` must contain exactly one `role="primary"`. `results` is optional
    and only used to print the control section; controls must not appear in
    any specification's units, and `to_synthesis_unit` already enforces that.
    """
    primary, others = _as_specs(specs)
    uncovered = list(context.uncovered)

    primary_res = _resolve(primary, rule, uncovered)
    other_res = [_resolve(s, rule, uncovered) for s in others]

    parts: list[str] = []
    parts += _sec_header(context, primary)
    parts += _sec_verdict(primary_res, rule, context)
    parts += _sec_prereg(rule, primary, context)
    parts += _sec_coverage(context)
    parts += _sec_evidence(primary_res)
    parts += _sec_independence(primary_res, rule)
    parts += _sec_sensitivity(primary_res, other_res)
    parts += _sec_comparability(comparability_report)
    parts += _sec_controls(results)
    parts += _sec_what_would_change(primary_res, rule, context)
    parts += _sec_caveats(primary_res, context)
    parts += _sec_provenance(context, primary_res)
    return "\n".join(parts).rstrip() + "\n"


def report_json(
    specs: Sequence[Specification],
    rule: DecisionRule,
    context: ReportContext,
    comparability_report: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """The same content as a dict, for anything downstream of the report.

    The primary specification's strict document is reproduced verbatim under
    `synthesis`, so a consumer that only understands `synthesis.schema.json`
    can read that key and ignore the rest.
    """
    primary, others = _as_specs(specs)
    uncovered = list(context.uncovered)
    primary_res = _resolve(primary, rule, uncovered)

    return {
        "hypothesis": context.hypothesis,
        "hypothesis_id": context.hypothesis_id,
        "prepared_on": context.prepared_on or date.today().isoformat(),
        "commit": context.commit,
        "decision_rule": rule.as_dict(),
        "primary_specification": {
            "spec_id": primary.spec_id,
            "label": primary.label,
            "prereg_basis": primary.prereg_basis,
        },
        "overall_verdict": primary_res["doc"]["overall_verdict"],
        "n_analysis_units": len(primary_res["units"]),
        "n_independent_groups": primary_res["n_groups"],
        "n_supporting_groups": len(primary_res["support_groups"]),
        "independence_groups": primary_res["groups"],
        "synthesis": primary_res["strict_doc"],
        "sensitivity": [
            {
                "spec_id": s.spec_id,
                "label": s.label,
                "overall_verdict": _resolve(s, rule, uncovered)["doc"]["overall_verdict"],
            }
            for s in others
        ],
        "not_analysed": [dict(e) for e in context.not_analysed],
        "uncovered": list(context.uncovered),
        "comparability": comparability_report,
        "limitations": list(context.limitations),
    }


def write_report(path: Path, *args: Any, **kw: Any) -> Path:
    """Render and write. Also writes `<path>.json` beside it."""
    path = Path(path)
    md = render_report(*args, **kw)
    path.write_text(md, encoding="utf-8")

    kw.pop("results", None)
    doc = report_json(*args, **kw)
    path.with_suffix(path.suffix + ".json").write_text(
        json.dumps(doc, indent=2, sort_keys=True), encoding="utf-8"
    )
    return path


# ---------------------------------------------------------------------------
# Worked example
# ---------------------------------------------------------------------------


# A frozen fixture, not live output. These are the numbers `run_yf17d.py`
# produced from the repo's committed tables on 2026-09-17, copied here so this
# module's example runs without pandas, scipy or the data directory.
#
# It is a fixture and is labelled as one. If the committed tables change, this
# will not notice. The live path is `run_yf17d.py --report`, which computes the
# numbers and renders the report in one pass; that is what belongs in a
# deliverable. Anything produced from the block below is a demonstration of the
# renderer, not a result.
_YF17D = {
    "day7_raw":      {"Trial 1": (+0.438, 0.103), "Trial 2": (+0.032, 0.930)},
    "delta":         {"Trial 1": (+0.645, 0.009), "Trial 2": (+0.681, 0.030)},
    "day7_adjusted": {"Trial 1": (+0.545, 0.044), "Trial 2": (-0.045, 0.907)},
}
_ARMS = {"Trial 1": ("SDY1264|ARM4368", "trial1", 15),
         "Trial 2": ("SDY1264|ARM4369", "trial2", 10)}


def _yf17d_units(spec: str) -> tuple[dict[str, Any], ...]:
    out = []
    for cohort, (group, slug, n) in _ARMS.items():
        effect, p = _YF17D[spec][cohort]
        out.append({
            "analysis_unit_id": f"GSE13485_{slug}",
            "gse_id": "GSE13485",
            "cohort": cohort,
            "platform": "GPL7567",
            "independence_group": group,
            "effect": effect,
            "p_value": p,
            "n": n,
        })
    return tuple(out)


def _yf17d_specs() -> list[Specification]:
    return [
        Specification(
            spec_id="day7_raw",
            label="EIF2AK4 at day 7, as measured",
            units=_yf17d_units("day7_raw"),
            role="primary",
            prereg_basis=(
                "Querec 2009 built its predictive signature from early "
                "expression as measured, and Ravindran 2014 puts the human "
                "signature peak at day 7. This is the specification that "
                "matches the published method, so it is the one that counts."
            ),
        ),
        Specification(
            spec_id="delta",
            label="EIF2AK4 day 7 minus day 0",
            units=_yf17d_units("delta"),
            note=(
                "A day 7 minus day 0 difference score. Where baseline is "
                "itself correlated with the outcome, a difference score "
                "inherits that correlation by construction, so a positive "
                "result here is not independent of the baseline association. "
                "In this fixture the baseline correlation is r = -0.826, "
                "p = 0.003 in Trial 2."
            ),
        ),
        Specification(
            spec_id="day7_adjusted",
            label="EIF2AK4 at day 7, baseline held constant",
            units=_yf17d_units("day7_adjusted"),
            note=(
                "Partial correlation, day 0 partialled out of both sides. This "
                "removes the baseline coupling the difference score inherits, "
                "and costs a degree of freedom doing so."
            ),
        ),
    ]


def yf17d_context() -> ReportContext:
    return ReportContext(
        hypothesis=(
            "Early expression of EIF2AK4 (GCN2) correlates with the magnitude "
            "of the YF-17D specific CD8+ T cell response."
        ),
        hypothesis_id="H1-CD8",
        prepared_by="evidence_rules/report.py",
        repositories=("ImmPort", "GEO"),
        search_terms=("YF-17D", "yellow fever vaccine", "GCN2", "EIF2AK4"),
        n_candidates=6,
        literature=(
            {"citation": "Querec et al. 2009, Nat Immunol 10(1):116-125",
             "pmid": "19029902",
             "claim": "a signature including EIF2AK4 correlated with and predicted "
                      "CD8+ T cell responses in an independent blinded trial"},
            {"citation": "Ravindran et al. 2014, Science 343(6168):313-317",
             "pmid": "24310610",
             "claim": "early GCN2 expression strongly correlates with the magnitude "
                      "of the later CD8+ T cell response; signature peaks at day 7"},
        ),
        not_analysed=(
            {"accession": "GSE13486", "stage": "retrieval",
             "reason": "SuperSeries containing GSE13485; skipped to avoid "
                       "duplicating the same expression data"},
            {"accession": "GSE125921", "stage": "eligibility",
             "reason": "baseline only, all 36 samples carry a _BL suffix, so no "
                       "post-vaccination predictor timepoint exists"},
            {"accession": "GSE13699", "stage": "analysis-ready preparation",
             "reason": "SDY1289 matrix extracted (22,184 x 126) but no design "
                       "file yet, so limma has not been run"},
            {"accession": "SDY1291", "stage": "parsing",
             "reason": "fetched and parsed, produces no linked GEO analysis unit"},
        ),
        limitations=(
            "One accession supplies both analysis units. The two trials were run "
            "a year apart with different vaccine lots, which is why they are "
            "treated as independent, but they share a lab, a platform and a "
            "protocol.",
            "The step between an eligible dataset and an analysis-ready matrix "
            "with a design file is manual. The stages either side are scripted "
            "and provenance-backed; this one is not yet.",
        ),
        command="python evidence_rules/run_yf17d.py --repo . --report evidence_report.md",
    )


def _demo() -> None:
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    comp = {
        "n_results": 3, "n_scoring": 2, "n_controls": 1,
        "estimands": ["EIF2AK4_day7->Act CD8 T Cell Response"],
        "effect_types": ["pearson_r"],
        "can_pool_magnitudes": True, "directions_comparable": True,
        "notes": [
            "1 result(s) carry a control role (positive_control) and are "
            "excluded from the evidence table. They test the pipeline, not the "
            "hypothesis.",
        ],
    }
    print(render_report(_yf17d_specs(), rule, yf17d_context(),
                        comparability_report=comp))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _self_test() -> None:
    n = 0
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    ctx = yf17d_context()
    specs = _yf17d_specs()

    # -- the primary specification must be declared -------------------------
    try:
        render_report([s for s in specs if s.role != "primary"], rule, ctx)
        raise AssertionError("a report with no primary specification was rendered")
    except ValueError as e:
        assert "exactly one" in str(e)
        n += 1

    two_primaries = [
        Specification(spec_id="a", label="A", units=_yf17d_units("delta"),
                      role="primary", prereg_basis="x"),
        Specification(spec_id="b", label="B", units=_yf17d_units("day7_raw"),
                      role="primary", prereg_basis="y"),
    ]
    try:
        render_report(two_primaries, rule, ctx)
        raise AssertionError("two primary specifications were accepted")
    except ValueError:
        n += 1

    # A primary without a stated basis is a post-hoc choice wearing a label.
    try:
        Specification(spec_id="a", label="A", units=_yf17d_units("delta"),
                      role="primary")
        raise AssertionError("primary accepted with no prereg_basis")
    except ValueError as e:
        assert "prereg_basis" in str(e)
        n += 1

    # -- the headline follows the declared primary, not the best result -----
    md = render_report(specs, rule, ctx)
    assert "# Evidence report" in md
    head = md.split("## What was pre-registered")[0]
    assert "INCONCLUSIVE" in head, "the primary spec is inconclusive; headline must say so"
    assert "SUPPORTIVE" not in head, (
        "the delta specification is supportive but is not the primary; it must "
        "not appear in the verdict section"
    )
    n += 1

    # ...and the specification that did come out supportive is still printed.
    assert "Sensitivity to the predictor specification" in md
    assert "specification-sensitive" in md
    assert "`delta`" in md and "supportive" in md
    n += 1

    # -- reordering the specs changes nothing that matters ------------------
    # Sensitivity rows follow the order the caller passed, which is the
    # caller's to choose. The verdict and the evidence must not.
    rev = render_report(list(reversed(specs)), rule, ctx)
    cut = "## Sensitivity to the predictor specification"
    assert rev.split(cut)[0] == md.split(cut)[0], (
        "the verdict or the evidence table changed when the specs were reordered"
    )
    assert report_json(list(reversed(specs)), rule, ctx)["overall_verdict"] == \
        report_json(specs, rule, ctx)["overall_verdict"]
    n += 1

    # -- coverage: what was dropped is in the report ------------------------
    for acc in ("GSE13486", "GSE125921", "GSE13699", "SDY1291"):
        assert acc in md, f"{acc} was dropped from the pipeline and not reported"
    n += 1

    # A dataset cannot be recorded as dropped without a stage and a reason.
    try:
        ReportContext(hypothesis="h",
                      not_analysed=({"accession": "GSE1", "stage": "eligibility"},))
        raise AssertionError("a drop with no reason was accepted")
    except ValueError as e:
        assert "reason" in str(e)
        n += 1

    # -- the rule is printed in full ---------------------------------------
    prereg = md.split("## What was pre-registered")[1].split("##")[0]
    for token in ("0.05", "BH", "up", "Minimum effect",
                  "Independent groups required"):
        assert token in prereg, f"{token} missing from the pre-registration section"
    assert "19029902" in md and "24310610" in md, "literature basis not printed"
    n += 1

    # -- independence ------------------------------------------------------
    assert "Independence accounting" in md
    assert "SDY1264|ARM4368" in md and "SDY1264|ARM4369" in md
    assert "2 independent group(s)" in md
    n += 1

    # Two units in one group must not read as two replications.
    one_group = tuple(
        {**u, "independence_group": "SDY1264|ARM4368"} for u in _yf17d_units("delta")
    )
    collapsed = render_report(
        [Specification(spec_id="delta", label="delta", units=one_group,
                       role="primary", prereg_basis="test")],
        rule, ReportContext(hypothesis="h"),
    )
    assert "1 independent group(s)" in collapsed
    assert "INCONCLUSIVE" in collapsed.split("## What was pre-registered")[0], (
        "two significant units in one independence group produced a supportive "
        "verdict against min_independent_groups=2"
    )
    n += 1

    # -- uncovered measurement short-circuits the verdict -------------------
    unc_ctx = ReportContext(hypothesis="h", uncovered=("protein abundance",))
    unc = render_report(
        [Specification(spec_id="delta", label="delta", units=_yf17d_units("delta"),
                       role="primary", prereg_basis="test")],
        rule, unc_ctx,
    )
    assert "INCONCLUSIVE" in unc.split("## What was pre-registered")[0]
    assert "inconclusive by construction" in unc
    assert "protein abundance" in unc
    # And the remedy section must ask for data, not for a better p-value.
    tail = unc.split("## What would have changed the verdict")[1]
    assert "protein abundance" in tail and "<=" not in tail
    n += 1

    # -- refutation is reported as such ------------------------------------
    flipped = tuple({**u, "effect": -abs(u["effect"])} for u in _yf17d_units("delta"))
    ref = render_report(
        [Specification(spec_id="flip", label="flipped", units=flipped,
                       role="primary", prereg_basis="test")],
        rule, ReportContext(hypothesis="h"),
    )
    assert "CONTRADICTORY" in ref.split("## What was pre-registered")[0]
    assert "opposite direction" in ref
    n += 1

    # -- the report must be printable on a cp1252 console -------------------
    # A Windows console encodes stdout as cp1252 by default, and printing a
    # character it has no mapping for raises UnicodeEncodeError rather than
    # degrading. An em dash or a <= sign in a table is not worth a crash on
    # half the team's machines, so the rendered report stays ASCII.
    # Checking only the rendered fixture is not enough: a non-ASCII character
    # on a branch the fixture never takes would pass and still crash a real
    # run. So the module's own source is checked instead, which covers every
    # literal whether or not this fixture reaches it.
    src = Path(__file__).read_text(encoding="utf-8")
    offenders = sorted({c for c in src if ord(c) > 127})
    assert not offenders, (
        "non-ASCII characters in report.py: "
        + ", ".join(f"U+{ord(c):04X} {c!r}" for c in offenders)
        + ". A Windows console encodes stdout as cp1252 and raises rather than "
        "degrading, so the report stays ASCII. Use <= for the operator and a "
        "plain hyphen for a dash."
    )
    for text in (md, json.dumps(report_json(specs, rule, ctx))):
        text.encode("cp1252")  # would raise UnicodeEncodeError
    n += 1

    # -- caveats: the rollup's own caveats reach the page -------------------
    cav = md.split("## Caveats")[1].split("##")[0]
    assert "re-analysis of public data" in cav, (
        "integrate() puts the re-analysis caveat in the narrative only; the "
        "report must surface it in the caveats section"
    )
    assert CONFIDENCE_NOTE not in cav, (
        "the confidence note is stated in the evidence section and must not be "
        "repeated as a caveat"
    )
    for lim in ctx.limitations:
        assert lim in cav, "a caller-supplied limitation was dropped"
    # Units collapsing to fewer groups must say so in the caveats, not only in
    # the independence table.
    assert "collapse to 1 independent" in collapsed
    n += 1

    # -- comparability notes reach the page --------------------------------
    comp = {"n_results": 2, "n_scoring": 2, "n_controls": 0,
            "estimands": ["a->b", "c->d"], "effect_types": ["logFC", "pearson_r"],
            "can_pool_magnitudes": False, "directions_comparable": False,
            "notes": ["2 different estimands present."]}
    with_comp = render_report(specs, rule, ctx, comparability_report=comp)
    assert "Comparability of the results pooled above" in with_comp
    assert "2 different estimands present." in with_comp
    assert "no averaged effect size" in with_comp
    n += 1

    # -- controls are printed, and outside the evidence table --------------
    try:
        try:
            from .analysis_result import from_group_contrast, from_regression
        except ImportError:
            from analysis_result import from_group_contrast, from_regression  # type: ignore
    except ImportError:
        from_group_contrast = None  # type: ignore

    if from_group_contrast is not None:
        control = from_group_contrast(
            analysis_unit_id="GSE13485_sexcheck", gse_id="GSE13485",
            effect=0.229, n=87,
            predictor={"name": "UTY"}, outcome={"name": "sex"},
            test_group="M", reference_group="F", test_is_higher_outcome=True,
            fdr=0.0047, role="positive_control",
        )
        scoring = from_regression(
            analysis_unit_id="GSE13485_trial1", gse_id="GSE13485",
            effect=0.438, n=15,
            predictor={"name": "EIF2AK4_day7"},
            outcome={"name": "Act CD8 T Cell Response"},
            p_value=0.103,
        )
        with_ctrl = render_report(specs, rule, ctx, results=[scoring, control])
        assert "## Controls" in with_ctrl
        assert "positive_control" in with_ctrl and "UTY" in with_ctrl
        ev = with_ctrl.split("## Evidence")[1].split("## Independence")[0]
        assert "UTY" not in ev, "a control leaked into the evidence table"
        n += 1

        # A scoring-only result set has no control section at all.
        assert "## Controls" not in render_report(specs, rule, ctx, results=[scoring])
        n += 1

    # -- the JSON twin agrees with the Markdown -----------------------------
    doc = report_json(specs, rule, ctx)
    assert doc["overall_verdict"] == "inconclusive"
    assert doc["primary_specification"]["spec_id"] == "day7_raw"
    assert doc["n_independent_groups"] == 2
    assert doc["n_analysis_units"] == 2
    assert len(doc["sensitivity"]) == 2
    assert {s["overall_verdict"] for s in doc["sensitivity"]} == {
        "supportive", "inconclusive"}
    assert validate_synthesis(doc["synthesis"], strict=True) == []
    assert len(doc["not_analysed"]) == 4
    n += 1

    # The embedded document must be the schema-valid strict one.
    assert "synthesis.schema.json document" in md
    embedded = json.loads(md.split("```json")[1].split("```")[0])
    assert validate_synthesis(embedded, strict=True) == [], (
        "the document embedded in the report is not schema-valid"
    )
    assert embedded["overall_verdict"] == "inconclusive"
    n += 1

    # -- a pipe in a field must not break the table ------------------------
    piped = tuple({**u, "cohort": "Trial 1 | lot A"} for u in _yf17d_units("delta"))
    pmd = render_report(
        [Specification(spec_id="p", label="p", units=piped, role="primary",
                       prereg_basis="test")],
        rule, ReportContext(hypothesis="h"),
    )
    assert "Trial 1 \\| lot A" in pmd
    for line in pmd.splitlines():
        if line.startswith("|") and "---" not in line:
            assert line.count("|") - line.count("\\|") >= 2
    n += 1

    # -- a report needs a hypothesis and a specification -------------------
    try:
        ReportContext(hypothesis="   ")
        raise AssertionError("a report context with no hypothesis was accepted")
    except ValueError:
        n += 1
    try:
        Specification(spec_id="x", label="x", units=())
        raise AssertionError("a specification with no units was accepted")
    except ValueError:
        n += 1
    try:
        render_report([], rule, ctx)
        raise AssertionError("a report with no specifications was rendered")
    except ValueError:
        n += 1

    # -- writing to disk -----------------------------------------------------
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = write_report(Path(td) / "evidence_report.md", specs, rule, ctx)
        assert p.exists() and p.read_text(encoding="utf-8") == md
        side = p.with_suffix(p.suffix + ".json")
        assert side.exists()
        assert json.loads(side.read_text(encoding="utf-8"))["overall_verdict"] == \
            "inconclusive"
        n += 1

    print(f"report.py: {n} assertion groups passed")


def main() -> int:
    # The report is ASCII and a self-test enforces that, but a console whose
    # encoding cannot take what we print should not crash the run. Windows
    # defaults to cp1252.
    try:
        sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
    except (AttributeError, OSError, ValueError):
        pass

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--self-test", action="store_true", help="run assertions only")
    ap.add_argument("-o", "--out", type=Path, help="write the report here")
    args = ap.parse_args()

    if args.self_test:
        _self_test()
        return 0
    if args.out:
        rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
        p = write_report(args.out, _yf17d_specs(), rule, yf17d_context())
        print(f"wrote {p} and {p.with_suffix(p.suffix + '.json')}")
        return 0
    _demo()
    _self_test()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
