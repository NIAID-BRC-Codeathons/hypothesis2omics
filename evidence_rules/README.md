# evidence_rules

The last stage of the pipeline: given results, what are we allowed to claim,
and what does the report have to disclose?

No dependencies beyond the standard library, plain dicts in and out, except
`run_yf17d.py` which needs pandas and scipy to read the data tables.

```bash
python evidence_rules/evidence_rules.py    # four case demo, 30 self-tests
python evidence_rules/synthesis.py         # schema-valid output, 16 test groups
python evidence_rules/adapters.py          # worked example, 12 self-tests
python evidence_rules/analysis_result.py   # the result contract, 20 test groups
python evidence_rules/limma_adapter.py     # reads a limma topTable, 17 test groups
python evidence_rules/cohorts.py           # cohort structure from sample titles
python evidence_rules/report.py            # the evidence report, 24 test groups
```

Each runs in about a second and needs nothing installed.

To produce the actual report from the repo's committed tables:

```bash
python evidence_rules/run_yf17d.py --repo . --report evidence_report.md
```

That writes the Markdown report and `evidence_report.md.json` beside it.

---

## What is here

| File | Does |
|---|---|
| `evidence_rules.py` | The decision logic. Frozen rule, independence counting, coverage gaps. |
| `analysis_result.py` | One result shape for both branches of the pipeline. Effect type, estimand, orientation, evidence role. |
| `limma_adapter.py` | Reads a limma `topTable` without letting it pick winners. |
| `cohorts.py` | Recovers cohort structure from GEO sample titles. |
| `synthesis.py` | Emits `schemas/synthesis.schema.json`, backed by the frozen rule. |
| `adapters.py` | Reads the eligibility engine's verdicts and answers the set-level question it cannot. |
| `report.py` | Renders the evidence report a human reads and grades. |
| `run_yf17d.py` | End-to-end run on the tables already in `data/`, and the report. |

`synthesis.py` implements `schemas/synthesis.schema.json`, which was specified
but had no implementation. It includes a validator, so nothing needs
`jsonschema` installed to check its own output.

---

## The three checks

**1. Lock the rule before you look.**

`DecisionRule` is frozen at construction. Build it during planning, pass it to
scoring afterwards. Nothing in between can widen alpha, flip the expected
direction, or drop the effect floor once numbers are in. `rule.as_dict()` gives
the report a copy of the criteria to print, so a reader can check that what was
promised is what was applied.

For the YF-17D benchmark the direction is not a guess. Querec 2009 and
Ravindran 2014 both report the GCN2 correlation as positive, so
`direction="up"` is pre-registerable from published work. See
`docs/ground_truth.md`.

**2. Count independent groups, not rows.**

Three results from three subseries of one study are one confirmation. This is
easy to get wrong because the accessions all look different, and it is the
first thing a reviewer pokes at.

Two live examples in our own candidate set:

- GSE125921 and GSE136163 are both the GEO face of SDY1529. One study, two
  accessions.
- GSE13485 is one accession holding two trials run a year apart with different
  vaccine lots, per its own Overall design field. Querec's blinded validation
  set is the second of those trials, so the split is the paper's own.

Grouping therefore happens at the `analysis_unit_id` grain, not the accession
grain. `cohort` separates units inside one GSE; `superseries`, `bioproject`,
`study` and `program` merge units across GSEs. An explicit
`independence_group` on a unit always wins.

**3. Say inconclusive when you mean it.**

A significant effect in the wrong direction is recorded as a refutation rather
than quietly dropped. A required measurement nobody actually made means the
question went unanswered, whatever the rows say. A hypothesis with no
supporting data must return `inconclusive`, and that is a test, not a
disclaimer.

---

## The report

`report.py` turns the synthesis document into Markdown a domain expert can
grade, which is what the project's evaluation criteria ask for. It refuses
three things, and the refusals are the interesting part.

**It will not headline a verdict without a declared primary specification.**
The benchmark is the reason. Three defensible definitions of "early EIF2AK4
expression" on the same 25 subjects give different overall verdicts. So
`render_report` requires exactly one `Specification` with `role="primary"` and
a `prereg_basis` saying why, and raises otherwise. Choosing the primary after
seeing the numbers is still possible; it just has to happen in the open.

For the benchmark the primary is `day7_raw`, expression as measured, because
that is what Querec's signature was built from. It is also the specification
that comes out inconclusive. The difference score, which reads as a clean
replication, is a sensitivity analysis with a known confound, and it is
reported as one.

**It will not hide the specifications that lost.** Every non-primary
specification appears in a sensitivity table with its own verdict, next to the
primary. When they disagree the report says so in those words.

**It will not report only what was analysed.** `ReportContext.not_analysed`
records every dataset that entered the pipeline and left before producing a
result, with the stage it left at and why. A report listing four analysed
datasets and nothing else implies four datasets were all there was. The
constructor rejects an entry missing a stage or a reason.

Controls get their own section. A result carrying `role="positive_control"` or
`"negative_control"` never appears in the evidence table — `to_synthesis_unit`
refuses it at the boundary — but a reviewer does want to see that the
sex-marker check came out right, so the report prints it, labelled as evidence
about the machinery rather than about the hypothesis.

The "what would have changed the verdict" section is read off the decision
rule, not written. It is not a plan for reaching a particular verdict.

```python
from report import ReportContext, Specification, write_report

write_report(
    "evidence_report.md",
    specs,          # exactly one role="primary", with a prereg_basis
    rule,           # the frozen DecisionRule
    context,        # hypothesis, literature, what was not analysed
)
```

`write_report` also writes `evidence_report.md.json`, whose `synthesis` key is
the schema-valid document verbatim, so anything downstream that only
understands `synthesis.schema.json` can read that key and ignore the rest.

---

## Two open questions for the team

**The synthesis schema has no field for independence.** `per_dataset` is a flat
array with `additionalProperties: false`, so an overall verdict computed from
it cannot express "these two units are one study". Under `strict=True` the
grouping is carried in `narrative` and in each unit's `caveats`, which keeps
output valid and keeps the information visible. `strict=False` adds
`independence_group` per unit plus `n_independent_groups` and `decision_rule`
at the top. That is the shape worth considering for schema v1.1.

**Confidence is required but deterministic rules do not produce one.** The
schema requires a float in [0, 1]. This module computes one from the rule and
the adjusted p-value and attaches a caveat to every unit saying it is not a
calibrated posterior. If we want a real confidence number that is a modelling
decision, not a formatting one.

---

## Deliberate limits

**Direction is never inferred from hypothesis text.** "Associated with" does
not state a sign, and guessing one is exactly the assumption this code exists
to prevent. Set it from the literature or leave it `"either"`.

**Scoring is per row and does no multiple testing across rows.** Testing many
genes? Correct first and pass adjusted values in as `fdr`.
`benjamini_hochberg()` is included if you want it without pulling in scipy.

**Magnitudes are never pooled across estimands.** `comparability()` requires
the same estimand and the same effect type before magnitudes may be combined.
Two logFC results, one for EIF2AK4 against CD8 response and one for UTY against
sex, share a unit and do not estimate the same quantity. When pooling is
refused, no averaged effect size appears anywhere in the report.

**The numbers in `run_yf17d.py` are real; the numbers in the other demos are
not.** `run_yf17d.py` reads the repo's committed tables and reports what it
finds. Every other module's demo uses invented or fixture values to show the
logic, and says so on its first line. `report.py` carries a frozen copy of the
`run_yf17d.py` output so its example runs without pandas; that block is
labelled a fixture, and the live path is `run_yf17d.py --report`.
