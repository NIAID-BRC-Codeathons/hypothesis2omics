# evidence_rules

The last stage of the pipeline: given results, what are we allowed to claim?

Three files, no dependencies beyond the standard library, plain dicts in and
out.

```bash
python evidence_rules/evidence_rules.py   # four case demo, 30 self-tests
python evidence_rules/synthesis.py        # schema-valid output, 14 test groups
python evidence_rules/adapters.py         # worked example, 12 self-tests
```

All three run in about a second and need nothing installed.

---

## What is here

| File | Does |
|---|---|
| `evidence_rules.py` | The decision logic. Frozen rule, independence counting, coverage gaps. |
| `synthesis.py` | Emits `schemas/synthesis.schema.json`, backed by the frozen rule. |
| `adapters.py` | Reads the eligibility engine's verdicts and answers the set-level question it cannot. |

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

## Using it

```python
from evidence_rules import DecisionRule
from synthesis import synthesize, validate_synthesis

# during planning, before any result exists
rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)

# after analysis
doc = synthesize(analysis_units, rule)
assert validate_synthesis(doc) == []

doc["overall_verdict"]   # supportive | contradictory | inconclusive
doc["per_dataset"]       # one analysis_unit_verdict each
doc["narrative"]         # rule, grouping and caveats in prose
```

Each analysis unit needs `analysis_unit_id`, `gse_id`, and an effect with
`fdr` or `p_value`. `cohort`, `platform`, `superseries`, `bioproject`, `study`
and `independence_group` are all optional and used when present. Missing
fields are skipped, never guessed.

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

**Every effect size in the demos is invented.** They exist to show the logic,
not a result. The output says so on its first line.
