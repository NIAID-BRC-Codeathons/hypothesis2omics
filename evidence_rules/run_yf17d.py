"""End to end: the team's own committed data to a schema-valid verdict.

    python run_yf17d.py --repo /path/to/hypothesis2omics

Reads the three normalized tables already in the repo, joins them per subject,
correlates EIF2AK4 against the CD8 response, and hands the result to the
pre-registered decision rule.

Needs pandas and scipy. The evidence_rules modules need neither.


What it joins
-------------

    data/validator_immport/sample_manifest.tsv      GSM -> subject, day, arm
    data/validator_geo/feature_expression.tsv       EIF2AK4 per GSM
    data/validator_immport/quantitative_outcome.tsv CD8 response per subject

`handoff_adapter.py` already proves these link. It reports counts and overlaps
rather than the joined table, so this builds the table and runs the number.


The cohort split
----------------

The manifest's `arm_name` column holds `Trial1` (15 subjects, ARM4368) and
`Trial2` (10 subjects, ARM4369). That is the same split GSE13485 states in its
Overall design field, and Querec's blinded validation set is Trial 2. So the
two arms are separate analysis units and the independence bar is met from one
accession, honestly.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evidence_rules import DecisionRule          # noqa: E402
from synthesis import synthesize, validate_synthesis  # noqa: E402


PREDICTOR_DAY = 7.0     # Ravindran 2014: human signature kinetics peak at day 7
BASELINE_DAY = 0.0


def load(repo: Path, study: str = "SDY1264") -> pd.DataFrame:
    man = pd.read_csv(repo / "data/validator_immport/sample_manifest.tsv", sep="\t")
    expr = pd.read_csv(repo / "data/validator_geo/feature_expression.tsv", sep="\t")
    out = pd.read_csv(repo / "data/validator_immport/quantitative_outcome.tsv", sep="\t")

    man = man[(man.study_accession == study) & (man.repository_name == "GEO")]
    link = man[[
        "repository_accession", "subject_accession", "study_time_collected", "arm_name"
    ]].rename(columns={
        "repository_accession": "sample_accession",
        "study_time_collected": "day",
        "arm_name": "cohort",
    })

    wide = (
        expr.merge(link, on="sample_accession")
            .pivot_table(index=["subject_accession", "cohort"],
                         columns="day", values="expression_value")
            .reset_index()
    )
    df = wide.merge(
        out[["subject_accession", "result_value"]].rename(columns={"result_value": "cd8"}),
        on="subject_accession",
    )
    df["delta"] = df[PREDICTOR_DAY] - df[BASELINE_DAY]
    return df


def partial_corr(x, y, z):
    """corr(x, y) with z partialled out of both, plus its two-sided p."""
    rxy = stats.pearsonr(x, y)[0]
    rxz = stats.pearsonr(x, z)[0]
    ryz = stats.pearsonr(y, z)[0]
    r = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    dof = len(x) - 3
    t = r * np.sqrt(dof / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), dof)


def units_for(df: pd.DataFrame, spec: str) -> list[dict]:
    """One analysis unit per arm, under the named predictor specification."""
    arms = {"Trial1": ("SDY1264|ARM4368", "Trial 1"),
            "Trial2": ("SDY1264|ARM4369", "Trial 2")}
    units = []
    for arm, (group, label) in arms.items():
        s = df[df.cohort == arm]
        if spec == "delta":
            effect, p = stats.pearsonr(s["delta"], s.cd8)
        elif spec == "day7_adjusted":
            effect, p = partial_corr(s[PREDICTOR_DAY].values, s.cd8.values,
                                     s[BASELINE_DAY].values)
        elif spec == "day7_raw":
            effect, p = stats.pearsonr(s[PREDICTOR_DAY], s.cd8)
        else:
            raise ValueError(spec)
        units.append({
            "analysis_unit_id": f"GSE13485_{arm.lower()}",
            "gse_id": "GSE13485",
            "cohort": label,
            "platform": "GPL7567",
            "independence_group": group,
            "effect": float(effect),
            "p_value": float(p),
            "n": int(len(s)),
        })
    return units


SPECS = [
    ("day7_raw",      "EIF2AK4 at day 7, as measured"),
    ("delta",         "EIF2AK4 day 7 minus day 0"),
    ("day7_adjusted", "EIF2AK4 at day 7, baseline held constant"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", required=True, type=Path,
                    help="checkout of NIAID-BRC-Codeathons/hypothesis2omics")
    args = ap.parse_args()

    df = load(args.repo)
    print("=" * 72)
    print("EIF2AK4 vs Act CD8 T Cell Response, day 15, SDY1264 / GSE13485")
    print("=" * 72)
    print(f"\n{len(df)} subjects with both a predictor and an outcome:")
    print(df.groupby("cohort").size().to_string())

    # Locked from the literature, before any of these numbers existed.
    rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
    print(f"\nPre-registered rule: {rule.as_dict()}")

    print("\nBaseline is not neutral here, which is the whole problem:")
    for arm in ("Trial1", "Trial2"):
        s = df[df.cohort == arm]
        r, p = stats.pearsonr(s[BASELINE_DAY], s.cd8)
        print(f"  {arm}: r(day 0, CD8) = {r:+.3f}, p = {p:.3f}")
    print("  A day 7 minus day 0 score inherits that by construction.")

    for spec, label in SPECS:
        units = units_for(df, spec)
        doc = synthesize(units, rule)
        errs = validate_synthesis(doc)
        assert errs == [], errs
        print(f"\n{label}")
        for u, raw in zip(doc["per_dataset"], units):
            print(f"  {u['cohort']:<9} r={raw['effect']:+.3f}  p={raw['p_value']:.3f}  "
                  f"n={raw['n']:<3} -> {u['verdict']}")
        print(f"  overall: {doc['overall_verdict'].upper()}")

    print("\n" + "=" * 72)
    print("Three specifications, three different overall verdicts, one dataset.")
    print("Which one counts has to be chosen before the numbers, not after.")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
