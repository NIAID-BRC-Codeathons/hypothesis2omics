"""H2, the antibody arm: TNFRSF17 against the YF-17D neutralizing antibody titer.

    python run_h2.py --repo . --exploratory          # numbers only, no verdict
    python run_h2.py --repo . --report docs/h2.md    # refuses until the probe is confirmed

Companion to `run_yf17d.py`, which does the CD8 arm. Same data, same frozen
rule, different hypothesis:

    H1  early EIF2AK4  ->  magnitude of the CD8+ T cell response   (run_yf17d.py)
    H2  early TNFRSF17 ->  neutralizing antibody titer             (this file)

Both come from Querec et al. 2009, which reported two *distinct* signatures.
EIF2AK4 being flat against antibody titer, and TNFRSF17 being flat against the
CD8 response, are the expected results rather than failures.


Why this refuses to print a verdict
-----------------------------------

The predictor is one probe on a UniGene Custom CDF array, and we do not have a
platform annotation file. `Hs.2556_at` is TNFRSF17 according to a gene-database
lookup, and it behaves the way a plasmablast marker should in this dataset, but
neither of those is an annotation. If the probe is not TNFRSF17 then every
number below is about some other gene, and a verdict would be a confident
claim resting on an unverified identifier.

So `FeatureIdentity.confirmed_by` must name the annotation source before a
report can be rendered. A biology consistency check is recorded as `evidence`
and is deliberately not sufficient. `--exploratory` prints the numbers without
a verdict, which is what they are worth until someone confirms the probe.

This is the same discipline as `report.py` requiring a declared primary
specification: the thing that would embarrass us is structurally blocked rather
than left to whoever is reading.


What the numbers already show
-----------------------------

The difference-score specification looks strong and negative in both arms,
p=0.005 and p=0.040, which a naive read would report as refuting Querec.
It does not survive:

  - baseline TNFRSF17 is POSITIVELY correlated with titer (+0.42, +0.34), so a
    day-7-minus-day-0 score subtracts that off and is dragged negative by
    construction
  - holding baseline constant, neither arm reaches alpha (p=0.054, p=0.097)
  - Trial 2 fails leave-one-out, going to p=0.140, and its outcome has four
    distinct values across ten people
  - the two arms measure titer at DIFFERENT timepoints, day 60 and day 90, so
    they are not estimates of the same quantity and their agreement is not
    replication

That last point is why this file sets `estimand_id` explicitly per arm rather
than letting it be derived from the outcome name. `comparability()` then
reports two estimands and refuses to pool, and the report prints that.

Needs pandas-free stdlib plus scipy and numpy for the correlations.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Optional

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analysis_result import comparability, from_regression       # noqa: E402
from evidence_rules import DecisionRule                          # noqa: E402
from report import ReportContext, Specification, write_report     # noqa: E402


STUDY = "SDY1264"
GSE = "GSE13485"
PLATFORM = "GPL7567"

# Querec 2009: a signature including TNFRSF17 predicted the neutralizing
# antibody response. Higher predictor, higher outcome. Locked here, before any
# number in this file has been computed.
RULE = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)

PRIMARY_SPEC = "d7"


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FeatureIdentity:
    """Which probe is which gene, and how we know.

    confirmed_by
        The annotation source, named. A platform annotation table, a GEO GPL
        record, or a person who looked it up in one. Empty means unconfirmed
        and blocks report generation.
    evidence
        Weaker support: consistency checks, expected kinetics, prior use in a
        paper. Recorded because it is worth recording, and deliberately not
        able to satisfy the gate. A gene that behaves the way you expect is not
        the same as a gene you have identified.
    """

    probe_id: str
    gene: str
    platform: str
    confirmed_by: Optional[str] = None
    evidence: tuple[str, ...] = ()

    @property
    def is_confirmed(self) -> bool:
        return bool((self.confirmed_by or "").strip())

    def caveat(self) -> str:
        if self.is_confirmed:
            return (f"{self.gene} is probe {self.probe_id} on {self.platform}, "
                    f"confirmed against {self.confirmed_by}.")
        return (f"The predictor is probe {self.probe_id} on {self.platform}, "
                f"taken to be {self.gene}. THIS IS NOT CONFIRMED against a "
                f"platform annotation. If the probe is not {self.gene}, every "
                f"number here is about a different gene.")


# Update confirmed_by when someone checks this against the GPL7567 annotation.
# Do not put the consistency check in that field: it is evidence, not an
# annotation, and it belongs below.
TNFRSF17 = FeatureIdentity(
    probe_id="Hs.2556_at",
    gene="TNFRSF17",
    platform=PLATFORM,
    confirmed_by=None,
    evidence=(
        "GPL7567 is a UniGene-based Custom CDF, so its probe names are UniGene "
        "cluster ids with an _at suffix. A gene-database lookup gives TNFRSF17 "
        "the cluster Hs.2556, which makes Hs.2556_at the expected probe. That "
        "lookup could not be corroborated against a second source.",
        "UniProt Q02223 gives TNFRSF17 GeneID 608, and Q9P2K8 gives EIF2AK4 "
        "GeneID 440275. Neither entry carries a UniGene cross-reference, "
        "because UniGene was retired in 2019 and UniProt dropped those links. "
        "So the identifier this platform uses has outlived the database that "
        "issued it, which is the root of this whole problem.",
        "Hs.2556_at is ranked 114 of 20077 probes by mean day-7-minus-day-0 "
        "change (+0.63 log2, 1.55x, top 0.57%), which is how a plasmablast "
        "marker should behave at the day-7 plasmablast peak after YF-17D. "
        "Many genes are induced at day 7, so this is consistent with the "
        "identification rather than proof of it.",
    ),
)


# Routes already tried for the annotation, so nobody spends the hour again:
#
#   GEO acc.cgi?acc=GPL7567          reCAPTCHA
#   GEO FTP over HTTPS               robots disallowed
#   PMC and Europe PMC, Querec 2009  reCAPTCHA / robots disallowed
#   mygene.info REST                 robots disallowed
#   two UniGene mapping tables       robots / connect timeout
#   CRAN + Bioconductor              proxy denies cloud.r-project.org
#
# What closes it, in rough order of effort:
#
#   1. A browser. Open the GPL7567 record on GEO, answer the CAPTCHA, download
#      the SOFT or annotation file, and look up Hs.2556 and Hs.412102. Two
#      minutes for a human, impossible for anything without a browser.
#   2. R with network access to Bioconductor:
#          BiocManager::install("org.Hs.eg.db")
#          library(org.Hs.eg.db)
#          as.list(org.Hs.egUNIGENE[c("608", "440275")])
#      608 is TNFRSF17 and 440275 is EIF2AK4, both from UniProt. If that
#      returns Hs.2556 and Hs.412102 the pairings are confirmed, and both
#      FeatureIdentity records here can name it.
#   3. data/geo_matrix_parse_module.py already fetches from GEO through a path
#      that works. Pulling the platform table alongside the series matrix would
#      put the annotation in the repo permanently, which is the real fix.

# Kept here for contrast, and because it turns out not to be the contrast it
# looks like. `feature_expression.tsv` has a `gene` column reading EIF2AK4
# next to `Hs.412102_at`, which looks like an annotation and is not: both
# values are hand-entered in `data/validator_handoff_config.json` under
# `geo.gene` and `geo.feature_id`, and the parser's own provenance records them
# under `selection`, meaning they were supplied rather than looked up.
#
# So H1's predictor rests on the same unverified identifier as H2's. The gene
# column is the config value copied through the pipeline, which is how a
# hand-typed value acquires the appearance of provenance. Nothing in the
# repository maps a probe to a gene.
EIF2AK4 = FeatureIdentity(
    probe_id="Hs.412102_at",
    gene="EIF2AK4",
    platform=PLATFORM,
    confirmed_by=None,
    evidence=(
        "Querec 2009 identified EIF2AK4 from this exact array, so a UniGene "
        "Custom CDF v9 U133 Plus 2.0 carries it by construction.",
        "data/validator_handoff_config.json pairs gene EIF2AK4 with feature_id "
        "Hs.412102_at. That pairing is a hand-entered configuration value, not "
        "an annotation lookup, and the parser provenance records it under "
        "`selection` as an input.",
    ),
)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def load(repo: Path, probe: str) -> dict[str, Any]:
    """Subject-level expression by day, plus the titer and arm per subject."""
    man: dict[str, dict[str, Any]] = {}
    with open(repo / "data/validator_immport/sample_manifest.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["study_accession"] != STUDY or r["repository_name"] != "GEO":
                continue
            man[r["repository_accession"]] = {
                "subject": r["subject_accession"],
                "day": float(r["study_time_collected"]),
                "arm": r["arm_accession"],
                "arm_name": r["arm_name"],
            }
    if not man:
        raise SystemExit(f"no {STUDY} GEO rows in the sample manifest")

    titer: dict[str, dict[str, Any]] = {}
    tpath = repo / "data/galaxy_file_input/ImmPort_neut_ab_titer_results.tsv"
    if not tpath.exists():
        raise SystemExit(
            f"{tpath} is missing. Regenerate it with "
            "`python data/immport_neut_ab_export.py`."
        )
    with open(tpath, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["study_accession"] != STUDY:
                continue
            v = r["value_preferred"].strip()
            if not v:
                continue
            titer[r["subject_accession"]] = {
                "titer": float(v),
                "arm": r["arm_accession"],
                "day": float(r["study_time_collected"]),
            }

    mpath = repo / f"data/galaxy_file_input/{GSE}_intensities.csv"
    if not mpath.exists():
        raise SystemExit(
            f"{mpath} is missing. Regenerate it with "
            f"`python data/extract_series_matrix.py --input-dir <geo cache>`."
        )
    expr = None
    with open(mpath, newline="") as fh:
        header = next(fh).rstrip("\r\n").split(",")
        for line in fh:
            if line.startswith(probe + ","):
                vals = line.rstrip("\r\n").split(",")
                expr = {header[i]: float(vals[i]) for i in range(1, len(header))}
                break
    if expr is None:
        raise SystemExit(f"probe {probe} is not in {mpath.name}")

    by_subject: dict[str, dict[float, float]] = {}
    for gsm, meta in man.items():
        if gsm in expr:
            by_subject.setdefault(meta["subject"], {})[meta["day"]] = expr[gsm]

    return {"expression": by_subject, "titer": titer, "manifest": man}


ARMS = (("ARM4368", "Trial 1", "trial1"), ("ARM4369", "Trial 2", "trial2"))

SPECS = (
    ("d7",       "TNFRSF17 at day 7, as measured"),
    ("d3",       "TNFRSF17 at day 3, as measured"),
    ("d0",       "TNFRSF17 at day 0, baseline"),
    ("FC_d7-d0", "TNFRSF17 day 7 minus day 0"),
    ("FC_d3-d0", "TNFRSF17 day 3 minus day 0"),
    ("d7_adj",   "TNFRSF17 at day 7, baseline held constant"),
)

SPEC_NOTES = {
    "d3": "The other early timepoint Querec's signature was built from.",
    "d0": "Pre-vaccination. Included because a baseline association is what "
          "makes a difference score untrustworthy, so it has to be visible.",
    "FC_d7-d0": "A difference score. Where baseline is itself correlated with "
                "the outcome, the difference inherits that correlation by "
                "construction, in the direction opposite to the baseline "
                "association.",
    "FC_d3-d0": "The same construction at the earlier timepoint.",
    "d7_adj": "Partial correlation, day 0 partialled out of both sides. This "
              "removes the coupling the difference score inherits, and costs a "
              "degree of freedom.",
}


def partial_corr(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    """corr(x, y) with z partialled out of both, and its two-sided p."""
    rxy = stats.pearsonr(x, y)[0]
    rxz = stats.pearsonr(x, z)[0]
    ryz = stats.pearsonr(y, z)[0]
    r = (rxy - rxz * ryz) / math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    dof = len(x) - 3
    t = r * math.sqrt(dof / (1 - r ** 2))
    return float(r), float(2 * stats.t.sf(abs(t), dof))


def arm_subjects(data: dict[str, Any], arm: str) -> list[str]:
    """Subjects in this arm with a titer and all three expression timepoints."""
    E, T = data["expression"], data["titer"]
    return sorted(
        s for s in E
        if s in T and T[s]["arm"] == arm and {0.0, 3.0, 7.0} <= set(E[s])
    )


def one_cell(data: dict[str, Any], arm: str, spec: str) -> dict[str, Any]:
    """One arm under one specification."""
    E, T = data["expression"], data["titer"]
    subs = arm_subjects(data, arm)
    y = np.array([math.log2(T[s]["titer"]) for s in subs])
    base = np.array([E[s][0.0] for s in subs])

    if spec == "d7_adj":
        r, p = partial_corr(np.array([E[s][7.0] for s in subs]), y, base)
        effect_type = "partial_r"
    else:
        if spec == "d7":
            x = np.array([E[s][7.0] for s in subs])
        elif spec == "d3":
            x = np.array([E[s][3.0] for s in subs])
        elif spec == "d0":
            x = base
        elif spec == "FC_d7-d0":
            x = np.array([E[s][7.0] - E[s][0.0] for s in subs])
        elif spec == "FC_d3-d0":
            x = np.array([E[s][3.0] - E[s][0.0] for s in subs])
        else:
            raise ValueError(spec)
        rs, ps = stats.spearmanr(x, y)
        r, p = float(rs), float(ps)
        effect_type = "spearman_rho"

    outcome_day = sorted({T[s]["day"] for s in subs})
    return {
        "effect": r,
        "p_value": p,
        "n": len(subs),
        "effect_type": effect_type,
        "outcome_day": outcome_day[0] if len(outcome_day) == 1 else None,
        "outcome_days": outcome_day,
        "subjects": subs,
    }


def leave_one_out(data: dict[str, Any], arm: str, spec: str) -> dict[str, Any]:
    """How much of a result rests on one person. Only meaningful for n under 30."""
    E, T = data["expression"], data["titer"]
    subs = arm_subjects(data, arm)
    if spec.startswith("FC_"):
        day = 7.0 if spec == "FC_d7-d0" else 3.0
        xs = [E[s][day] - E[s][0.0] for s in subs]
    elif spec in ("d0", "d3", "d7"):
        day = {"d0": 0.0, "d3": 3.0, "d7": 7.0}[spec]
        xs = [E[s][day] for s in subs]
    else:
        return {}
    ys = [math.log2(T[s]["titer"]) for s in subs]
    worst_p, worst_r, dropped = -1.0, None, None
    for i in range(len(subs)):
        keep = [j for j in range(len(subs)) if j != i]
        r, p = stats.spearmanr([xs[j] for j in keep], [ys[j] for j in keep])
        if p > worst_p:
            worst_p, worst_r, dropped = float(p), float(r), subs[i]
    return {"worst_p": worst_p, "worst_r": worst_r, "dropped": dropped,
            "survives": worst_p <= RULE.alpha}


# ---------------------------------------------------------------------------
# Wiring into the evidence layer
# ---------------------------------------------------------------------------


def units_for(data: dict[str, Any], spec: str) -> tuple[dict[str, Any], ...]:
    out = []
    for arm, label, slug in ARMS:
        cell = one_cell(data, arm, spec)
        out.append({
            "analysis_unit_id": f"{GSE}_{slug}",
            "gse_id": GSE,
            "cohort": label,
            "platform": PLATFORM,
            "independence_group": f"{STUDY}|{arm}",
            "effect": cell["effect"],
            "p_value": cell["p_value"],
            "n": cell["n"],
        })
    return tuple(out)


def results_for(data: dict[str, Any], spec: str,
                identity: FeatureIdentity) -> list[Any]:
    """AnalysisResult objects, with the estimand named explicitly per arm.

    The two arms measure the outcome at different timepoints. Leaving the
    estimand to be derived from the outcome name would make them look like two
    estimates of one quantity, which is exactly the mistake this is here to
    stop. `comparability()` respects an explicit `estimand_id`.
    """
    out = []
    for arm, label, slug in ARMS:
        cell = one_cell(data, arm, spec)
        day = cell["outcome_day"]
        outcome_name = (f"YF-17D neutralizing antibody titer, day "
                        f"{int(day) if day is not None else 'mixed'}")
        out.append(from_regression(
            analysis_unit_id=f"{GSE}_{slug}",
            gse_id=GSE,
            effect=cell["effect"],
            n=cell["n"],
            predictor={"name": f"{identity.gene}_{spec}",
                       "probe_id": identity.probe_id,
                       "identity_confirmed": identity.is_confirmed},
            outcome={"name": outcome_name, "day": day, "scale": "log2"},
            effect_type=cell["effect_type"],
            covariates=("day 0 expression",) if spec == "d7_adj" else (),
            p_value=cell["p_value"],
            cohort=label,
            platform=PLATFORM,
            independence_group=f"{STUDY}|{arm}",
            hypothesis_id="H2-antibody",
            estimand_id=(f"{identity.gene}_{spec}->nab_titer_day"
                         f"{int(day) if day is not None else 'mixed'}"),
            software={"scipy": stats.__name__},
            provenance={"probe_id": identity.probe_id,
                        "probe_identity_confirmed": identity.is_confirmed},
        ))
    return out


def specifications(data: dict[str, Any]) -> list[Specification]:
    specs = []
    for spec, label in SPECS:
        primary = spec == PRIMARY_SPEC
        specs.append(Specification(
            spec_id=spec,
            label=label,
            units=units_for(data, spec),
            role="primary" if primary else "sensitivity",
            prereg_basis=(
                "Querec 2009 built its antibody signature from early expression "
                "as measured, not from a difference score, and TNFRSF17 "
                "induction in this dataset peaks at day 7. This is the "
                "specification that matches the published method."
            ) if primary else None,
            note=None if primary else SPEC_NOTES.get(spec),
        ))
    return specs


def context(data: dict[str, Any], identity: FeatureIdentity,
            repo: Path, command: str) -> ReportContext:
    coupling = "; ".join(
        f"{label} r={one_cell(data, arm, 'd0')['effect']:+.3f} "
        f"p={one_cell(data, arm, 'd0')['p_value']:.3f}"
        for arm, label, _ in ARMS
    )
    days = {label: one_cell(data, arm, "d7")["outcome_days"]
            for arm, label, _ in ARMS}

    limitations = [
        identity.caveat(),
        "The two arms measure the outcome at different timepoints "
        + ", ".join(f"{k} at day {int(v[0])}" for k, v in days.items() if v)
        + ". They are not estimates of the same quantity, so their agreement "
          "is not replication. The estimand is recorded per arm and pooling is "
          "refused.",
        "Baseline TNFRSF17 is itself positively associated with the outcome "
        f"({coupling}). Any day-7-minus-day-0 score therefore inherits that "
        "association with the sign reversed, which is why the difference-score "
        "specifications look strong and negative.",
        "One accession supplies both analysis units. The two trials ran a year "
        "apart with different vaccine lots, which is why they are treated as "
        "separate groups, but they share a lab, a platform and a protocol.",
        "Querec's antibody signature was a multi-gene classifier. A single gene "
        "taken from it need not correlate on its own, so a null here is not by "
        "itself a failure to reproduce that paper.",
    ]
    for e in identity.evidence:
        limitations.append("Probe identity evidence, not confirmation: " + e)

    return ReportContext(
        hypothesis=(
            "Early expression of TNFRSF17 predicts the magnitude of the "
            "neutralizing antibody response to YF-17D vaccination."
        ),
        hypothesis_id="H2-antibody",
        prepared_by="evidence_rules/run_h2.py",
        repositories=("ImmPort", "GEO"),
        search_terms=("YF-17D", "yellow fever vaccine", "TNFRSF17", "BCMA",
                      "neutralizing antibody"),
        n_candidates=4,
        literature=(
            {"citation": "Querec et al. 2009, Nat Immunol 10(1):116-125",
             "pmid": "19029902",
             "claim": "a signature including TNFRSF17 predicted the neutralizing "
                      "antibody response with up to 100 percent accuracy"},
        ),
        not_analysed=(
            {"accession": "GSE13699", "stage": "feature extraction",
             "reason": "SDY1289 has day-28 titers for 16 subjects and an "
                       "extracted GPL6104 matrix, but the TNFRSF17 probe on "
                       "GPL6104 has not been identified. No platform annotation "
                       "has been supplied, so the second cohort cannot be "
                       "tested and this hypothesis rests on one study."},
            {"accession": "GSE82152", "stage": "feature extraction",
             "reason": "SDY1294 has titers at days 28 and 84 for 21 subjects, "
                       "but the GPL21975 expression matrix has not been "
                       "extracted into the validator input"},
            {"accession": "GSE125921", "stage": "eligibility",
             "reason": "baseline only, all 36 samples carry a _BL suffix, so no "
                       "post-vaccination predictor timepoint exists"},
            {"accession": "GSE136163", "stage": "feature extraction",
             "reason": "the post-vaccination half of SDY1529, days 3/7/14/84, "
                       "with day-0 and day-84 titers for 36 subjects, but not "
                       "yet extracted into the validator input"},
        ),
        limitations=tuple(limitations),
        command=command,
        commit=_commit(repo),
    )


def _commit(repo: Path) -> Optional[str]:
    import subprocess
    try:
        out = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def print_table(data: dict[str, Any], identity: FeatureIdentity) -> None:
    print("=" * 78)
    print(f"H2: {identity.gene} ({identity.probe_id}) vs YF-17D neutralizing "
          f"antibody titer")
    print(f"    {STUDY} / {GSE} / {PLATFORM}")
    print("=" * 78)
    if not identity.is_confirmed:
        print()
        print("  PROBE IDENTITY NOT CONFIRMED")
        print("  " + identity.caveat())
    print()
    for arm, label, _ in ARMS:
        cell = one_cell(data, arm, "d7")
        T = data["titer"]
        titers = [T[s]["titer"] for s in cell["subjects"]]
        print(f"  {label} ({arm}): n={cell['n']}, titer at day "
              f"{', '.join(str(int(d)) for d in cell['outcome_days'])}, "
              f"range {min(titers):.0f}-{max(titers):.0f}, "
              f"{len(set(titers))} distinct values")
    print()
    print(f"  Pre-registered rule: {RULE.as_dict()}")
    print(f"  Pre-registered primary specification: {PRIMARY_SPEC}")
    print()

    head = "  " + "specification".ljust(42)
    for _, label, _ in ARMS:
        head += label.rjust(24)
    print(head)
    for spec, label in SPECS:
        row = "  " + label.ljust(42)
        for arm, _, _ in ARMS:
            c = one_cell(data, arm, spec)
            cell = "r={:+.3f} p={:.3f} n={}".format(
                c["effect"], c["p_value"], c["n"])
            marker = " *" if c["p_value"] <= RULE.alpha else "  "
            row += (cell + marker).rjust(24)
        print(row)
    print()
    print("  * reaches the pre-registered alpha. Sign matters: the rule expects "
          f"'{RULE.direction}'.")

    print()
    print("  Leave-one-out on the specifications that reach alpha:")
    any_shown = False
    for spec, label in SPECS:
        for arm, arm_label, _ in ARMS:
            c = one_cell(data, arm, spec)
            if c["p_value"] > RULE.alpha:
                continue
            loo = leave_one_out(data, arm, spec)
            if not loo:
                continue
            any_shown = True
            verdict = "survives" if loo["survives"] else "DOES NOT SURVIVE"
            print(f"    {spec:<10} {arm_label:<9} worst p={loo['worst_p']:.3f} "
                  f"when dropping {loo['dropped']} -> {verdict}")
    if not any_shown:
        print("    nothing reached alpha")


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
    except (AttributeError, OSError, ValueError):
        pass

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", required=True, type=Path,
                    help="checkout of NIAID-BRC-Codeathons/hypothesis2omics")
    ap.add_argument("--report", type=Path, metavar="PATH",
                    help="write the evidence report here, plus PATH.json")
    ap.add_argument("--exploratory", action="store_true",
                    help="print the numbers without a verdict, which is what "
                         "they are worth while the probe is unconfirmed")
    args = ap.parse_args()

    identity = TNFRSF17
    data = load(args.repo, identity.probe_id)
    print_table(data, identity)

    if not args.report:
        if not args.exploratory:
            print()
            print("  No report requested. Use --report PATH to write one, or")
            print("  --exploratory to acknowledge these are numbers, not a verdict.")
        return 0

    if not identity.is_confirmed:
        print()
        print("=" * 78)
        print("REFUSING TO WRITE A REPORT")
        print("=" * 78)
        print()
        print(f"  {identity.caveat()}")
        print()
        print("  What is recorded as evidence for the identification:")
        for e in identity.evidence:
            print(f"    - {e}")
        print()
        print("  A verdict on an unconfirmed probe is a confident claim about a")
        print("  gene we have not identified. Set FeatureIdentity.confirmed_by")
        print("  in this file once the GPL7567 annotation has been checked, and")
        print("  name the annotation source there. Do not put the consistency")
        print("  check in that field.")
        print()
        print("  Meanwhile: python run_h2.py --repo . --exploratory")
        return 2

    specs = specifications(data)
    ctx = context(data, identity, args.repo,
                  f"python evidence_rules/run_h2.py --repo . --report {args.report}")
    comp = comparability(results_for(data, PRIMARY_SPEC, identity))
    p = write_report(args.report, specs, RULE, ctx, comparability_report=comp)
    print()
    print(f"wrote {p}")
    print(f"wrote {p.with_suffix(p.suffix + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
