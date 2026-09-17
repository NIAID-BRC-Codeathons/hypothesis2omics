# Evidence report

**Hypothesis.** Early expression of EIF2AK4 (GCN2) correlates with the magnitude of the YF-17D specific CD8+ T cell response.

|  |  |
|---|---|
| Hypothesis ID | H1-CD8 |
| Primary specification | `day7_raw`: EIF2AK4 at day 7, as measured |
| Prepared by | evidence_rules/report.py |
| Prepared on | 2026-09-17 |
| Pipeline commit | `2ba6aaf` |

## Verdict

**INCONCLUSIVE**

2 analysis unit(s) in 2 independent group(s). 0 support, 0 contradict, 2 inconclusive. Support came from 0 independent group(s); 2 were required.

Inconclusive is a verdict, not a failure to produce one. It says the pre-registered criteria were not met, and it is not the same claim as a refutation.

## What was pre-registered

Set during planning, before any of the numbers below existed.

| Criterion | Value | What it means |
|---|---|---|
| Significance | 0.05 | a unit must reach this to count either way |
| Multiplicity | BH | correction applied across features |
| Expected direction | up | a significant effect the other way scores as a refutation, not as nothing |
| Minimum effect | 0.0 | magnitude floor; a p-value on a trivial effect is still trivial |
| Independent groups required | 2 | how many separate dataset groups must support before the overall verdict is allowed to be supportive |

**Primary specification.** `day7_raw`: EIF2AK4 at day 7, as measured. Querec 2009 built its predictive signature from early expression as measured, and Ravindran 2014 puts the human signature peak at day 7. This is the specification that matches the published method, so it is the one that counts. The two difference-score specifications below are sensitivity analyses.

**Where the expected direction comes from.**

| Source | Claim relied on |
|---|---|
| Querec et al. 2009, Nat Immunol 10(1):116-125 (PMID 19029902) | a signature including EIF2AK4 correlated with and predicted CD8+ T cell responses in an independent blinded trial |
| Ravindran et al. 2014, Science 343(6168):313-317 (PMID 24310610) | early GCN2 expression strongly correlates with the magnitude of the later CD8+ T cell response; signature peaks at day 7 |

The direction is taken from these, not inferred from the wording of the hypothesis.

## What was searched, and what was not analysed

Repositories: ImmPort, GEO.
Search terms: `YF-17D`, `yellow fever vaccine`, `GCN2`, `EIF2AK4`.
Candidate datasets considered: 6.

Datasets that entered the pipeline and produced no result:

| Accession | Left at stage | Reason |
|---|---|---|
| GSE13486 | retrieval | SuperSeries containing GSE13485; skipped to avoid duplicating the same expression data |
| GSE125921 | eligibility | baseline only, all 36 samples carry a _BL suffix, so no post-vaccination predictor timepoint exists |
| GSE13699 | analysis-ready preparation | SDY1289 matrix extracted (22,184 x 126) but no design file yet, so limma has not been run |
| SDY1291 | parsing | fetched and parsed, produces no linked GEO analysis unit |

These are listed so the analysed set is not mistaken for everything that was available.

## Evidence

Specification `day7_raw`: EIF2AK4 at day 7, as measured.

| Analysis unit | GSE | Cohort | Platform | n | Effect | Significance | Verdict | Conf. | Independence group |
|---|---|---|---|--:|--:|--:|---|--:|---|
| GSE13485_trial1 | GSE13485 | Trial 1 | GPL7567 | 15 | +0.438 | p 0.103 | inconclusive | 0.40 | SDY1264\|ARM4368 |
| GSE13485_trial2 | GSE13485 | Trial 2 | GPL7567 | 10 | +0.032 | p 0.93 | inconclusive | 0.07 | SDY1264\|ARM4369 |

Confidence is a deterministic function of the pre-registered alpha and the adjusted p-value. It is not a calibrated posterior probability and should not be read as one.

Why each unit scored as it did:

- **GSE13485_trial1**: p=0.103 is above 0.05.
- **GSE13485_trial2**: p=0.93 is above 0.05.

## Independence accounting

2 analysis unit(s) resolve to 2 independent group(s). Replication is counted in groups, not rows.

| Independence group | Units | Counts as |
|---|---|---|
| SDY1264\|ARM4368 | GSE13485_trial1 | 1 replication |
| SDY1264\|ARM4369 | GSE13485_trial2 | 1 replication |

Support was required from at least 2 group(s) and came from 0.

## Sensitivity to the predictor specification

The same subjects, analysed under other defensible definitions of the predictor. The primary row is the pre-registered one.

| Specification | Predictor | Per unit | Overall |
|---|---|---|---|
| `day7_raw` **(primary)** | EIF2AK4 at day 7, as measured | Trial 1 +0.438 p=0.103; Trial 2 +0.032 p=0.93 | inconclusive |
| `delta` | EIF2AK4 day 7 minus day 0 | Trial 1 +0.645 p=0.00943; Trial 2 +0.681 p=0.03 | supportive |
| `day7_adjusted` | EIF2AK4 at day 7, baseline held constant | Trial 1 +0.545 p=0.0438; Trial 2 -0.045 p=0.907 | inconclusive |

**The verdict is specification-sensitive.** 2 different overall verdicts arise on the same data (inconclusive, supportive). The primary specification is the one reported above, because it was declared first. Nothing here licenses reporting whichever came out best.

- `delta`: A day 7 minus day 0 difference score. Where baseline is itself correlated with the outcome, a difference score inherits that correlation by construction, so a positive result here is not independent of the baseline association. Baseline against outcome: Trial 1 r = -0.439, p = 0.102; Trial 2 r = -0.826, p = 0.003.
- `day7_adjusted`: Partial correlation, day 0 partialled out of both sides. This removes the baseline coupling the difference score inherits, and costs a degree of freedom doing so.

## What would have changed the verdict

- `GSE13485_trial1` reaching p <= 0.05. It came out at 0.103 on n=15.
- `GSE13485_trial2` reaching p <= 0.05. It came out at 0.93 on n=10.

These are read off the pre-registered rule, not chosen. They are not a plan to reach a particular verdict.

## Caveats

- This is a re-analysis of public data. Agreeing with a published finding shows the pipeline reproduces it, which is not independent confirmation of the underlying biology.
- One accession supplies both analysis units. The two trials were run a year apart with different vaccine lots, which is why they are treated as independent, but they share a lab, a platform and a protocol.
- The step between an eligible dataset and an analysis-ready matrix with a design file is manual. The stages either side are scripted and provenance-backed; this one is not yet.

## Reproducing this

```
python evidence_rules/run_yf17d.py --repo . --report docs\evidence_report_yf17d.md
```

This report is generated, not written. Re-running the command above on the same inputs reproduces it byte for byte apart from the date line.

The machine-readable twin of this report is the `synthesis.schema.json` document below, which is what downstream tooling should consume. The Markdown above is a rendering of it and adds no numbers of its own.

<details><summary>synthesis.schema.json document</summary>

```json
{
  "narrative": "Overall verdict: inconclusive, because no result met the pre-registered criteria.\n\nPre-registered before any result existed: alpha 0.05, BH correction, expected direction 'up', minimum effect 0.0, and at least 2 independent dataset group(s) required.\n\n2 analysis unit(s) resolve to 2 independent group(s):\n  SDY1264|ARM4368: GSE13485_trial1\n  SDY1264|ARM4369: GSE13485_trial2\n\nCaveats:\n  - This is a re-analysis of public data. Agreeing with a published finding shows the pipeline reproduces it, which is not independent confirmation of the underlying biology.",
  "overall_verdict": "inconclusive",
  "per_dataset": [
    {
      "analysis_unit_id": "GSE13485_trial1",
      "caveats": [
        "Confidence is derived deterministically from the pre-registered rule and the adjusted p-value. It is not a calibrated posterior probability."
      ],
      "cohort": "Trial 1",
      "confidence": 0.4,
      "evidence_summary": "effect +0.438, p 0.103, n = 15",
      "gse_id": "GSE13485",
      "platform": "GPL7567",
      "rationale": "p=0.103 is above 0.05",
      "verdict": "inconclusive"
    },
    {
      "analysis_unit_id": "GSE13485_trial2",
      "caveats": [
        "Confidence is derived deterministically from the pre-registered rule and the adjusted p-value. It is not a calibrated posterior probability."
      ],
      "cohort": "Trial 2",
      "confidence": 0.07,
      "evidence_summary": "effect +0.032, p 0.93, n = 10",
      "gse_id": "GSE13485",
      "platform": "GPL7567",
      "rationale": "p=0.93 is above 0.05",
      "verdict": "inconclusive"
    }
  ]
}
```

</details>
