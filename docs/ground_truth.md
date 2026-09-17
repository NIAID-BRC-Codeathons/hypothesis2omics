# Ground truth for the YF-17D benchmark

Everything below is from the primary sources, checked against GEO directly.
Where something is unverified it says so.

## Why this document exists

Three of us are filling in dataset metadata by hand and hitting the same
`unknown` fields. Most of those unknowns are already answered in the papers
that produced the data. This is that answer set, so nobody guesses twice.

---

## The two benchmark hypotheses

Both come from Querec et al. 2009, Nat Immunol 10(1):116-125, PMID 19029902.
That paper is the source of GSE13485.

**H1 (CD8 arm).** Early expression of EIF2AK4 (GCN2) correlates with and
predicts the magnitude of the later YF-17D specific CD8+ T cell response.

**H2 (antibody arm).** Early expression of TNFRSF17 predicts the neutralizing
antibody response.

Querec's own words: a gene signature including complement C1qB and EIF2AK4
"correlated with and predicted YF-17D CD8+ T cell responses with up to 90%
accuracy in an independent, blinded trial", and a distinct signature including
TNFRSF17 "predicted the neutralizing antibody response with up to 100%
accuracy".

**Direction is positive.** Confirmed independently by Ravindran et al. 2014,
Science 343(6168):313-317, PMID 24310610, PMC4048998: YF-17D induced early
GCN2 expression in blood "strongly correlates with the magnitude of the later
CD8+ T cell response". That paper also supplies causality in mice, not just
correlation: GCN2 knockouts show reduced IFN-gamma+ CD8 in liver and lung, the
effect is dendritic cell intrinsic (bone marrow chimeras plus
GCN2flox/flox x CD11c-cre), and the mechanism is autophagy driving antigen
presentation.

So `direction = "up"` is defensible and pre-registerable. It is not being
inferred from the hypothesis text, it is taken from two published studies.

**Predictor timepoint is day 7.** Ravindran notes the human signature kinetics
peak at day 7, not day 3 and not day 15.

> Unverified: the proposal cites a "Rodrigues-Coffinet 2025 Science
> Immunology" paper. It did not turn up in any search. Somebody should confirm
> it exists before it goes in a report.

---

## GSE13485, what the record actually says

Pulled from the GEO accession page, not from a secondary source.

| Field | Value |
|---|---|
| Samples | 87 |
| Subjects | 25 |
| Tissue | PBMC |
| Platform | GPL7567, Human Genome U133 Plus 2.0 Custom CDF v9, 20,077 genes |
| Timepoints | pre-vaccination, day 1, day 3, day 7, day 21 |
| BioProject | PRJNA113999 |
| Parent | SubSeries of GSE13486 |
| Lab | Pulendran, Emory Vaccine Research Center |

Three things follow from that table.

**EIF2AK4 is measurable here.** The `unknown` in our metadata can be closed.
Querec identified EIF2AK4 from this exact array, so a Custom CDF v9 U133 Plus
2.0 carries it by construction. Same for TNFRSF17.

**There is no day 15.** GSE13485 samples at day 21. Any spec that names day 15
as the outcome timepoint needs correcting, or it will silently match nothing.

**Participant IDs are in the sample titles.** They look like
`Trial1 Subject ID 1901 Day 0`. So subject and timepoint are parseable from
GEO metadata alone, which means within-GEO participant linkage is solvable
without ImmPort. That closes a second `unknown`.

### The split that matters

GSE13485's own Overall design field: "This study is a combination of 2
separate YF-17D vaccination trials. Trial 1 has 15 subjects, and Trial 2 has
10 subjects. The two trials occurred 1 year apart and used different lots of
the vaccine."

The "independent, blinded trial" in Querec's abstract is Trial 2, and it is
inside GSE13485. So the faithful reproduction is fit on Trial 1, test on
Trial 2. One accession, two independent groups.

This breaks an assumption in `evidence_rules.independence_group()`, which
treats an accession as the smallest unit and cannot split inside one. Noted as
a real limitation, not a design choice. Fix is to allow an explicit
`group_override` per row.

### Sibling series

GSE13486 (SuperSeries, BioProject PRJNA110163) contains:

- GSE13484, 2 donors, in-vitro YF-17D stimulation of PBMC, platform GPL7566
- GSE13485, the 25 recipient vaccination time course

GSE13484 is an in-vitro stimulation experiment, not a vaccination trial. It
should not be counted as a second independent confirmation of an in-vivo
hypothesis. It is the same lab, same submission, same BioProject.

---

## The outcome problem, stated precisely

Expression is in GEO. The outcome is not.

GEO holds the transcriptome. It does not hold per participant CD8+ T cell
response magnitude or neutralizing antibody titer. Those live in ImmPort and
are surfaced by ImmuneSpace. That is why `outcome_data_accessible` is
`unknown` across SDY1264, SDY1529 and GSE13485 at once. It is one problem, not
three, and it is structural rather than a metadata gap.

Likely route: ImmuneSpace, `ImmuneSpaceR::getDataset("fcs_analyzed_result")`
for the CD8 arm and the neutralizing antibody assay table for H2, joined on
`participant_id`. Worth confirming with Yaphet, who has been posting the
ImmPort API docs.

Until somebody pulls actual participant level outcome values, every hypothesis
we ask correctly returns `inconclusive`. That is the pipeline working, but it
is not a result.

---

## Dataset notes

**SDY1529** is an antibody focused study. It will likely fail a CD8 outcome
criterion by design, which is correct behaviour rather than a bug. It is a
good fit for H2, the TNFRSF17 arm.

**GSE125921 and GSE136163** are both the GEO face of SDY1529. They are one
study, not two. Counting them separately would inflate independent group count
by one, which is exactly the failure mode the grouping code exists to catch.

**GSE13699** (Lausanne) is a genuinely separate programme and is the most
useful candidate for a second independent group on the CD8 arm.

---

## What this gives the pipeline

Pre-registerable, before any numbers exist:

```python
rule = DecisionRule(
    alpha=0.05,
    direction="up",              # from Querec 2009 and Ravindran 2014
    min_independent_groups=2,
)
```

Predictor: EIF2AK4 at day 7. Outcome: YF-17D specific CD8+ T cell response
magnitude. Expected sign: positive. Replication target: Trial 1 fit, Trial 2
test, plus GSE13699 as a separate programme.

And the honest negative control, unchanged: a hypothesis with no supporting
data must come back `inconclusive`, not `refutes`.

---

## Sources

- Querec TD, Akondy RS, Lee EK, Cao W, et al. Systems biology approach
  predicts immunogenicity of the yellow fever vaccine in humans. Nat Immunol
  2009;10(1):116-125. PMID 19029902.
- Ravindran R, Khan N, Nakaya HI, et al. Vaccine activation of the nutrient
  sensor GCN2 in dendritic cells enhances antigen presentation. Science
  2014;343(6168):313-317. PMID 24310610, PMC4048998.
- GEO GSE13485, GSE13486, GSE13484 accession records, read directly.
