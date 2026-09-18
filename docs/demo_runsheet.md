# Project 10 demo run-sheet

**Hypothesis2Omics: Agentic Multi-Omics Validation** | NIAID-BRC AI Codeathon 2.0

Every number below was produced by a command in this sheet, on the data committed to the repo. Nothing here is illustrative.

This is a script for the presenter, not something you execute. The fenced blocks are the commands to type; the prose around them is what to say.

---

## The one sentence

> We built a pipeline that tells you when you don't know.

If you only get thirty seconds, that is the sentence. Everything else supports it.

---

## Before you start

```powershell
cd $env:USERPROFILE\Documents\h2o-team
git status --short
git rev-parse --abbrev-ref HEAD
python evidence_rules\report.py --self-test
python evidence_rules\run_h2.py --repo . --exploratory | Select-Object -Last 4
```

Expect `report.py: 25 assertion groups passed` and the leave-one-out block from `run_h2.py`. If either fails, fix it now rather than discovering it live.

**Which branch.** Everything demoed here lives on `evidence-report` until PR #10 merges. Check with `git rev-parse --abbrev-ref HEAD` and stay there. Do not `git pull --ff-only origin main` mid-demo; if #10 has merged, switch to `main` beforehand and re-run the pre-flight.

**Step 2 writes into `docs/`.** Running it regenerates `docs/evidence_report_yf17d.md` and its `.json`, so `git status` goes dirty and the diff is the commit hash the report stamps. Harmless. Do not commit it during the demo, and `git checkout -- docs/` afterwards if you want a clean tree.

**Fallback if the live runs misbehave:** `docs/evidence_report_yf17d.md` is already committed. Open it and walk through it instead. Decide before you start whether you are running live or reading the committed report, and do not switch mid-demo.

---

## Step 1 - What goes in and what comes out (60 seconds, no commands)

A plain-English hypothesis goes in. A search over ImmPort and GEO, an eligibility judgment on every candidate, a real analysis on the eligible ones, and an evidence report come out.

Say the scale once, then move on:

- 4 ImmPort studies, 5 GEO series, **518 linked samples**, 6 study/platform analysis units
- 7 modules in the evidence layer, standard library only, each with its own self-tests
- A human reads and corrects the parsed hypothesis before anything proceeds. **Nothing runs unattended end to end, by design.**

Owners, if asked: Melanie (hypothesis to test spec), Yaphet (ImmPort search, MCP registration), Yijun (retrieval, parsing, normalization), Amar (eligibility, GEO matrix extraction), Archit and Slim (analysis execution), Rishi (decision rules, independence, synthesis, report).

---

## Step 2 - The benchmark, and the result nobody wants to report

```powershell
python evidence_rules\run_yf17d.py --repo . --report docs\evidence_report_yf17d.md
```

The hypothesis is Querec 2009: early EIF2AK4 expression correlates with the magnitude of the CD8+ T-cell response to YF-17D.

**What to say while it runs:** three defensible definitions of "early expression", same 25 subjects, and they do not agree.

| Specification | Trial 1 | Trial 2 | Verdict |
|---|---|---|---|
| **day 7 as measured (pre-registered)** | +0.438 p=0.103 | +0.032 p=0.930 | **inconclusive** |
| day 7 minus day 0 | +0.645 p=0.009 | +0.681 p=0.030 | supportive |
| day 7, baseline held constant | +0.545 p=0.044 | -0.045 p=0.907 | inconclusive |

**The line that matters:**

> The middle row is positive and significant in both trials. That reads as a replication, and it is the row we would have reported if we had picked after seeing the numbers. We didn't. The pre-registered specification is day 7 as measured, because that is what Querec's signature was built from, and it comes out inconclusive.

Then the why:

```
r(day 0, CD8):  Trial 1  -0.439 p=0.102     Trial 2  -0.826 p=0.003
```

> Baseline EIF2AK4 is itself correlated with the outcome. A day-7-minus-day-0 score inherits that by construction. Hold baseline constant and Trial 2 goes to -0.045.

---

## Step 3 - The report refuses to flatter itself (90 seconds)

Open `docs/evidence_report_yf17d.md`. Show three things, in this order:

**1. The rule, printed before the results.** Alpha, direction, minimum effect, replication threshold, all fixed during planning. A reader can check that what was promised is what was applied.

**2. The sensitivity table.** All three specifications side by side with their verdicts, and the sentence "the verdict is specification-sensitive" in those words. The generator requires exactly one specification marked primary with a stated reason and **raises if you don't give it one**. You cannot render a verdict without saying which specification counts.

**3. What was *not* analysed.** GSE13486, GSE125921, GSE13699, SDY1291, each with the stage it left at and why. The constructor rejects an entry missing a reason.

> A report listing four analysed datasets and nothing else implies four datasets were all there was.

---

## Step 4 - The system refusing to answer (the memorable part)

```powershell
python evidence_rules\run_h2.py --repo . --report docs\h2.md
```

This **refuses**, writes no file, and exits 2.

> Second hypothesis from the same paper: TNFRSF17 predicts the neutralizing antibody response. We have the data. We have the numbers. The pipeline will not report them, because nobody has confirmed which probe on this array is TNFRSF17.

Then the uncomfortable part, which you should say rather than wait to be asked:

> We found that the gene-to-probe pairing our own benchmark rests on is a hand-entered value in a config file. Nothing in the repository maps a probe to a gene. It travelled downstream into a column that looks like provenance. So the same gap applies to the result in step 2.

```powershell
python evidence_rules\run_h2.py --repo . --exploratory
```

| Specification | Trial 1 (n=15) | Trial 2 (n=10) |
|---|---|---|
| **day 7 as measured (pre-registered)** | -0.433 p=0.107 | -0.337 p=0.341 |
| day 0 baseline | +0.419 p=0.120 | +0.337 p=0.341 |
| day 7 minus day 0 | **-0.679 p=0.005** | **-0.655 p=0.040** |
| day 7, baseline held constant | -0.525 p=0.054 | -0.587 p=0.097 |

> The difference score is significant and negative in both arms. A naive read reports that we refuted Querec. It doesn't survive: baseline is *positively* correlated here, so the difference gets dragged negative by construction. Hold baseline constant, neither arm reaches alpha. Trial 2 fails leave-one-out at p=0.140.

**The payoff line:**

> Same artifact, caught twice, on two different hypotheses, in opposite directions. That is what freezing the rule first buys you.

---

## Step 5 - Close (30 seconds)

> The pipeline's honest output on this benchmark is: inconclusive, inconclusive, and a refusal. We think that is the right answer, and we think a system that produces it is harder to build than one that produces a result.

---

## Questions you should expect

**"So you didn't find anything."**
We reproduced the specification-sensitivity that makes this class of result fragile, and we built the thing that catches it. The finding is about the method. At n=15 and n=10 on a four-value dilution ladder, a definitive answer would have been the suspicious outcome.

**"Why is inconclusive a result?"**
Because the alternative is a false claim. A hypothesis whose predictor was never measured, or whose timepoint doesn't exist in the data, returns inconclusive rather than refutes. The question went unasked.

**"Is this just a wrapper around a t-test?"**
The statistics are ordinary. The refusals are the contribution: a frozen rule, independence counted in groups rather than rows, controls barred from evidence tables, estimands that will not pool, and a probe identity that gates report generation.

**"What about the genome-wide analysis?"**
Archit ran 12 specifications across two studies, GSE13485 and GSE13699. Zero probes survive FDR at q<0.05 in any of them, which is the correct result for a discovery screen at this sample size. It is a different estimand from our benchmark, antibody titer rather than CD8 response, and is reported as exploratory rather than folded into the verdict.

**"What's not done?"**
The eligible-to-analysis-ready step is only partly automated: `data/extract_series_matrix.py` produces the matrix, the design file is still prepared by hand, and limma needs both. The second cohort for H2 is blocked on a platform annotation. And no study in the antibody set has a usable baseline titer, so a fold-change antibody endpoint isn't available anywhere.

---

## Do not say

- "We replicated Querec." We did not, on the pre-registered specification.
- "We refuted Querec." We did not, and the specification that looks like it does is an artifact.
- Any number for TNFRSF17 without the words "the probe is not confirmed" attached.
- "The pipeline is fully automated." One stage is manual and the README says so.

---

## Timing

| Step | Minutes |
|---|---:|
| 1 What it does | 1.0 |
| 2 The benchmark | 1.5 |
| 3 The report | 1.5 |
| 4 The refusal | 2.0 |
| 5 Close | 0.5 |
| | **6.5** |

Leaves room in a 10-minute slot for one demo hiccup and two questions.
