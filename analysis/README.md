# analysis

Stage 6 of the pipeline: executing the analysis on an eligible dataset.

Owner: Archit Vasan. The code in this directory is his; this README documents
what it needs and what is not here yet.

---

## What is here

`run_galaxy_correlation.py` runs the YF-17D probe-versus-titer Spearman
correlation on Galaxy, using the `Feature-wise Correlation Tests` tool rather
than local scipy, and compares the result probe-by-probe against the local
reference tables.

```bash
export GALAXY_URL=https://usegalaxy.org
export GALAXY_API_KEY=<your key>          # environment only, never in a file
python analysis/run_galaxy_correlation.py "A__FC_d7-d0__rank"
```

It writes `galaxy_out/perfeature_*.tsv`, `galaxy_out/summary_*.tsv` and a
provenance JSON recording the tool id, every tool parameter, the Galaxy history
and job ids, and SHA256 of both inputs and the per-feature output.

The recorded cross-check: Galaxy job `bbd44e69cb8906b536b291633e51df67` agreed
with local scipy to about 5e-07 on rho, raw p and the BH-adjusted q. The
comparison checks all three deliberately, because agreeing on rho alone would
not establish that the two implementations adjust the same way.

That job id is not the one carried by the artifacts now committed under
`galaxy_out/`, which record `bbd44e69cb8906b5315143fa638bf584` and
`bbd44e69cb8906b58409527efaaaf5f4`. Either this line refers to a run whose
outputs are not in the repository, or it is stale. Left as written pending
confirmation from the author.

`verify_galaxy_agreement.py` re-derives that 5e-07 figure from the committed
tables rather than from the runner's own `comparison_vs_local` block:

```bash
python analysis/verify_galaxy_agreement.py
```

It checks each per-feature table against the `per_feature_sha256` recorded in
its provenance JSON, then recomputes feature-set identity, the maximum absolute
rho, p and q differences, the top-1000 overlap and the q<0.05 counts. It exits
non-zero on any failure and needs nothing installed. See
`galaxy_out/README.md` for the coverage it currently spans, which is two of the
twelve specifications.

---

## What it cannot do yet

**This script will not run from a fresh checkout.** It requires
`analysis/galaxy_in/matA_<key>.tsv` and `matB_<key>.tsv`, which are written by
`export_for_galaxy.py`. That script is not in this repository. Running
`run_galaxy_correlation.py` without it exits immediately with a `missing ...`
message naming the file it wanted.

It also requires `analysis/out/corr_<key>.tsv`, the local reference table, for
the comparison step. Without it the Galaxy job still runs and the provenance is
still written, but `comparison_vs_local` is null.

Needed to make this directory self-contained:

- `export_for_galaxy.py`, which builds the two input matrices
- `build_cohorts.py`, which produces the local reference tables
- `bioblend` added to the project dependencies; it is imported here and is not
  in `pyproject.toml`

---

## What runs where, and why it matters

From the script's own docstring, and worth repeating because it bounds what the
Galaxy integration demonstrates:

| Half | Does |
|---|---|
| Galaxy | the per-feature Spearman and the BH adjustment, that is, the arithmetic |
| local | subject-level collapse, the day-7-minus-day-0 fold change, within-arm rank de-confounding, the log2 = 7.0 floor filter |

The local half is not reproducible in Galaxy because no tool implements it, and
every consequential design decision lives there. This script makes the Galaxy
half reproducible and verifies it. It does not make the preparation
reproducible, and it should not be described as though it does.

One consequence worth stating plainly: the choices that made our benchmark
verdict specification-sensitive are all in the unreproducible half. Provenance
here is strongest where the decisions are fewest.

---

## The encoding, and its cost

The Galaxy tool correlates feature *f* in matrix A against feature *f* in
matrix B, pairing by name. To express "every probe against one titer vector",
matrix B is written with the same titer column repeated once per probe.
Feature-wise correlation then reduces exactly to `Spearman(probe_i, titer)`.

It works, and it is the right way to drive this tool. It also means matrix B
carries roughly 500,000 cells for a 20,077-probe run across 25 subjects, all
copies of one column, about 4 MB of redundancy per upload. Worth knowing before
running many of them.

---

## Reading the output honestly

The script prints this when nothing survives correction, and it is correct to:

> zero probes at FDR<0.05 -- agreement confirms the arithmetic, not a finding.
> A reproduced null is still null.

Across the twelve specifications run so far, over two studies, no probe reaches
`q < 0.05`. That is the expected result for a genome-wide screen at n = 13 to
25 and is reported as exploratory rather than folded into the benchmark
verdict, because the outcome is antibody titer while the primary benchmark
hypothesis concerns the CD8 response. Different estimand, reported separately.
