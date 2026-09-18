# galaxy_out

Provenance, summary and per-feature artifacts written by
`run_galaxy_correlation.py`. See `analysis/README.md` for what the script does,
what it needs, and the boundary between the Galaxy half and the local half.

## Coverage

Two of the twelve specifications recorded in `analysis/out/correlation_provenance.json`
have been rerun on Galaxy:

| analysis_key | n subjects | Galaxy job id |
|---|---|---|
| `A__FC_d7-d0__rank` | 25 | `bbd44e69cb8906b5315143fa638bf584` |
| `A_d60arm__FC_d7-d0__log2titer` | 15 | `bbd44e69cb8906b58409527efaaaf5f4` |

The other ten have local results only.

Note: `analysis/README.md` cites job `bbd44e69cb8906b536b291633e51df67` for the
cross-check. That id matches neither job above, so it refers to a run whose
artifacts are not here.

## Checking the agreement yourself

```bash
python analysis/verify_galaxy_agreement.py
```

It verifies each per-feature table against the `per_feature_sha256` recorded in
its provenance JSON, then recomputes the comparison against the local reference
tables rather than reading the `comparison_vs_local` block that the runner
wrote. Exits non-zero if anything fails. Standard library only.

Expected output:

| analysis_key | feature sets | max abs rho diff | max abs p diff | top-1000 overlap |
|---|---|---|---|---|
| `A__FC_d7-d0__rank` | identical, 20,077 | 4.99831e-07 | 4.99378e-08 | 1000/1000 |
| `A_d60arm__FC_d7-d0__log2titer` | identical, 20,077 | 4.90143e-07 | 4.98952e-08 | 1000/1000 |

Zero features reach q < 0.05 in either implementation for either
specification, consistent with `n_q05: 0` across all twelve specifications in
`correlation_provenance.json`. Agreement confirms the arithmetic, not a
finding. A reproduced null is still null.

## Files

- `galaxy_run_<key>.json` — job and history ids, input SHA256s, tool id and all
  tool parameters, and the runner's own comparison against the local reference.
- `summary_<key>.tsv` — the Galaxy tool's run summary.
- `perfeature_<key>.tsv` — per-probe rho, p, BH-adjusted q and significance
  flag, 20,077 rows. Committed so the agreement check above can be re-run.

The matching local reference tables are `analysis/out/corr_<key>.tsv`. Only the
two with a Galaxy counterpart are committed; the other ten local tables are
summarized in `correlation_provenance.json` and have nothing to be checked
against.
