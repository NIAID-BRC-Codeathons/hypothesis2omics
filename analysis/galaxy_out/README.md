# galaxy_out

Provenance and summary artifacts written by `run_galaxy_correlation.py`.
See `analysis/README.md` for what the script does, what it needs, and the
boundary between the Galaxy half and the local half.

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

## Agreement with the local implementation

Recomputed from the Galaxy `perfeature_*.tsv` tables against the local
`out/corr_*.tsv` tables, rather than read from the provenance JSON's own
`comparison_vs_local` block:

| analysis_key | feature sets | max abs rho diff | max abs p diff | top-1000 overlap |
|---|---|---|---|---|
| `A__FC_d7-d0__rank` | identical, 20,077 | 4.99831e-07 | 4.99378e-08 | 1000/1000 |
| `A_d60arm__FC_d7-d0__log2titer` | identical, 20,077 | 4.90143e-07 | 4.98952e-08 | 1000/1000 |

Zero features reach q < 0.05 in either implementation for either
specification, consistent with `n_q05: 0` across all twelve specifications in
`correlation_provenance.json`. Agreement confirms the arithmetic, not a finding.

## Files

- `galaxy_run_<key>.json` — job and history ids, input SHA256s, tool id and all
  tool parameters, and the script's own comparison against the local reference.
- `summary_<key>.tsv` — the Galaxy tool's run summary.

The `perfeature_<key>.tsv` tables (about 1.6 MB each) are not committed. They
are derived and regenerable from the recorded job ids.
