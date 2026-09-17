# File contract

What each stage reads and writes. Read this when a stage reports a missing input, or when
you need to know whether an artifact already exists before paying to regenerate it.

`h2o_preflight.py` checks this contract mechanically. Prefer running it over reading this
file to answer "what state am I in".

## The chain

```
STAGE                          WRITES
-----                          ------
1  parse_hypothesis            <outdir>/01_parsed.yaml
       |                            ^ HUMAN GATE. A person reads and corrects this.
2  build_test_spec             <outdir>/02_test_spec.json
                               <outdir>/03_search_spec.yaml
3  search_spec                 nothing (returns a list of accessions)
3b study_readiness_module      <out>/study_readiness.tsv
                               <out>/ready_studies.txt
                               <out>/excluded_studies.tsv
                               <out>/study_readiness.provenance.json
4  fetch_immport_studies       data/immport_cache/<SDY>/...
                               data/immport_cache/manifest.json
                               data/immport_cache/provenance_log.jsonl
5  parse_immport_studies       data/immport_cache/<SDY>/parsed/sample_manifest.tsv
                               data/immport_cache/parsed/sample_manifest.tsv      <-- KEY
                               data/immport_cache/parsed/batch_manifest.json
6  plan_geo_retrieval          data/geo_cache/plan/geo_download_plan.tsv          <-- KEY
                               data/geo_cache/plan/geo_download_plan.provenance.json
                               data/geo_cache/plan/gsm_resolution_cache.json
                               data/geo_cache/plan/gse_metadata_cache.json
7  fetch_planned_geo           data/geo_cache/<GSE>/...
                               data/geo_cache/manifest.json
                               data/geo_cache/provenance_log.jsonl
8  parse_geo_matrices          data/geo_cache/parsed/geo_matrix_parse_manifest.json  <-- KEY
                               data/geo_cache/parsed/<SDY>/<EXP>/<GSE>_<GPL>/
                                   expression.tsv.gz
                                   samples.tsv
                                   geo_sample_metadata_long.tsv
                                   provenance.json
9  validator_handoff_parse     <config.output_dir>/sample_manifest.tsv
                               <config.output_dir>/feature_expression.tsv
                               <config.output_dir>/quantitative_outcome.tsv
                               <config.output_dir>/validator_input_manifest.json
                               <config.output_dir>/*.provenance.json
10 handoff_adapter             <--output>.json   (criterion-level evidence)
11 eligibility_engine          stdout, and <--output>.json when given
12 analysis                    Galaxy limma outputs, or an in-process regression
13 synthesize                  a synthesis.schema.json document, in memory
```

The three paths marked KEY are the joints. If one is missing, every stage after it is
blocked, and that is what a "blocked, missing ..." line from the preflight means.

## Hazard 1 — the validator bundle exists twice

Two directories hold the same three filenames, and two different consumers read them:

| Directory | Written by | Read by |
|---|---|---|
| `data/validator_input/` | `data/validator_handoff_parse_module.py` (its config's `output_dir`) | `scientific_validator/handoff_adapter.py` (its `--input-dir` default) |
| `data/validator_immport/` + `data/validator_geo/` | committed to the repository | `evidence_rules/run_yf17d.py`, which hardcodes both paths |

`run_yf17d.py` does **not** read `data/validator_input/`. Check which directory a script
wants before pointing it anywhere.

Do not unify them without asking. They currently hold the same study, but the split is what
`run_yf17d.py` depends on, and merging them is a behavior change rather than a cleanup.

## Hazard 2 — `derive_search_spec` ignores your output directory

It always writes `03_search_spec.yaml` into the parent of the `02_test_spec.json` you passed.
`outdir` is not a parameter of that tool. `build_test_spec` does accept `outdir` and writes
both files there.

## Provenance

Every stage above writes a sibling `*.provenance.json` recording source, output hashes,
versions, counts, and timing — `AGENTS.md` makes this mandatory: *"Any data retrieval,
dataset selection, workflow execution, or eligibility decision must record: source/query
used, timestamp, tool/version, and inputs → outputs. No unlogged side effects."*

So an artifact with **no** provenance sibling was produced outside the pipeline. That is the
check `h2o_preflight.py` runs, and it is the mechanical signature of hand-rolled work.

Known exception, recorded so the check does not cry wolf: `data/extract_series_matrix.py`
writes `data/galaxy_file_input/<GSE>_intensities.csv` with no provenance file. That is a gap
in the script, not a licence to skip provenance elsewhere.

## Caches and what is committed

`data/immport_cache/` and `data/geo_cache/` are gitignored — a fresh clone has neither, and
the preflight will correctly report stages 4–8 as not run even though the validator bundles
downstream of them are present. Those bundles, `data/galaxy_file_input/`, and
`scientific_validator/output/` are committed so the benchmark can be re-run without
credentials.

Repeat fetches are cache-satisfied: `fetch_immport_studies`, `fetch_planned_geo` and
`plan_geo_retrieval` all pass `force=False` and return status `cached`. Only the CLIs take
`--force`. Re-running a fetch is cheap; re-running it with `--force` is not.
