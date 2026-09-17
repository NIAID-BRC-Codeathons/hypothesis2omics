# Tool surface

Exact names, parameters, defaults, credentials, and the per-tool traps. Read the row for a
tool before calling it for the first time in a session.

MCP tool names are `<mcp.json key>_<tool>`, hyphens becoming underscores. `mcp/register.mjs`
sets the key from the manifest's `"name"`, defaulting to the directory name.

| Manifest directory | mcp.json key | Tool prefix |
|---|---|---|
| `mcp/parse-hypo/` | `hypothesis-parser` | `hypothesis_parser_` |
| `mcp/keyword-to-immport/` | `keyword-to-immport` | `keyword_to_immport_` |
| `mcp/data_normalizer/` | `hypothesis2omics` | `hypothesis2omics_` |

None of the three servers expose MCP **prompts** or **resources**. Tools only.

---

## `hypothesis_parser_*` — needs an LLM gateway

Credentials: `OPENAI_API_KEY` or `ARGO_USER`. Optional `OPENAI_BASE_URL`, `OPENAI_MODEL`,
`OPENAI_FALLBACK_MODEL`, `H2O_TIMEOUT_S`. Defaults resolve to Argo
(`https://apps.inside.anl.gov/argoapi/v1`, model `claudeopus5`), which needs the
Argonne-auth network.

> **Trap.** Setting `OPENAI_BASE_URL` away from Argo without also setting `OPENAI_MODEL`
> leaves the default model id `claudeopus5` — an Argo id that means nothing at another
> gateway. `h2o_preflight.py` fails on this.

### `parse_hypothesis` — **costs an LLM call**
`hypothesis_text: str`, `outdir: str`, `model: str = DEFAULT_MODEL`, `max_tokens: int = 200000`
→ `{path, parsed, warnings, yaml, provenance, next_step}`. Writes `<outdir>/01_parsed.yaml`.

> **Trap.** `max_tokens` defaults to `200000` as an MCP tool and `4000` in the
> `parse_hypothesis.py` CLI. They are not the same call.

**Stop after this.** The human gate is here.

### `build_test_spec` — **costs an LLM call plus an ImmPort read**
`parsed_path: str`, `outdir: str | None = None` (defaults to the parsed file's parent),
`model: str = DEFAULT_MODEL`, `hypothesis_id: str | None = None`, `max_tokens: int = 200000`
→ `{test_spec_path, search_spec_path, test_spec, search_spec, search_spec_yaml, provenance}`.
Writes `02_test_spec.json` and `03_search_spec.yaml`. Records the parsed file's sha256 prefix.

It reads ImmPort's facet vocabulary to ground the spec; on failure it warns to stderr and
proceeds ungrounded.

### `derive_search_spec` — free, deterministic
`test_spec_path: str`, `poc_format: bool = False` → `{path, search_spec, yaml}`.

> **Trap.** Always writes `03_search_spec.yaml` **beside the `02_test_spec.json` you passed**.
> There is no `outdir` parameter. Passing one has no effect because it is not accepted.

### `audit_search_terms` — free of charge, expensive in requests
`search_spec_path: str` → `{terms: [{group, term, required_group, hits, unique}], no_hits,
widens_required_set}`. Writes nothing.

> **Trap.** Issues **one ImmPort query per term**. On a spec with thirty terms that is thirty
> round trips. Use it once to audit a spec, not in a loop.

---

## `keyword_to_immport_*` — no credentials

Public endpoint `https://www.immport.org/data/query/api/search/study`. All three are
network reads that write nothing.

### `search_studies`
`keywords: list[str]` (must be non-empty), `limit: int = 10` (clamped 1–1000),
`exact_phrases: bool = False`, `filters: dict | None = None`
→ `{total, studies: [{accession, score, ...}]}`.

`sourceFields` is fixed at `study_accession,brief_title,condition_studied,research_focus`.

### `search_spec`
`spec: str | dict` — a YAML string or an object. Honors `search_terms` (group → terms),
`require`, `min_groups`, `filters`, `limit` (default 100), `per_term` (default 200).
`repositories` and `required_features` are parsed and ignored.

> **Trap.** Returns a **bare `list[str]`** of accessions, best-first. Not a dict like its
> siblings. Do not index it by key.

### `list_facet_values`
`facet: str | None = None` → `{facetName: [values...]}`, most-used first. `gender` is
renamed `sex`. Valid facets: `studyAccession subjectAccession programName conditionOrDisease
researchFocus clinicalTrial sex race ethnicity species minAge ageRange assayMethod
biosampleType hasAssessment hasLabTest searchFields`.

---

## `hypothesis2omics_*` — the data pipeline

Credentials: `IMMPORT_API_KEY` or `IMMPORT_API_KEY_FILE`, for the ImmPort fetch only.
Roots: `HYPOTHESIS2OMICS_DATA_ROOT` (default `<repo>/data`),
`HYPOTHESIS2OMICS_GEO_PARSED_ROOT` (default `<data_root>/geo_cache/parsed`).

> **Trap.** MCP callers cannot supply arbitrary paths to any of these. Every path derives
> from the configured data root. If you need a different location, set the environment
> variable before the server starts — not per call.

### Pure reads — free, no network, no writes

| Tool | Parameters |
|---|---|
| `list_analysis_units` | none → run status, parser version, counts, and every unit |
| `get_analysis_unit` | `unit_id: str` → provenance, sample fields, outputs with path/size/sha256 |
| `get_sample_metadata` | `unit_id`, `gsm_accessions: list \| None`, `attributes: list \| None`, `limit: int = 25` (1–100), `offset: int = 0` |
| `get_expression_info` | `unit_id: str` → id column, sample accessions, probe count, sha256 |

`unit_id` is `"{study_accession}/{experiment_accession}/{gse_accession}_{gpl_accession}"`.

### Downloads

- **`fetch_immport_studies(study_accessions: list[str], max_files_per_study: int | None = None)`**
  Accessions must match `^SDY\d+$`; duplicates rejected. Needs an ImmPort credential.
  Writes into `data/immport_cache/` plus `provenance_log.jsonl`.
- **`plan_geo_retrieval(batch_size: int = 100, timeout: int = 60)`**
  NCBI E-utilities. Writes `data/geo_cache/plan/`.
  > **Trap.** `data/geo_plan_module.py` honors `NCBI_API_KEY`, but
  > `mcp/data_normalizer/mcp-server.json` does not forward it. Through MCP you get the
  > unauthenticated rate limit. Use the CLI if that matters.
- **`fetch_planned_geo()`** — no parameters, no credential. Downloads only plan rows whose
  `download_recommendation == "download"`. To fetch an arbitrary GSE you must use
  `python data/geo_fetch_module.py GSE...`.

> **Trap.** All three pass `force=False` internally. A repeat call is satisfied from cache
> and returns status `cached`. Only the CLIs expose `--force`.

### Local writes, no network

`parse_immport_studies()`, `parse_geo_matrices()`, `inventory_downloaded_files()` — all take
no parameters and derive every path from the data root.

---

## CLI-only stages

No MCP tool exists for any of these. They are the back half of the pipeline.

```bash
# readiness gate: drop studies with no GEO linkage before downloading
python mcp/parse-hypo/study_readiness_module.py SDY1264 SDY1289 \
  --api-key-file data/immport_cache/<key>.json --output-dir data/readiness

# validator-input bundle
python data/validator_handoff_parse_module.py data/validator_handoff_config.json

# criterion-level evidence
python scientific_validator/handoff_adapter.py --input-dir data/validator_input \
  --study-accession SDY1264 --gene EIF2AK4 \
  --output scientific_validator/output/SDY1264_validator_input_evidence.json

# eligibility
python scientific_validator/eligibility_engine.py \
  --spec scientific_validator/examples/yf17d_test_spec.json \
  --dataset scientific_validator/output/SDY1264_validator_input_evidence.json

# Galaxy intensity matrices
python data/extract_series_matrix.py --input-dir /path/to/GEO/datasets

# the committed YF-17D benchmark (needs pandas + scipy)
python evidence_rules/run_yf17d.py --repo .
```

`handoff_adapter.py` defaults `--gene` to `EIF2AK4`; `--study-accession` and `--output` are
required. `eligibility_engine.py` always prints its verdict to stdout and writes a file only
when given `--output`.

## Library, not CLI

```python
from evidence_rules import DecisionRule          # frozen at construction
from synthesis import synthesize, validate_synthesis
from adapters import eligible_datasets, groups_from_verdicts, rule_from_spec

rule = DecisionRule(alpha=0.05, direction="up", min_independent_groups=2)
doc = synthesize(units, rule)
assert validate_synthesis(doc) == []
```

`evidence_rules/`, `synthesis.py` and `adapters.py` are standard library only and each runs
its own self-tests via `python <file>.py`.

## Smoke test

`python mcp/parse-hypo/hypothesis_mcp.py --check` runs stages 1 → 2 → 3 → audit against a
fixed hypothesis. It is a live check of the LLM gateway and the ImmPort endpoint in one call.
