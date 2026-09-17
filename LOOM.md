# hypothesis2omics

A biological hypothesis in; a search of public transcriptomic and immune-system
repositories, an eligibility judgment on each candidate dataset, an analysis on the
eligible ones, and an evidence report out. Every stage records provenance. Nothing here
runs unattended from end to end — a person reads and corrects the parsed hypothesis first.

## There is a skill for this repository

`.claude/skills/hypothesis2omics/SKILL.md` carries the full route: which tool or command
owns each of the thirteen stages, the file contract between them, which calls cost money
or download gigabytes, and the traps. Its two reference files go deeper on the tool
parameters and the path contract.

Reading it before calling a `hypothesis2omics_*`, `hypothesis_parser_*` or
`keyword_to_immport_*` tool saves reconstructing the pipeline from five separate READMEs.

## The preflight tells you where you are

```bash
uv run .claude/skills/hypothesis2omics/scripts/h2o_preflight.py --repo .
```

Standard library only, so it runs before `uv sync` has happened. It reports whether the
environment can run a stage, which stages are runnable now with their exact command lines,
and whether every artifact on disk carries provenance. Exit 0 proceed, 1 blocking, 2 usage.

## What is already installed

These exist as tools. There is no need to write equivalents, and the equivalents would not
write the provenance files the rest of the pipeline depends on.

| Tool | Does |
|---|---|
| `hypothesis_parser_parse_hypothesis` | hypothesis text → named components (LLM) |
| `hypothesis_parser_build_test_spec` | components → test spec + search spec (LLM) |
| `hypothesis_parser_derive_search_spec` | test spec → search spec, deterministic |
| `hypothesis_parser_audit_search_terms` | one ImmPort query per term; reports dead terms |
| `keyword_to_immport_search_studies` | keywords → ImmPort studies |
| `keyword_to_immport_search_spec` | a search spec → a bare list of accessions |
| `keyword_to_immport_list_facet_values` | the controlled vocabulary for filters |
| `hypothesis2omics_fetch_immport_studies` | downloads ImmPort Tab archives |
| `hypothesis2omics_parse_immport_studies` | Tab archives → sample manifests |
| `hypothesis2omics_plan_geo_retrieval` | resolves GSM → GSE via NCBI, writes a plan |
| `hypothesis2omics_fetch_planned_geo` | downloads the series the plan marks `download` |
| `hypothesis2omics_parse_geo_matrices` | series matrices → per-unit expression tables |
| `hypothesis2omics_list_analysis_units` | what was parsed, with counts |
| `hypothesis2omics_get_analysis_unit` | one unit's provenance and outputs |
| `hypothesis2omics_get_sample_metadata` | paged sample metadata for a unit |
| `hypothesis2omics_get_expression_info` | probe count, sample list, sha256 |

The back half has no MCP tools and is CLI only: the readiness gate
(`mcp/parse-hypo/study_readiness_module.py`), the validator bundle
(`data/validator_handoff_parse_module.py`), criterion evidence
(`scientific_validator/handoff_adapter.py`), eligibility
(`scientific_validator/eligibility_engine.py`), and synthesis
(`evidence_rules/`, imported as a library). The skill lists the exact command lines.

## Cost classes

- **LLM gateway**: the two `hypothesis_parser_*` parse and build steps.
- **Downloads**: `fetch_immport_studies`, `fetch_planned_geo`, `plan_geo_retrieval`, and
  `study_readiness_module.py`. All are cache-satisfied on a repeat call; only the CLIs
  expose `--force`.
- **Expensive read**: `audit_search_terms` issues one ImmPort request per term.
- Everything else is local or a free public read.

## Two places the paths are easy to get wrong

The validator bundle exists twice under different directory names.
`scientific_validator/handoff_adapter.py` reads `data/validator_input/`;
`evidence_rules/run_yf17d.py` hardcodes `data/validator_immport/` and `data/validator_geo/`.
Same three filenames, different directories, two different consumers.

`hypothesis_parser_derive_search_spec` has no output-directory parameter. It always writes
`03_search_spec.yaml` beside the `02_test_spec.json` it was given.

## Conventions this repository follows

`AGENTS.md` is the full set. The four that come up most:

- Existing infrastructure before new scripts — Galaxy, BRC Analytics, the published MCP
  tools. Differential expression, GEO matrix parsing and ImmPort search are all already
  solved here.
- Provenance is mandatory. Every stage records source, timestamp, version, inputs → outputs.
  An artifact without a sibling `*.provenance.json` came from outside the pipeline.
- Failures are surfaced, not worked around. A hand-rolled fallback after a tool error
  returns a plausible number with nothing to trace it to.
- Scientific thresholds — alpha, expected direction, effect floors, predictor timepoints,
  eligibility criteria — are the user's to set, not the agent's to infer.

## Setup notes

`uv` has to be on the PATH at launch; Orbit bundles its own, the CLI uses yours. The
`hypothesis2omics` server runs out of this repo's `pyproject.toml`, which pulls roughly
230 MB of scientific stack — installing that on first launch outruns the MCP startup window
and the server just looks broken, so `uv sync` once up front.

`IMMPORT_API_KEY` is needed only for the ImmPort fetch; `OPENAI_API_KEY` or `ARGO_USER`
only for the two parser steps. Everything else runs without credentials.

Loom's web and remote shell will not expose these servers: `web-mode-gate` is default-deny
with an allowlist covering `galaxy_`, `brc_analytics_`, `gtn_`, `notebook_` and
`skills_fetch`. CLI and Orbit are unaffected.

## Making the tool-first rule binding

Loom injects this file as *data, not instructions* — it is told the file may have shipped
with the folder rather than been written by you, so it carries no authority. The global
`~/.pi/agent/LOOM.md` is different: it rides in the cached system prompt under "follow it as
if the user had restated it at the top of this conversation."

`LOOM.global.md` in this directory is a short snippet for that file. **Append it, do not
copy over the file** — that path is global to every project you open, and a `cp` destroys
whatever standing instructions you already had:

```bash
cat LOOM.global.md >> ~/.pi/agent/LOOM.md     # or: loom, then /instructions init, then paste
```
