---
name: hypothesis2omics
description: Drive the hypothesis2omics pipeline — a plain-English biological hypothesis through ImmPort and GEO retrieval, normalization, scientific eligibility, analysis, and evidence synthesis — using the MCP tools and CLI scripts this repository already ships, rather than writing replacement code. Covers (1) which tool or command owns each of the thirteen stages, (2) the file contract between stages and the two places it is easy to get wrong, (3) which calls cost money, download gigabytes, or are free reads, and (4) what must never be hand-rolled. Use whenever working in the hypothesis2omics repository or with its tools — parsing a hypothesis, searching ImmPort, fetching or parsing GEO series, building a validator-input bundle, assessing eligibility, or synthesizing a verdict. Trigger even when the request names no tool: "find datasets for this hypothesis", "is this study usable", "run the YF-17D benchmark", "why is the validator saying uncertain", and "extract the expression matrix" all qualify. It drives existing tools and does not choose scientific thresholds.
---

# hypothesis2omics

## What this produces

A route through the pipeline, and a preflight verdict on whether the next stage can run.
It does not produce a scientific result — the tools do that.

## First action, always

```bash
uv run "$SKILL_DIR/scripts/h2o_preflight.py" --repo .
```

`$SKILL_DIR` is the directory holding this file, known from the skill's location. Do not
search the filesystem for it. Without `uv`: `python3 "$SKILL_DIR/scripts/h2o_preflight.py"`
— the script is standard library only, on purpose, so it runs before `uv sync` has happened.

It prints three sections and the stages that are runnable now.

- **Exit 0** — proceed.
- **Exit 1** — something named is blocking. Fix what it names, then re-run. A non-zero exit
  is an instruction to go fix, not a reason to stop and not a reason to work around it.
- **Exit 2** — you pointed `--repo` at something that is not this checkout.

Re-run it after any stage that writes files. It is the cheapest way to find out that a
stage silently did nothing.

## Do not write code that a tool already runs

This is the failure this skill exists for. In testing, an agent with the tools present in
its list wrote its own code instead of calling them. The code ran, produced a table, and
looked like success — which is why it has to be caught deliberately.

`AGENTS.md` already states the rule: *"Do not reimplement functionality (e.g., differential
expression, groupby, ontology lookups) that an existing library or service already provides
correctly."* Made specific, keyed to what you are about to type:

| About to write | Call instead |
|---|---|
| `pd.read_csv` on a `*_series_matrix.txt` or `.txt.gz` | `hypothesis2omics_parse_geo_matrices()`; `python data/extract_series_matrix.py --input-dir DIR` for Galaxy CSVs |
| An HTTP request to `immport.org` | `keyword_to_immport_search_studies` / `keyword_to_immport_search_spec` |
| `wget` / `curl` / `urllib` against `ftp.ncbi.nlm.nih.gov/geo` | `hypothesis2omics_fetch_planned_geo()` |
| An E-utilities query to resolve GSM → GSE | `hypothesis2omics_plan_geo_retrieval()` |
| Inline differential expression | Galaxy limma. The result enters as an `AnalysisResult` (`evidence_rules/analysis_result.py`) |
| A hand-rolled verdict from a p-value | `evidence_rules.DecisionRule.score(row)`, then `synthesis.synthesize(units, rule)` |
| A prose judgment that a dataset is eligible | `python scientific_validator/eligibility_engine.py --spec ... --dataset ...` |
| A chosen alpha, direction, effect floor, or predictor timepoint | Nothing. `AGENTS.md` reserves these for explicit user sign-off. Ask. |

And the one that produced the observed failure:

> *"The tool errored, so I did it manually."*

`AGENTS.md`: **fail loudly, not silently.** A hand-rolled fallback is the silent failure —
it returns a plausible number with no provenance attached, and nothing downstream can tell
it apart from a real one. Report the tool error and stop.

**How this gets caught.** Every scripted stage writes a sibling `*.provenance.json`.
Hand-rolled code does not. `h2o_preflight.py` fails on any stage output that has no
provenance sibling. If you see that finding, an artifact was produced outside the pipeline
— find out which, and regenerate it with the stage that owns it.

## The human gate

`README.md`: *"A human reads and corrects the parsed hypothesis before anything proceeds.
Nothing in this pipeline runs unattended from end to end, by design."*

After stage 1, stop. Show `01_parsed.yaml` and wait. Do not call stage 2 on a parsed file a
person has not seen. In particular, `directionality` is never inferred from hypothesis text —
"associated with" states no sign, and guessing one is the failure the whole evidence layer
exists to prevent.

## The stages

Cost classes: **$** calls an LLM gateway · **↓** downloads · **·** free local or read-only.

| # | Stage | Call | Cost |
|---|---|---|---|
| 1 | hypothesis → parsed components | `hypothesis_parser_parse_hypothesis(hypothesis_text, outdir)` | $ |
| — | **human gate on `01_parsed.yaml`** | stop and show it | |
| 2 | parsed → test spec + search spec | `hypothesis_parser_build_test_spec(parsed_path)` | $ |
| 3 | search spec → ImmPort accessions | `keyword_to_immport_search_spec(spec)` | · |
| 4 | fetch ImmPort studies | `hypothesis2omics_fetch_immport_studies(study_accessions)` | ↓ |
| 5 | parse ImmPort Tab archives | `hypothesis2omics_parse_immport_studies()` | · |
| 6 | plan GEO retrieval | `hypothesis2omics_plan_geo_retrieval()` | ↓ |
| 7 | fetch planned GEO series | `hypothesis2omics_fetch_planned_geo()` | ↓ |
| 8 | parse GEO series matrices | `hypothesis2omics_parse_geo_matrices()` | · |
| 9 | build validator-input bundle | `python data/validator_handoff_parse_module.py data/validator_handoff_config.json` | · |
| 10 | extract criterion evidence | `python scientific_validator/handoff_adapter.py --input-dir data/validator_input --study-accession SDY... --gene ... --output ...` | · |
| 11 | assess eligibility | `python scientific_validator/eligibility_engine.py --spec ... --dataset ...` | · |
| 12 | analysis | Galaxy limma, or a per-subject regression | · |
| 13 | synthesize a verdict | `evidence_rules.synthesize(units, rule)` + `validate_synthesis(doc)` | · |

Stages 9–13 have **no MCP tools**. They are Bash and Python imports. Do not look for a
`hypothesis2omics_*` tool for them — there isn't one, and inventing a replacement is the
exact failure above.

There is an optional readiness gate between 3 and 4 that drops studies with no GEO linkage
before you download anything:

```bash
python mcp/parse-hypo/study_readiness_module.py SDY1264 SDY1289 \
  --api-key-file data/immport_cache/<key>.json --output-dir data/readiness
```

## Two contract hazards

**The validator bundle exists twice.** `handoff_adapter.py` reads `data/validator_input/`.
`evidence_rules/run_yf17d.py` hardcodes `data/validator_immport/` and `data/validator_geo/`.
Same three filenames, different directories. Check which one a script wants before pointing
it anywhere, and do not repoint either script to unify them without asking — that is a
behavior change, not a cleanup.

**`derive_search_spec` ignores your output directory.** It always writes `03_search_spec.yaml`
beside the `02_test_spec.json` you handed it. `outdir` is not one of its parameters.

## Reference files

Read on demand, not up front.

| File | When |
|---|---|
| `references/tool-surface.md` | before calling a tool you have not called in this session — exact names, parameters, defaults, credentials, and the per-tool traps |
| `references/file-contract.md` | when a stage says its input is missing, or when you need to know what a stage writes |

## Scope

This skill routes work to existing tools. It does not choose alpha, expected direction,
effect floors, predictor timepoints, eligibility criteria, or cohort groupings. `AGENTS.md`
reserves those for explicit user sign-off: *"Do not silently decide scientific thresholds,
calculation rules, schemas … these require explicit user sign-off, not agent judgment."*
When one is missing, ask. Do not pick a defensible default and proceed.

## Available scripts

- `scripts/h2o_preflight.py` — environment, pipeline state, and provenance check. Run it
  first and after any stage that writes files. `uv run "$SKILL_DIR/scripts/h2o_preflight.py"
  --repo .`, or `--json` for a machine-readable report. Exit 0 proceed, 1 blocking, 2 usage.
