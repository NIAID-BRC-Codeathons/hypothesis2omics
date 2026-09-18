# parse-hypo

A plain-English biological hypothesis in, a ranked list of ImmPort study accessions
out, with a reviewable artifact at every step.

```
"GCN2/EIF2AK4 activity is associated with the magnitude of the
 CD8+ T-cell response following YF-17D vaccination."
        |
        |  STEP 1  (LLM)    decompose into named components
        v
  01_parsed.yaml        intervention / predictor / cell_type / outcome /
        |               relationship_type / directionality + warnings
        |
        |  <-- THE HUMAN GATE: a person reads this file and fixes what
        |      the parser got wrong. Nothing proceeds on its own.
        v
        |  STEP 2  (LLM)    expand each component into synonyms, candidate
        |                   measurements, and checkable eligibility criteria
        v
  02_test_spec.json     also consumed by scientific_validator/eligibility_engine.py
        |
        |  STEP 3  (deterministic, no model)   collect per-component search_terms
        v
  03_search_spec.yaml   consumed by keyword-to-immport's search_spec
        |
        v
  ["SDY1264", "SDY1289", "SDY1294", ...]
        |
        |  STEP 4  (ImmPort, no model)   per-study GEO-readiness gate
        v
  study_readiness.tsv
        |         \
        v          v
  ready_studies.txt   excluded_studies.tsv
  -> data retrieval    -> kept for later review
     (immport_fetch_module.py)
```

Search terms are never invented in step 3. They are lifted from the component that
owns them in step 2, so every term traces back to a named part of the specification.

## Install and run

```sh
pip install mcp openai httpx pyyaml
export OPENAI_API_KEY=...             # see the credential note below
```

Steps 1 and 2 call an **OpenAI-compatible gateway**, configured the way any OpenAI
SDK user would expect:

| variable | default | |
|---|---|---|
| `OPENAI_API_KEY` | `ac.jdoe` | key, or an Argo username. `ARGO_USER` still works as a fallback. |
| `OPENAI_BASE_URL` | Argo | `OPENAI_URL` is accepted as an alias. |
| `OPENAI_MODEL` | `claudeopus5` | **Set this whenever you change the base URL** — the default is an Argo model id. |
| `OPENAI_FALLBACK_MODEL` | `claudesonnet5`, off once `OPENAI_MODEL` is set | second try when the first model returns nothing. |

The default is Argo, so an existing `ARGO_USER` setup needs no change. Against
the real thing:

```sh
export OPENAI_API_KEY=sk-...
export OPENAI_BASE_URL=https://api.openai.com/v1
export OPENAI_MODEL=gpt-4o
```

Two front doors onto the same code.

**CLI** — the reference implementation:

```sh
python3 parse_hypothesis.py --hypothesis "..." --outdir runs/HYP001
# read runs/HYP001/01_parsed.yaml, correct anything wrong
python3 build_test_spec.py runs/HYP001/01_parsed.yaml --audit --search
```

**MCP** — steps 1-4 as tools, for a client that would rather call them than shell
out:

```sh
python3 hypothesis_mcp.py --check     # live smoke test, no client needed
python3 hypothesis_mcp.py             # serve over stdio
```

ImmPort search itself lives in [`../keyword-to-immport/`](../keyword-to-immport/),
a separate server registered alongside this one, so a client can go from hypothesis
to accessions without leaving the session. Steps 2-3 import that module directly by
path -- there is one copy of the facet list and the spec search, not two.

Step 1 needs only `openai` and `pyyaml`. Steps 2-3 additionally need `httpx` and
`mcp`, because they import `keyword-to-immport/server.py`.

## The five files

```
h2o_common.py   ../keyword-to-immport/server.py   <- leaves: import nothing local
      |          |    |
parse_hypothesis.py   |            <- step 1, no ImmPort dependency at all
      |     |         |
   build_test_spec.py              <- steps 2-3
      |     |
study_readiness_module.py          <- step 4
      |     |
  hypothesis_mcp.py                <- MCP front door, wraps steps 1-4
```

| file | role |
|---|---|
| `h2o_common.py` | Gateway calls, JSON extraction, YAML scalars, provenance. Repository-agnostic. |
| `parse_hypothesis.py` | Step 1. Runs and tests with ImmPort unreachable. |
| `build_test_spec.py` | Steps 2 and 3. Grounds filters in ImmPort's controlled vocabulary. |
| `study_readiness_module.py` | Step 4. Gates candidates on sample-level GEO links. |
| `hypothesis_mcp.py` | MCP front door for steps 1-4; imports the modules above. |

That `parse_hypothesis.py` has no ImmPort dependency is a property worth keeping. It
means step 1 can be developed, tested and run without the repository being reachable,
and `hypothesis_mcp.py` defers its `build_test_spec` import to first use so a client
that only ever parses never loads the ImmPort layer either.

### h2o_common.py

Shared machinery. Nothing here knows what a hypothesis is.

| | |
|---|---|
| `client()` | An OpenAI client for the configured gateway. Imported lazily so `--help` works offline. |
| `extract_json(text)` | Pull one JSON object out of a model response, fenced or not. |
| `EmptyCompletion` | Raised when the gateway returns success with no content. Carries `completion_tokens`, which separates the two failure modes below. |
| `_complete(...)` | One chat completion. Raises `EmptyCompletion` rather than returning `""`. |
| `_complete_json(...)` | `_complete` plus JSON extraction, retried once with a doubled budget. |
| `_ask(...)` | The workhorse: one call, a model fallback on an empty completion, and one repair round-trip if the validator objects. Returns the object and what it took to get it. |
| `_scalar` / `_header` | YAML rendering: quote only when a bare word would misparse. |
| `_provenance(...)` | The provenance block. Records the model that *actually* answered, never the one requested. |

**Two empty-completion failure modes, and they need opposite fixes.**
`completion_tokens == 0` means the model generated nothing at all — a larger budget
will not help, and a different model is the only way forward. Non-zero means internal
reasoning ate the text budget, and a larger `max_tokens` fixes it. `_complete_json`
grows the budget for the second; `_ask` falls back to another model for the first.

### parse_hypothesis.py

| | |
|---|---|
| `PARSE_SYSTEM` | The step 1 prompt. Faithful and literal over creative — a later step does the expanding. |
| `validate_parsed(obj)` | Returns a list of problems. Empty means valid. |
| `parsed_to_yaml(...)` | Renders `01_parsed.yaml`, carrying `hypothesis_text` so step 2 need not reconstruct it. |
| `load_parsed(path)` | Reads the file back, hand-edits and all. Raises `ValueError` **listing every problem** if an edit broke the schema, so a bad edit stops here rather than surfacing three steps downstream. |
| `parse_hypothesis(...)` | The step itself. Returns the components and a provenance block. |

`warnings` is the parser naming what it inferred rather than read — an inferred
`cell_type`, an outcome with no measurable quantity. Nothing downstream branches on
it. It exists for whoever is standing at the gate.

### build_test_spec.py

| | |
|---|---|
| `_vocabulary()` | Fetches ImmPort's facet values so the model can only filter on strings that exist. |
| `_vocabulary_block(vocab)` | Renders them into the prompt. **This is the pattern the ontology change should follow.** |
| `SPEC_SYSTEM` / `SPEC_SCHEMA` | The step 2 prompt and its output shape. |
| `validate_spec(obj, vocab)` | Checks search groups, facet values against the live vocabulary, criterion ids, and that every scientific criterion has both a `description` and an `exclusion`. |
| `add_derived_fields(spec)` | Computes `inclusion_criteria`, `exclusion_criteria` and `primary_statistical_test.example_model` from fields the model already wrote, so the two renderings cannot drift. |
| `to_search_spec(spec)` | **Step 3.** Deterministic: collect `search_terms` per component, dedupe, keep order. |
| `search_spec_to_yaml(...)` | Renders `03_search_spec.yaml`. |
| `audit_terms(spec)` | One ImmPort query per term. `hits` is what a term returns alone; `unique` is what only that term found. |

`02_test_spec.json` is a superset of two contracts written independently: Step 2 of
`Example_input_output.md`, and what `eligibility_engine.py` reads. Both are in the
file, and the document's keys are *derived* from the engine's rather than asked for
twice.

### hypothesis_mcp.py

| tool | cost |
|---|---|
| `parse_hypothesis(hypothesis_text, outdir)` | 1 gateway call |
| `build_test_spec(parsed_path, outdir=None)` | 1 gateway call + an ImmPort vocabulary read |
| `derive_search_spec(test_spec_path)` | none — re-runs step 3 after a hand edit |
| `audit_search_terms(search_spec_path)` | 1 ImmPort query per term |
| `check_study_readiness(...)` | at most 1 cached Tab ZIP download per study |

Every tool writes its numbered artifact **and** returns the content, so an MCP run
leaves the same trail on disk as a CLI run. Errors raise: an empty completion, a
broken hand edit, or a spec the validator rejects surfaces as a failed tool call,
never as a plausible-looking empty result.

`H2O_TIMEOUT_S` caps a single gateway attempt (default 100s; the SDK retries, so a hard
hang costs a small multiple).

**On the gate.** `build_test_spec` takes a filesystem path, never an inline object,
so `01_parsed.yaml` has to exist before step 2 can run. That is a speed bump, not a
fence — an agent can write the file and call step 2 in the same turn. Provenance
records the input's sha256 and makes no claim that anybody read it. Evidence of
review has to come from outside this package.

## Step 4: study_readiness_module.py

Step 3 returns accessions ImmPort will hand back for a keyword match; it says
nothing about whether a study has the specific linkage data this pipeline needs
downstream. `study_readiness_module.py` checks, per accession, before the full
(multi-GB) retrieval runs:

- Downloads only the study's small `*_Tab.zip` release file (a few MB, not the full
  archive) and looks inside it.
- **`ready`** = it contains `expsample_public_repository.txt` with at least one row
  where `REPOSITORY_NAME == GEO` — i.e. the study has sample-level data actually
  linked to a GEO accession, not just a topical keyword match.
- `study_link.txt`, where present, is recorded (its links to ClinicalTrials.gov, a
  publication, GEO, etc.) but never gates the verdict: it is a study-level table of
  links to *any* external resource, and a study can have GEO-linked sample data
  without one (SDY63) or have one pointing somewhere else entirely, like
  ClinicalTrials.gov, while having no GEO-linked samples at all (SDY1479).

Measured 2026-09-17 against the 7 accessions `03_search_spec.yaml` returns for the
YF-17D hypothesis: the 4 hand-curated candidates (SDY1264, SDY1289, SDY1294,
SDY1529) plus SDY1291 came back ready; SDY15 and SDY271 — the 2 keyword-only hits
noted elsewhere as having no recall denominator — came back not ready, with no
`expsample_public_repository.txt` in either Tab archive at all. That's a partial
answer to that open question: those two extras are not usable by this pipeline
regardless of topical relevance.

```sh
python3 study_readiness_module.py SDY63 SDY1479 \
    --api-key-file /path/to/immport-key.json \
    --output-dir runs/READINESS01
```

Requires its own ImmPort API key (`--api-key-file`, with `browse` scope) — separate
from the LLM gateway credential steps 1-2 use, and from `keyword-to-immport`'s
public, keyless search API. Imports `../../data/immport_fetch_module.py` by path,
the same way `build_test_spec.py` imports `keyword-to-immport/server.py`, so there
is one source of truth for how ImmPort's manifest and download endpoints are
called.

Outputs: `study_readiness.tsv` (one row per study checked), `ready_studies.txt`
(plain accession list, the intended input to
`immport_fetch_module.fetch_immport_datasets`), `excluded_studies.tsv` (not-ready
studies with a reason, kept rather than discarded), and a provenance JSON.

The same operation is exposed by `hypothesis_mcp.py` as
`check_study_readiness(sdy_ids, outdir, cache_dir=None)`. The MCP tool resolves
credentials from `IMMPORT_API_KEY_FILE` or `IMMPORT_API_KEY`; it does not accept a
credential argument. The standalone CLI remains available for direct use with
`--api-key-file`.

## Rules that are load-bearing

Each was established by measuring the live services, not by reading documentation.
Removing one reintroduces a failure that has already happened.

- **`temperature=1`, stated explicitly.** Argo rejects every other value for
  `claudesonnet5` and the gpt5 family, and the SDK injects a default when the
  argument is absent. Omitting it is not the same as not sending it.
- **Argo model IDs are Argo's, not the vendor's** — `claudeopus5`, not
  `claude-opus-5`. `GET /v1/models` is the only source of truth. This is why
  `OPENAI_MODEL` has to be set whenever `OPENAI_BASE_URL` points somewhere else:
  the default id is meaningless off Argo.
- **ImmPort ANDs every word of a query and has no OR operator.** The pipeline runs
  one query per term and unions the results, which is why terms are grouped.
- **Never decompose a compound identifier.** Measured 2026-09-16: `17D` alone
  returned 14 studies of which 9 were unrelated (influenza, HIV, kidney transplant),
  because short fragments collide with arm and protocol codes. `YF-17D`, `YF17D`,
  `YFV` and `Stamaril` were all clean.
- **ImmPort indexes study titles, conditions and research focus — not molecular
  content.** Measured 2026-09-16: `GCN2`, `EIF2AK4`, `ATF4` and every other gene
  symbol returned zero studies. Predictor terms are still populated (they cost
  nothing and document intent) but `predictor` must never appear in
  `search.require`, which `validate_spec` enforces.
- **On Argo the key is an identifier, not a secret** — but it is a credential and
  it lands in shell history. Keep it in a file you `source`, not in `~/.bashrc` and not
  in anything committed. The network password is network-only and belongs in no
  script.

## Planned change: ontology-grounded synonym generation

**Status: specified, not implemented. The open decisions below need sign-off before
any code is written.**

### The problem

Step 2 currently asks the model to generate `search_terms` from its own knowledge.
That is the weakest link in the pipeline, for three reasons:

1. **It generates identifiers.** Asking a model for gene aliases and vaccine
   synonyms is precisely the situation where a fabricated identifier is plausible
   and unfalsifiable at a glance. A wrong alias does not error — it silently returns
   the wrong studies, or none.
2. **It is not reproducible.** Identical input produces a slightly different term
   list run to run, which changes the union and reorders the results. Observed: the
   same hypothesis returned the same seven accessions across runs, but not always in
   the same rank order.
3. **No term carries provenance.** AGENTS.md requires source, version and
   inputs-to-outputs for everything else in this pipeline. Search terms are the one
   thing that currently arrive from nowhere.

### The change

Ground `search_terms` in an ontology the same way `search.filters` are already
grounded in ImmPort's controlled vocabulary. `_vocabulary()` and
`_vocabulary_block()` are the existing pattern: fetch the authority, inject it into
the prompt, validate the output against it. The new step mirrors them.

```
parsed component  ->  ontology lookup  ->  candidate synonyms (with accessions)
                                                |
                                                v
                                   injected into SPEC_SYSTEM as a
                                   grounded candidate list, exactly
                                   as _vocabulary_block does today
```

Proposed shape, in `build_test_spec.py`:

| new | mirrors |
|---|---|
| `_synonyms(parsed) -> dict[str, list[dict]]` | `_vocabulary()` |
| `_synonym_block(terms) -> str` | `_vocabulary_block()` |
| a `source` check in `validate_spec` | the facet-value check already there |

**Three ways the terms could enter, and why the recommendation is the third.**

- **A. Replace.** Ontology lookup produces `search_terms` outright; no model
  involvement. Fully reproducible and fully traceable, but ImmPort indexes free text
  written by submitters, and ontology-correct labels are not always what a submitter
  typed.
- **B. Ground.** Retrieved synonyms go into the prompt as candidates and the model
  selects among them. Constrains fabrication without losing judgment, but the output
  is still a model's choice and still varies.
- **C. Union and label (recommended).** Ontology terms and model terms both go in,
  each tagged with its source. Then use `audit_search_terms` — which already
  measures exactly this — to find out which source earns its terms, and decide
  between A and B from data rather than from argument.

C is recommended because this codebase already owns the instrument to settle the
question, and neither A nor B should be chosen on intuition.

### Schema impact

`02_test_spec.json` carries structured terms with provenance:

```json
"intervention": {
  "name": "YF-17D vaccination",
  "search_terms": [
    {"term": "YF-17D",   "source": "<ontology>:<accession>", "version": "<release>"},
    {"term": "Stamaril", "source": "llm"}
  ]
}
```

(Placeholders deliberately: the authority is decision 1 below, and inventing a
plausible-looking accession here is the exact failure this section exists to
prevent.)

`03_search_spec.yaml` stays **flat lists of strings**, so `keyword-to-immport`'s `search_spec`
is unchanged. Step 3 flattens; provenance lives in step 2's artifact. This keeps the
search layer untouched by the change.

### Candidate authorities

Starting point for evaluation, **not a measured result** — none of these have been
called from this project, and the per-component mapping is a proposal:

| component | candidate | notes |
|---|---|---|
| `predictor` (genes) | HGNC, NCBI Gene | Authorities for human gene symbols and recorded aliases. |
| `intervention` (vaccines) | MeSH entry terms, Vaccine Ontology (VO) | MeSH carries trade and colloquial names; VO to be verified. |
| `outcome`, `cell_type` | Cell Ontology (CL), NCI Thesaurus | |
| `assay` | ImmPort's own `assayMethod` facet | Already fetched by `_vocabulary()`; reuse it. |

Access paths to compare: EBI OLS (no key), BioPortal (free key required), UMLS
(license required), or pinned local OBO/OWL files.

### Open decisions — these need sign-off

1. **Which authority per component.** A scientific choice, not an engineering one.
2. **Live API or pinned snapshot.** AGENTS.md puts reproducibility first, which
   argues for a pinned, versioned local file plus a cache. A live API is easier and
   makes a run un-reproducible.
3. **Failure behavior.** When a lookup returns nothing — fall back to the model,
   emit no term, or fail the step?
4. **Whether to ground `predictor` at all.** Gene symbols return zero ImmPort
   studies. Those terms are documentation and ranking only, so grounding them is
   correctness for its own sake. That may still be worth it; it should be a decision.
5. **Whether step 1 also gets ontology normalization.** Out of scope as written here.

### Acceptance criteria for the change

- Every term in `02_test_spec.json` traces to an ontology accession and version, or
  is explicitly tagged `source: llm`.
- Identical input yields an identical ontology-derived term set.
- `audit_search_terms` shows no loss of relevant studies against the current
  baseline for the YF-17D hypothesis (7 accessions: SDY1264, SDY1289, SDY1294,
  SDY1291, SDY1529, SDY271, SDY15).
- The ImmPort layer and the `03_search_spec.yaml` schema are unchanged.
- Lookups are cached; the pipeline does not re-query an authority per run.

## Not included here

- `validate/test_spec_shape.py` — offline gate that checks a `02_test_spec.json`
  against both contracts it has to satisfy. Worth bringing along as dev tooling.
- `validate_parse.py`, `validate_judge_selftest.py` — the LLM-judge harness for
  step 1, where the judge must quote a verbatim span and the span is checked in code.

## Known limitations

- **Criterion ids are model-generated and unstable.** The eligibility engine keys
  its assessment on them, and the hand-written reference spec overlapped ours on
  1 of 8. Until ids come from a fixed list, engine output is largely `unknown`.
  Team decision, unresolved.
- **Step 1 free-text fields vary run to run.** Categorical fields
  (`relationship_type`, `directionality`) are stable; `intervention` took four
  distinct values on identical input. This is the same root cause the ontology
  change addresses.
- **The CLI does not record the parsed file's sha256; the MCP does.** Two front
  doors, slightly different provenance. Worth aligning.
- **`hypothesis_mcp.py` has never run under a real MCP client** — only via `--check`
  and direct calls. Its tool registration is inferred from `keyword-to-immport/server.py`'s idiom.
- **`audit_search_terms` has no overall timeout.** `H2O_TIMEOUT_S` covers gateway calls only;
  ImmPort calls use `keyword-to-immport/server.py`'s httpx defaults.
- **Modules use `print` to stderr, not `logging`.** Safe for MCP (the protocol owns
  stdout) but not what AGENTS.md asks for.
