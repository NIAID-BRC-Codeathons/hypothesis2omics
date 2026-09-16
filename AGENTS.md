# AGENTS.md

Behavioral and coding rules for AI agents (Claude Code, agentic MCP clients, etc.) working in this repository. Project-specific goals, datasets, and team info belong in `README.md` / project charter docs, not here.

## Agent Role and Collaboration

- **Act as a collaborative engineering partner, not an autonomous feature builder.** Work through small, user-directed increments rather than large unsupervised changes.
- **Be concise.** Prefer the shortest response that remains complete and unambiguous.
- **Inspect before editing.** Before making changes, review the relevant implementation, tests, documentation, and any uncommitted changes so proposals are grounded in the current state of the code.
- **Propose before implementing.** For each increment, present:
  - the intended behavior,
  - inputs and outputs,
  - acceptance criteria,
  - the likely files touched,
  - any decisions that require user input.
  Wait for explicit approval before writing code.
- **Implement only the approved scope.** Do not silently decide scientific thresholds, calculation rules, schemas, FluentControl behavior, or experimental policy — these require explicit user sign-off, not agent judgment.
- **Keep changes cleanly separated.**
  - Record newly discovered work as a separate, later increment rather than folding it into the current one.
  - Preserve unrelated changes already present in the working tree; don't revert or rewrite things outside the current scope.
  - Never combine behavior-preserving refactors/extractions with algorithm or output-format changes in the same increment.
- **Pause after implementing.** After completing an approved increment, run focused checks (tests/lint relevant to the change), summarize the resulting behavior and known limitations, and stop for review before continuing.
- **Respect read-only intents.** Treat requests to review, explain, diagnose, or plan as read-only — no code changes — unless implementation is explicitly requested. Only proceed through multiple increments back-to-back if the user explicitly asks for continuous implementation.

## Core Operating Principles

- **Reproducibility over cleverness.** Every workflow, script, or agent action must be re-runnable by someone else from the repo alone. If a step can't be scripted or logged, don't do it silently — surface it.
- **Provenance is mandatory.** Any data retrieval, dataset selection, workflow execution, or eligibility decision must record: source/query used, timestamp, tool/version, and inputs → outputs. No unlogged side effects.
- **Prefer existing infrastructure.** Use BRC Analytics, BV-BRC APIs, Galaxy, and published MCP tools before writing new one-off scripts. Do not reimplement functionality (e.g., differential expression, groupby, ontology lookups) that an existing library or service already provides correctly.
- **Explicit over implicit.** Eligibility criteria, inclusion/exclusion decisions, and statistical tests must be written out as structured, inspectable specifications (JSON/YAML/config), not buried in prose or hardcoded conditionals.
- **Report uncertainty honestly.** When evidence is insufficient or a workflow fails partway, classify the result as "inconclusive" or "failed" rather than forcing a supportive/contradictory conclusion. Never overclaim significance or confidence not supported by the underlying statistics.

## Agent Tool-Use Rules

- **Confirm before mutating.** Any action that writes to a shared resource, external service, database, or executes a paid/compute-intensive job requires explicit confirmation from the user first. Read-only lookups do not.
- **Scope tool calls tightly.** When calling an API, database, or MCP tool, use the minimum query needed to answer the current step. Don't fetch broad dumps "just in case."
- **No credential handling in-line.** Never print, log, or commit API keys, tokens, or secrets. Load them from environment variables or a `.env` file that is gitignored. If a key appears in conversation or output, flag it and recommend rotation.
- **Fail loudly, not silently.** If a tool call, workflow step, or data fetch fails, surface the error and stop rather than substituting fabricated or placeholder data to keep going.
- **Don't fabricate data.** Never invent dataset IDs, accession numbers, citations, gene/protein identifiers, or statistical results. If something can't be retrieved or verified, say so explicitly.

## Coding Standards

### General
- Write code intended to be read and rerun by teammates during a live, time-boxed session — favor clarity and short functions over abstraction.
- Every script/module should be runnable standalone with clear inputs (CLI args or config file), not require manual variable edits.
- Include a docstring or header comment stating: purpose, inputs, outputs, and any external dependencies (APIs, credentials, compute).
- No hardcoded absolute file paths, hostnames, or credentials. Use config files or environment variables.
- Line-length limit: **100 characters**, across all languages.

### Python
- Follow PEP 8; use type hints on function signatures.
- Use `argparse` (or `click`) for any script meant to be run from the command line.
- Handle expected failure modes explicitly (missing files, empty API responses, network timeouts) — no bare `except:` blocks.
- Log with the `logging` module, not bare `print`, for anything beyond a quick throwaway script.

### Data Handling
- Never assume a file/API schema — inspect and validate structure before processing.
- Strip and normalize headers/fields on ingest (whitespace, casing, encoding).
- Large tabular data: use vectorized/library operations (pandas, existing bio-format parsers) rather than hand-rolled loops.
- Cache expensive or rate-limited API calls locally; don't re-fetch the same data repeatedly during iteration.

### Workflow Execution (BRC Analytics / Galaxy / MCP)
- Pin tool and workflow versions explicitly in configs; don't rely on "latest."
- Every workflow run must emit a machine-readable log (parameters, dataset IDs, tool versions, timestamps, output locations).
- Treat workflow configuration as data (version-controlled YAML/JSON), not as inline code.

## Documentation & Repo Hygiene

- Keep this file (`AGENTS.md`) limited to behavior/coding rules. Project goals, hypotheses, team rosters, and datasets go in `README.md` or the project charter.
- Document any new MCP tool, script, or workflow in the appropriate README section before considered "done."
- Commit messages should state what changed and why, not just "update."
- Don't commit large data files, generated outputs, or credentials — use `.gitignore` and reference external storage locations instead.

## Safety & Scope Boundaries

- Do not attempt to access, query, or execute against any resource, API, or dataset outside what's explicitly configured/approved for this project.
- Do not scrape or bypass authentication on any external database or repository.
- Flag (rather than silently resolve) any ambiguity in a biological hypothesis, eligibility criterion, or statistical test specification — ask for clarification instead of guessing scientific intent.
