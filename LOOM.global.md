## Working in a repository that ships its own tools

When a repository registers MCP servers or ships analysis scripts, those are the
implementation. Call them. Do not write code that reproduces what one of them does — a
replacement usually runs and returns a plausible number, which is why this is worth stating
rather than leaving to judgment.

This applies specifically to: parsing repository data formats, querying a data repository's
API, downloading from a public archive, differential expression, and any eligibility or
statistical verdict the repository already has a script for.

When a tool fails, report the failure and stop. Do not substitute a hand-rolled fallback to
keep going. A fallback result carries no provenance, and nothing downstream can distinguish
it from a real one.

Do not choose a scientific threshold on my behalf — significance level, expected direction,
effect-size floor, timepoint, or inclusion criterion. Ask me. A defensible default chosen
silently is the failure these pipelines are built to prevent.
