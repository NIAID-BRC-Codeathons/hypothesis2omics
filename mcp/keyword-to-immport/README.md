# keyword-to-immport

A mini MCP server that turns keywords (or a whole search spec) into ImmPort study
accession numbers. It wraps exactly one endpoint — `GET /api/search/study`
([study_search](https://www.immport.org/data/query/swagger-ui/index.html#/Study%20Search/study_search))
— which is public, so no ImmPort account or API token is needed.

Everything lives in `server.py`. Dependencies are declared inline
([PEP 723](https://peps.python.org/pep-0723/)), so `uv` installs them on first run —
there is no virtualenv to create and no `requirements.txt` to keep in sync.

## Requirements

- [uv](https://docs.astral.sh/uv/getting-started/installation/) (`curl -LsSf https://astral.sh/uv/install.sh | sh`)
- Python 3.11+ (uv will fetch one if you don't have it)
- Network access to `www.immport.org`

## Quick check

Run the built-in self-check, which hits the live API:

```bash
uv run server.py --check
# ok SDY1411 ['SDY1264', 'SDY1289', 'SDY271', ...]
```

## Add it to Claude Code

From this directory:

```bash
claude mcp add immport -- uv run --directory "$PWD" server.py
```

Scopes (`-s`, default `local`):

| Scope | Flag | Stored in | Use it for |
|---|---|---|---|
| local | *(default)* | your user settings, this project only | personal use |
| project | `-s project` | `.mcp.json` in the repo (committed) | sharing with the team |
| user | `-s user` | your user settings, all projects | use it everywhere |

Verify and inspect:

```bash
claude mcp list
claude mcp get immport
claude mcp remove immport
```

Inside a Claude Code session, `/mcp` shows connection status and the tools the
server exposes.

### Or configure it by hand

`.mcp.json` (project scope) or `~/.claude.json` (user scope):

```json
{
  "mcpServers": {
    "immport": {
      "command": "uv",
      "args": ["run", "--directory", "/absolute/path/to/keyword-to-immport", "server.py"]
    }
  }
}
```

`--directory` matters: without it, `uv` resolves `server.py` against whatever
directory Claude Code happens to start the process in.

Same JSON works for any other MCP client (Claude Desktop, Cursor, ...); only the
file it goes in differs.

## Add it to Loom

Loom has no user-installable extensions -- the ones in its README are Pi.dev
modules hardcoded into `bin/loom.js`. It does pick up MCP servers, though: it
merges `~/.pi/agent/mcp.json` on every launch and leaves keys it doesn't own
alone. One line registers every server in this repo there:

```bash
curl -fsSL https://raw.githubusercontent.com/NIAID-BRC-Codeathons/hypothesis2omics/loom-mcp-install/mcp/install.sh | bash
```

It clones this repo to `~/.loom/mcp/hypothesis2omics` and registers the servers;
re-run it to update. Override the location with `$H2O_MCP_DIR` and the branch
with `$H2O_MCP_BRANCH`. Flags pass through:

```bash
... | bash -s -- --dry-run   # print the resulting mcp.json, write nothing
... | bash -s -- --remove    # unregister
```

Restart `loom` (or Orbit) and run `/mcp` to confirm. The tools appear prefixed
with the directory name: `keyword_to_immport_search_spec`,
`keyword_to_immport_search_studies`, `keyword_to_immport_list_facet_values`.

Already have the repo cloned? Skip the curl and run the registration step
directly:

```bash
node mcp/register.mjs           # also takes --dry-run / --remove / --self-test
```

Two caveats: `uv` must be on Loom's PATH at launch (Orbit bundles its own, the
CLI uses yours), and Loom's **web/remote shell won't expose these tools** -- its
`web-mode-gate` is default-deny with an allowlist that covers only Galaxy, BRC
Analytics, GTN, and notebook tools. CLI and Orbit are unaffected.

## Tools

### `search_spec(spec)` — the main one

Takes a YAML or JSON spec and returns **accession IDs only**, best match first.

```yaml
repositories: [ImmPort]          # optional, ignored (ImmPort only)
search_terms:                    # group name -> terms; terms within a group are ORed
  intervention: [YF-17D, YF17D, yellow fever vaccine, yellow fever vaccination]
  outcome: [CD8, CD8+ T cell, T-cell response, vaccine response]
  assay: [transcriptomics, RNA-seq, microarray, gene expression]
require: [intervention, assay]   # optional: groups a study MUST match
min_groups: 2                    # optional: min groups matched (default: len(require) or 1)
filters:                         # optional facet filters, ANDed with every term query
  species: Homo sapiens
limit: 100                       # optional, max accessions returned (default 100, max 1000)
per_term: 200                    # optional, hits fetched per term (default 200, max 1000)
```

```
["SDY1264", "SDY1289", "SDY271", "SDY1294", "SDY1291", "SDY1529"]
```

Ranking: most groups matched first, then summed relevance score.

Group names are yours — `intervention`, `assay`, `whatever` — they only control
grouping and the `require` / `min_groups` logic.

### `search_studies(keywords, limit=10, exact_phrases=False, filters=None)`

Single query, returns accession + title + research focus + score. Use it for
one-off lookups; use `search_spec` when you have more than one concept.

### `list_facet_values(facet=None)`

Valid values for `filters`, most-used first — `species`, `assayMethod`,
`conditionOrDisease`, `researchFocus`, `biosampleType`, `programName`, `sex`,
`race`, `ethnicity`, `clinicalTrial`, `hasLabTest`, `hasAssessment`. These are
controlled vocabulary: `assayMethod: RNA sequencing` works, `RNA-seq` does not.

## Two ImmPort quirks worth knowing

1. **There is no OR operator.** `term` ANDs every word it is given, so
   `influenza malaria` returns 0 hits. That is why `search_spec` issues one
   request per term and unions the results itself — a 14-term spec is 14 requests.
2. **Only study metadata is indexed, not gene lists.** `GCN2` and `EIF2AK4` both
   return 0 hits. Gene-level predictors have to be filtered downstream from the
   expression data, not here.

## Not implemented

Pagination past 1000 hits per term, subject search, TSV output, other
repositories (GEO/SRA), result caching.
