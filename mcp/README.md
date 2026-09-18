# MCP servers

This repo exposes several MCP servers. They do not live in one place — one is a
self-contained script, another is a package that shares the repo's `pyproject.toml`
— so instead of a fixed directory, **a server opts in by dropping an
`mcp-server.json` next to its code.**

| Server | Manifest | What it does |
|---|---|---|
| `hypothesis-parser` | `parse-hypo/` | build specs; gate ImmPort studies |
| `keyword-to-immport` | `keyword-to-immport/` | find ImmPort accessions |
| `hypothesis2omics` | `data_normalizer/` | ingest and inspect linked data |

The current tool chain is:

```text
hypothesis_parser_parse_hypothesis
    -> human review of 01_parsed.yaml
    -> hypothesis_parser_build_test_spec
    -> keyword_to_immport_search_spec
    -> hypothesis_parser_check_study_readiness
    -> hypothesis2omics_fetch_immport_studies
    -> hypothesis2omics_parse_immport_studies
    -> hypothesis2omics_fetch_linked_geo
    -> hypothesis2omics_parse_linked_geo_matrices
    -> hypothesis2omics_list_analysis_units
```

The older plan-based GEO path remains available through `plan_geo_retrieval`,
`fetch_planned_geo`, and `parse_geo_matrices`.

## Install into Loom

Loom has no user-installable extensions — the Pi.dev extensions in its README are
hardcoded into `bin/loom.js`. It does pick up MCP servers, though: it merges
`~/.pi/agent/mcp.json` on every launch and leaves keys it doesn't own alone. One
line registers every server in this repo there:

```bash
curl -fsSL https://raw.githubusercontent.com/NIAID-BRC-Codeathons/hypothesis2omics/main/mcp/install.sh | bash
```

It clones this repo to `~/.loom/mcp/hypothesis2omics` and registers the servers;
re-run it to update. Override the location with `$H2O_MCP_DIR` and the branch with
`$H2O_MCP_BRANCH`. Flags pass through:

```bash
... | bash -s -- --dry-run   # print the resulting mcp.json, write nothing
... | bash -s -- --remove    # unregister
```

Already have the repo cloned? Skip the curl:

```bash
node mcp/register.mjs        # also takes --dry-run / --remove / --self-test
```

Restart `loom` (or Orbit) and run `/mcp` to confirm. Tools are prefixed with the
server name, hyphens becoming underscores — `keyword_to_immport_search_spec`,
`hypothesis_parser_check_study_readiness`, `hypothesis2omics_fetch_linked_geo`,
and so on.

### Before first use

- **`uv` must be on Loom's PATH at launch.** Orbit bundles its own; the CLI uses yours.
- **Pre-warm the Python environment.** `hypothesis2omics` runs out of this repo's
  `pyproject.toml`, which pulls the full scientific stack (~230 MB). Installing that
  on first launch outruns the agent's MCP startup window and the server just looks
  broken. Do it once, up front:

  ```bash
  uv sync --directory ~/.loom/mcp/hypothesis2omics
  ```

- **ImmPort credentials** are required for
  `hypothesis_parser_check_study_readiness` and
  `hypothesis2omics_fetch_immport_studies`. The readiness tool accepts credentials
  from `IMMPORT_API_KEY_FILE` or `IMMPORT_API_KEY`; the data-normalizer manifest
  passes `IMMPORT_API_KEY`. Export the applicable variable before starting Loom.
- **`OPENAI_API_KEY`** is required only for
  `hypothesis_parser_parse_hypothesis` and
  `hypothesis_parser_build_test_spec`. They default to Argo, where the "key" is your
  `ac.yourname` username (`ARGO_USER` still works) and you need the Argonne-auth
  network. For another gateway, also export `OPENAI_BASE_URL` and `OPENAI_MODEL`;
  the default model ID is Argo-specific. The deterministic parser tools do not call
  the LLM gateway.
- **Loom's web/remote shell will not expose these tools.** Its `web-mode-gate` is
  default-deny with an allowlist covering only `galaxy_`, `brc_analytics_`, `gtn_`
  and `notebook_`. CLI and Orbit are unaffected.

## Adding a server

Drop an `mcp-server.json` in its directory. The file *is* the `mcp.json` entry —
there is no separate schema to learn — with two placeholders and one extra key:

| | |
|---|---|
| `${dir}` | the directory holding the manifest |
| `${root}` | the repo root |
| `"name"` | the `mcp.json` key, and so the tool prefix. Defaults to the directory name. |

Any other `${VAR}` is left alone for `pi-mcp-adapter` to interpolate from the
environment when it spawns the server — that is how credentials get in without
being written to disk.

A self-contained [PEP 723](https://peps.python.org/pep-0723/) script:

```json
{
  "command": "uv",
  "args": ["run", "--directory", "${dir}", "server.py"],
  "directTools": true
}
```

A package that shares the repo's `pyproject.toml`:

```json
{
  "name": "hypothesis2omics",
  "command": "uv",
  "args": [
    "run",
    "--directory",
    "${root}",
    "python",
    "${root}/mcp/data_normalizer/server.py"
  ],
  "directTools": true,
  "env": { "IMMPORT_API_KEY": "${IMMPORT_API_KEY}" }
}
```

An HTTP server needs only `{"url": "..."}` instead of `command`/`args`.

`node mcp/register.mjs --self-test` checks that every manifest parses, that its
placeholders resolve to paths that exist, and that registering leaves Loom's own
`galaxy` / `brc-analytics` entries untouched.

## Other MCP clients

These servers are not Loom-specific. For Claude Code, see
[`keyword-to-immport/README.md`](keyword-to-immport/README.md); the same JSON works
for Claude Desktop, Cursor and the rest, only the file it goes in differs.
