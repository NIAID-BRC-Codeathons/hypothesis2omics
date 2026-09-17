# MCP servers

This repo exposes several MCP servers. They do not live in one place — one is a
self-contained script, another is a package that shares the repo's `pyproject.toml`
— so instead of a fixed directory, **a server opts in by dropping an
`mcp-server.json` next to its code.**

| Server | Manifest | What it does |
|---|---|---|
| `hypothesis-parser` | `mcp/parse-hypo/` | a plain-text hypothesis → a test spec and an ImmPort search spec (needs an LLM gateway) |
| `keyword-to-immport` | `mcp/keyword-to-immport/` | keywords or a search spec → ImmPort study accessions (public API, no credentials) |
| `hypothesis2omics` | `mcp/data_normalizer/` | run ingestion and inspect GEO units |

Together they chain: parse a hypothesis into a search spec, run that spec against
ImmPort for accessions, then ingest and inspect them.

## Install into Loom

Loom's **Pi.dev extensions** are hardcoded as `-e` args in `bin/loom.js`, so one of
those cannot be installed from outside. MCP servers are a different mechanism and are
fully open: Loom merges `$PI_CODING_AGENT_DIR/mcp.json` (default `~/.pi/agent/mcp.json`)
on every launch, sets only its own `galaxy` and `brc-analytics` keys, and leaves foreign
keys alone. One line registers every server in this repo there:

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
`hypothesis2omics_list_analysis_units`, and so on.

### Before first use

- **`uv` must be on Loom's PATH at launch.** Orbit bundles its own; the CLI uses yours.
- **Pre-warm the Python environment.** `hypothesis2omics` runs out of this repo's
  `pyproject.toml`, which pulls the full scientific stack (~230 MB). Installing that
  on first launch outruns the agent's MCP startup window and the server just looks
  broken. Do it once, up front:

  ```bash
  uv sync --directory ~/.loom/mcp/hypothesis2omics
  ```

- **`IMMPORT_API_KEY`** must be exported in the environment you start Loom from if you
  want `hypothesis2omics_fetch_immport_studies`. Everything else works without it.
- **`OPENAI_API_KEY`** likewise, for `hypothesis_parser_*` — those two steps call an
  LLM. It defaults to Argo, where the "key" is your `ac.yourname` username (`ARGO_USER`
  still works) and you need the Argonne-auth network. To use anything else, export
  `OPENAI_BASE_URL` **and** `OPENAI_MODEL` before starting Loom; the default model id
  is Argo's and means nothing elsewhere. The other servers need none of this.
- **Loom's web/remote shell will not expose these tools.** Its `web-mode-gate` is
  default-deny with an allowlist covering only the prefixes `galaxy_`, `brc_analytics_`,
  `gtn_`, `notebook_`, plus the exact tool `skills_fetch`. Those prefixes are the gate's
  own allowlist, not the `-e` extension list. CLI and Orbit are unaffected.

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

## Instructions, not just tools

Registering the servers gives Loom the tools. It does not tell it how they chain, which
calls cost money, or which stages have no tool at all. Two files in the repo root cover that:

- **`LOOM.md`** — Loom discovers it by walking the cwd's ancestors and injects it
  automatically, so a clone needs no install step. It is injected as *data, not
  instructions*: Loom tells the model the file may have shipped with the folder rather than
  been written by the user, so imperative text in it carries no authority. Written as a map
  and a tool inventory for that reason. Caps are 8 KB and 200 lines, silently truncated past
  either.
- **`LOOM.global.md`** — a short snippet for `~/.pi/agent/LOOM.md`, which *does* ride in the
  cached system prompt with real authority. Opt-in, and **appended** rather than copied:
  that path is global to every project the user opens.

The full route lives in `.claude/skills/hypothesis2omics/SKILL.md`, which both files point
at, together with `scripts/h2o_preflight.py` — a standard-library preflight that reports
environment readiness, which stage is runnable, and whether every artifact on disk carries
provenance.

Loom also has a real skills mechanism — `skills_fetch` plus the Claude Code `SKILL.md`
convention, with progressive disclosure into the system prompt. It is not usable from here:
`shared/loom-config.js` pins `ALLOWED_SKILLS_PREFIX` to `https://github.com/galaxyproject/`
and enforces it twice, because fetched skill text is treated as authoritative instructions
and an arbitrary repo is a prompt-injection vector. Reaching it would mean landing the skill
in a `galaxyproject/*` repository.

## Other MCP clients

These servers are not Loom-specific. For Claude Code, see
[`keyword-to-immport/README.md`](keyword-to-immport/README.md); the same JSON works
for Claude Desktop, Cursor and the rest, only the file it goes in differs.
