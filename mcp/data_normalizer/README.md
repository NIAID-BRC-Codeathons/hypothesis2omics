# data_normalizer

MCP server name: **`hypothesis2omics`**

Runs the ingestion stages and exposes the parsed GEO analysis units for
inspection. Registered by `mcp/register.mjs` alongside the other two servers;
see `mcp/README.md` for installation.

```bash
python mcp/data_normalizer/server.py       # stdio, what the manifest launches
node mcp/register.mjs --self-test          # checks this manifest parses and resolves
```

The server is `fastmcp` over stdio. It does not listen on a port.

---

## What it is for

The other two servers turn a hypothesis into a search. This one runs the
retrieval and parsing that the search implies, and then lets a client look at
what came back without loading expression matrices into the conversation.

Every executable tool returns the stage's provenance record rather than a bare
success flag. The server's own instructions tell the client to run stages one
at a time and read each status before continuing, because a failed fetch that
looks like a success produces an analysis on missing data.

---

## Tools

### Inspecting parsed output

| Tool | Returns |
|---|---|
| `list_analysis_units()` | every parsed GEO analysis unit with its sample and probe counts |
| `get_analysis_unit(unit_id)` | one unit's provenance, output hashes, and which metadata fields it has |
| `get_sample_metadata(unit_id, gsm_accessions, attributes, limit, offset)` | a filtered, paginated page of linked ImmPort and GEO sample metadata |
| `get_expression_info(unit_id)` | matrix dimensions, sample list, path and hash, without any values |

`get_expression_info` deliberately returns no expression values. A U133 Plus 2.0
unit is 20,077 probes by 87 arrays; handing that to a language model is not
inspection, it is a denial of service on the context window. Read the matrix
from the path it gives you.

### Running the pipeline

| Tool | Stage |
|---|---|
| `fetch_immport_studies(study_accessions, max_files_per_study)` | fetch named ImmPort studies, server-side credentials, cached |
| `parse_immport_studies()` | parse every cached study and combine their sample linkage |
| `inventory_downloaded_files()` | list what is cached and which files a parser can use |
| `plan_geo_retrieval(batch_size, timeout)` | resolve linked GSMs to series and platforms, write the plan |
| `fetch_planned_geo()` | fetch only the GSEs the plan marks `download` |
| `parse_geo_matrices()` | parse downloaded matrices into bounded analysis units |

Run them in that order. `plan_geo_retrieval` is metadata-only and marks explicit
SuperSeries records `skip_superseries`, so `fetch_planned_geo` does not download
the same expression data twice.

---

## Configuration

| Variable | Default | What it sets |
|---|---|---|
| `IMMPORT_API_KEY` | none | ImmPort credentials, passed through by the manifest |
| `HYPOTHESIS2OMICS_DATA_ROOT` | `data` | where executable stages read and write |
| `HYPOTHESIS2OMICS_GEO_PARSED_ROOT` | `data/geo_cache/parsed` | where the inspection tools look |

`create_server(parsed_root, data_root)` takes both as arguments too, which is
how to point a test at a fixture directory instead of the real cache.

**Never commit an API key.** `.gitignore` blocks `**/immport-key-*.json`. Each
person generates their own; see the credentials section in the root `README.md`.

---

## Notes

The server adds `mcp/` and the repository root to `sys.path` at import time so
`data_normalizer.*` and the `data/` modules resolve whichever directory it is
launched from. That is why the manifest can use `${root}` rather than requiring
a particular working directory.

The command-line equivalents of every executable tool are documented in
`data/README.md`, and they write the same provenance files. The MCP tools are a
second entry point to the same code, not a reimplementation.
