# data_normalizer

MCP server name: **`hypothesis2omics`**

Runs the ingestion stages and exposes the parsed GEO analysis units for
inspection. Registered by `mcp/register.mjs` alongside the other two servers;
see `mcp/README.md` for installation.

```bash
uv run --directory . python mcp/data_normalizer/server.py  # stdio
node mcp/register.mjs --self-test                       # validate manifests
```

The server is `fastmcp` over stdio. It does not listen on a port.

---

## What it is for

The other two servers turn a hypothesis into a search. This one runs the
retrieval and parsing that the search implies, and then lets a client look at
what came back without loading expression matrices into the conversation.

Every pipeline tool returns a bounded status summary with artifact or provenance
paths where applicable. The underlying data modules write the detailed provenance
records. The server's instructions tell the client to run stages one at a time and
read each status before continuing, because later stages must not treat a partial
or failed fetch as complete input.

---

## Tools

### Inspecting parsed output

| Tool | Returns |
|---|---|
| `list_analysis_units()` | units with sample and probe counts |
| `get_analysis_unit(unit_id)` | provenance, hashes, and metadata fields |
| `get_sample_metadata(...)` | filtered, paginated ImmPort and GEO metadata |
| `get_expression_info(unit_id)` | matrix dimensions, samples, path, and hash |

`get_expression_info` deliberately returns no expression values. A U133 Plus 2.0
unit is 20,077 probes by 87 arrays; handing that to a language model is not
inspection, it is a denial of service on the context window. Read the matrix
from the path it gives you.

### Running the pipeline

| Tool | Stage |
|---|---|
| `fetch_immport_studies(...)` | fetch each study's newest parser-required Tab ZIP |
| `parse_immport_studies()` | build combined sample and structured GSE links |
| `inventory_downloaded_files()` | identify cached files supported by parsers |
| `fetch_linked_geo()` | validate links and fetch deduplicated GSEs |
| `parse_linked_geo_matrices()` | derive linked units and parse cached matrices |
| `plan_geo_retrieval(...)` | resolve linked GSMs through the legacy plan path |
| `fetch_planned_geo()` | fetch only the GSEs the plan marks `download` |
| `parse_geo_matrices()` | parse the legacy plan into bounded units |

The current linked-data path is `fetch_immport_studies`, `parse_immport_studies`,
`fetch_linked_geo`, then `parse_linked_geo_matrices`. The ImmPort parser writes the
combined `sample_manifest.tsv` and `geo_series_links.tsv`. Before any GEO request,
the linked fetch validates that every study with GEO-linked samples has a structured
GSE link, deduplicates the GSE accessions, and hashes both selection inputs.

Matrix parsing preserves GEO expression values as submitted. It does not normalize,
annotate, or decide whether a dataset is scientifically eligible.

The plan-based tools remain available for compatibility. That path is
`plan_geo_retrieval`, `fetch_planned_geo`, then `parse_geo_matrices`.

---

## Configuration

| Variable | Default | What it sets |
|---|---|---|
| `IMMPORT_API_KEY` | none | ImmPort credentials, passed through by the manifest |
| `HYPOTHESIS2OMICS_DATA_ROOT` | `data` | where executable stages read and write |
| `HYPOTHESIS2OMICS_GEO_PARSED_ROOT` | `<data root>/geo_cache/parsed` | inspection root |

`create_server(parsed_root, data_root)` takes both as arguments too, which is
how to point a test at a fixture directory instead of the real cache.
Executable MCP callers cannot supply filesystem paths; all pipeline paths are
derived from the configured data root.

**Never commit an API key.** `.gitignore` blocks `**/immport-key-*.json`. Each
person generates their own; see the credentials section in the root `README.md`.

---

## Notes

The server adds `mcp/` and the repository root to `sys.path` at import time so
`data_normalizer.*` and the `data/` modules resolve whichever directory it is
launched from. That is why the manifest can use `${root}` rather than requiring
a particular working directory.

The adapters call the existing functions under `data/` for fetching and parsing.
They add bounded path selection, result summaries, and—for linked GEO fetching—an
MCP-specific validation and selection-provenance layer. The underlying command-line
programs remain available for direct, explicitly configured runs.
