# Dataset acquisition

Run both fetchers from the repository root after installing dependencies with `uv sync`. Generated
data and provenance are cached locally in directories excluded from Git.

## Default candidates

| Source | Accessions used when none are supplied |
|---|---|
| ImmPort | `SDY1529`, `SDY1264`, `SDY1294`, `SDY1289` |
| GEO | `GSE125921`, `GSE136163`, `GSE13485`, `GSE82152`, `GSE13699` |

Passing accessions on the command line replaces the corresponding default list.

## ImmPort acquisition

ImmPort requires a personal API key with `browse` and `download` scopes. Create one on the
[ImmPort API Keys page](https://www.immport.org/auth/api/keys) and save the downloaded JSON file
under `data/immport_cache/`. Never commit or share it.

```bash
# Fetch each study's newest Tab ZIP required by the ImmPort parser.
uv run python data/immport_fetch_module.py SDY1529 SDY1264 \
  --api-key-file data/immport_cache/immport-key-REPLACE_ME.json

# Download every file listed by ImmPort only when the full dataset is needed.
uv run python data/immport_fetch_module.py SDY1529 \
  --api-key-file data/immport_cache/immport-key-REPLACE_ME.json \
  --all-files
```

Always pass the downloaded JSON file with `--api-key-file`, as shown above. Replace
`immport-key-REPLACE_ME.json` with the downloaded filename; do not paste the raw key into a command.
The complete remote file manifest is retained in both modes. By default, the fetcher selects the
highest-release `<SDY_ID>-DR<n>_Tab.zip`; MySQL archives and result files are not downloaded.
`--max-files-per-study` can optionally cap the files selected for processing, which is most useful
with `--all-files`.

## Parse ImmPort sample links

Point the parser at an extracted ImmPort `Tab/` directory:

To discover and parse the newest
extracted Tab directory or Tab ZIP for every study under a cache root, run:

```bash
uv run python data/immport_batch_parse.py --cache-root data/immport_cache
```
A combined `sample_manifest.tsv`, structured `geo_series_links.tsv`, and `batch_manifest.json` are
written under `immport_cache/parsed/`. The series-link table extracts validated `GSE<number>`
accessions from each study's `study_link.txt` for downstream GEO retrieval.

## GEO acquisition

```bash
# Fetch the GSE accessions listed in the parsed ImmPort GEO links.
uv run python data/geo_fetch_module.py

# Override the parsed links and fetch selected studies only.
uv run python data/geo_fetch_module.py GSE13699 GSE125921
```

With no positional accessions, the fetcher reads
`data/immport_cache/parsed/geo_series_links.tsv`, validates and deduplicates its
`gse_accession` values, and fails before downloading if that input is missing or malformed. Use
`--links-file` to select a different structured link table. Each run writes
`geo_fetch_selection.provenance.json` alongside the existing GEO manifest and provenance log.

## Parse linked GEO datasets

After downloading GEO, parse the cached matrices using the combined ImmPort manifest and
structured GSE links:

```bash
uv run python data/geo_matrix_parse_module.py
```

The parser derives study/experiment/GSE/GPL units by matching each ImmPort-linked GSM to exactly
one cached series matrix and its submitted platform metadata. Missing and ambiguous matches fail
before parsing. The derived selection is saved as `geo_matrix_parse_selection.tsv` under
`data/geo_cache/parsed/`. A legacy `geo_download_plan.tsv` can still be supplied as the optional
positional argument.

## Build the validator input bundle

After the configured ImmPort and GEO inputs are available, create the scientific-validator input
bundle:

```bash
uv run python data/validator_handoff_parse_module.py \
  data/validator_handoff_config.json
```

`validator_handoff_config.json` explicitly selects the normalized ImmPort manifest, GEO series
matrix, study, platform, molecular feature, ImmPort Tab ZIP, quantitative outcome names, and output
directory. Review those selections before running the command; the program validates and applies
them but does not infer scientific choices.

The configured sources are converted into these primary files under `data/validator_input/`:

- `sample_manifest.tsv`: a verbatim copy of the normalized ImmPort sample manifest.
- `feature_expression.tsv`: values for the configured GEO molecular feature.
- `quantitative_outcome.tsv`: configured ImmPort outcomes joined to subjects and timepoints.

Each TSV receives a provenance JSON file. `validator_input_manifest.json` records the
configuration, output paths, row counts, hashes, parser version, and timestamps for the complete
bundle.

The MCP data pipeline can automate configuration assembly with its `build_validator_input` tool.
Callers must provide the study, GSE, GPL, gene, feature ID, outcome names, and timepoint unit. The
tool resolves the combined ImmPort manifest, matching cached GEO matrix, and newest cached ImmPort
Tab ZIP; writes `data/validator_input/validator_handoff_config.resolved.json`; and runs the same
handoff parser. Scientific selections remain explicit rather than being inferred by the pipeline.


## Export neutralizing-antibody results for Galaxy

Export the raw `neut_ab_titer_result.txt` rows from each discovered study's newest Tab source:

```bash
uv run python data/immport_neut_ab_export.py
```

This writes `ImmPort_neut_ab_titer_results.tsv` and its provenance JSON under
`data/galaxy_file_input/`. Source columns are retained with lowercase headers. Preferred and
reported values remain separate, so censored values such as `<10` are not substituted for blank
`value_preferred` fields.



Each study/experiment/GSE/GPL unit receives a compressed probe-by-sample expression matrix,
linked ImmPort sample rows, lossless long-form GEO sample metadata, and provenance. Values are
preserved as submitted; this step performs no normalization, annotation, or eligibility filtering.

## Python API

```python
from data.immport_fetch_module import fetch_immport_datasets
from data.geo_fetch_module import fetch_geo_datasets

immport_records = fetch_immport_datasets(
    ["SDY1529"],
    api_key_file="data/immport_cache/immport-key-REPLACE_ME.json",
    provenance_log_path="data/immport_cache/provenance_log.jsonl",
)
geo_records = fetch_geo_datasets(
    ["GSE13699"],
    provenance_log_path="data/geo_cache/provenance_log.jsonl",
)
```

Each call returns one structured record per accession, including selection mode and release, status,
source and local paths, file counts, sizes, SHA-256 checksums, errors, and timing. Pass
`all_files=True` to the Python fetcher only when every remote file is needed.

## Outputs and options

```text
data/
├── immport_cache/
│   ├── <SDY_ID>/
│   ├── manifest.json
│   └── provenance_log.jsonl
└── geo_cache/
    ├── <GSE_ID>/
    ├── manifest.json
    ├── geo_fetch_selection.provenance.json
    └── provenance_log.jsonl
```

Command-line runs write `manifest.json` and append to `provenance_log.jsonl`. Existing valid files
are reused. Use `--force` to download them again or `--destdir` to select another cache root. Run
either script with `--help` for all options.
