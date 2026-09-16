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
# Fetch selected studies.
uv run python data/immport_fetch_module.py SDY1529 SDY1264 \
  --api-key-file data/immport_cache/immport-key-REPLACE_ME.json

# Limit a smoke test to one file per study.
uv run python data/immport_fetch_module.py SDY1529 \
  --api-key-file data/immport_cache/immport-key-REPLACE_ME.json \
  --max-files-per-study 1
```

You may instead set `IMMPORT_API_KEY_FILE` to the JSON path or `IMMPORT_API_KEY` to the raw key.
Avoid placing a raw key in a command because shell history may retain it.

## Parse ImmPort sample links

Point the parser at an extracted ImmPort `Tab/` directory:

```bash
uv run python data/immport_parse_module.py \
  data/immport_cache/SDY1529/SDY1529-DR58_Tab/Tab \
  --output-dir data/immport_cache/SDY1529/parsed
```

The parser writes `sample_manifest.tsv` (one row per experimental sample) and
`sample_manifest.provenance.json` (source and output hashes, versions, counts, and timing). It
keeps samples without public-repository links and does not select assays or timepoints.

The same command accepts an unextracted `*_Tab.zip` file. To discover and parse the newest
extracted Tab directory or Tab ZIP for every study under a cache root, run:

```bash
uv run python data/immport_batch_parse.py --cache-root data/immport_cache
```

When both forms exist for the same release, the batch parser uses the extracted directory. It can
also read required tables directly from a ZIP without extracting it and does not use MySQL ZIPs.
Per-study outputs remain under `<SDY_ID>/parsed/`. A combined
`sample_manifest.tsv` and `batch_manifest.json` are written under `immport_cache/parsed/`.

## Plan GEO retrieval from ImmPort links

Resolve the combined GSM list into parent GEO series and platforms before downloading data:

```bash
uv run python data/geo_plan_module.py \
  data/immport_cache/parsed/sample_manifest.tsv \
  --output-dir data/geo_cache/plan
```

This metadata-only step writes `geo_download_plan.tsv`, a reusable GSM resolution cache, and a
provenance JSON file. Requests are grouped for NCBI E-utilities, and unresolved accessions remain
visible in the plan. GSE metadata is cached separately; explicit SuperSeries records are marked
`skip_superseries` to prevent redundant downloads. Use `--force` only when cached metadata needs
refreshing.

## GEO acquisition

GEO requires internet access but no account or credentials.

```bash
# Fetch the default candidates.
uv run python data/geo_fetch_module.py

# Fetch selected studies only.
uv run python data/geo_fetch_module.py GSE13699 GSE125921
```

## Parse planned GEO matrices

Parse cached matrices for plan rows marked `download`:

```bash
uv run python data/geo_matrix_parse_module.py \
  data/geo_cache/plan/geo_download_plan.tsv \
  --immport-manifest data/immport_cache/parsed/sample_manifest.tsv \
  --geo-cache-root data/geo_cache \
  --output-dir data/geo_cache/parsed
```

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
    max_files_per_study=1,
    provenance_log_path="data/immport_cache/provenance_log.jsonl",
)
geo_records = fetch_geo_datasets(
    ["GSE13699"],
    provenance_log_path="data/geo_cache/provenance_log.jsonl",
)
```

Each call returns one structured record per accession, including status, source and local paths,
file sizes, SHA-256 checksums, errors, and timing.

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
    └── provenance_log.jsonl
```

Command-line runs write `manifest.json` and append to `provenance_log.jsonl`. Existing valid files
are reused. Use `--force` to download them again or `--destdir` to select another cache root. Run
either script with `--help` for all options.
