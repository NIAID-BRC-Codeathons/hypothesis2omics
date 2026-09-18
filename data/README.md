# Data pipeline and analysis inputs

Run both fetchers from the repository root after installing dependencies with `uv sync`.

Fetchers and parsers write under `data/`. Most stages emit provenance, but
`extract_series_matrix.py` currently writes intensity CSVs without a provenance file. This
checkout contains selected committed benchmark cache artifacts; treat new raw downloads as
generated data and do not commit them.

## Input selection

| Source | Behavior when accessions are omitted |
|---|---|
| ImmPort | `SDY1529`, `SDY1264`, `SDY1294`, `SDY1289` |
| GEO | read GSEs from `immport_cache/parsed/geo_series_links.tsv` |

Explicit positional accessions replace the corresponding automatic selection.

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

Discover and parse the newest extracted Tab directory or Tab ZIP for every study under the
cache root:

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

Each study/experiment/GSE/GPL unit receives a compressed probe-by-sample expression matrix,
linked ImmPort sample rows, lossless long-form GEO sample metadata, and provenance. Values are
preserved as submitted; this step performs no normalization, annotation, or eligibility filtering.

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

The config supports a `datasets` array. Each entry identifies its study, GSE, platform, gene,
feature ID, outcome names, timepoint unit, source matrix, and ImmPort Tab ZIP. Exact overlapping
rows are deduplicated when several entries share a study or outcome.

**The gene and feature ID are entered by hand.** Nothing in this repository maps a probe to a
gene, so `gene: EIF2AK4` beside `feature_id: Hs.412102_at` is an assertion, not an annotation
lookup. The provenance records both values as supplied selections. Confirm these pairings against
the GPL platform annotations before interpreting the gene-labelled output.

## Export neutralizing-antibody results for Galaxy

Export raw `neut_ab_titer_result.txt` rows from each applicable study's newest Tab source:

```bash
uv run python data/immport_neut_ab_export.py
```

This writes `ImmPort_neut_ab_titer_results.tsv` and its provenance JSON under
`data/galaxy_file_input/`. Source columns are retained with lowercase headers. Preferred and
reported values remain separate, so censored values such as `<10` are not substituted for blank
`value_preferred` fields. Studies without a neutralizing-antibody result table are skipped and
listed in the provenance file; malformed or duplicate tables still stop the export.

## Extract Galaxy-ready intensity matrices

Convert GEO series-matrix files into the intensity CSVs the analysis stage consumes:

```bash
# every GSE* folder under a directory, recursively
uv run python data/extract_series_matrix.py --input-dir data/geo_cache
```

In directory mode the extractor finds folders whose names begin with `GSE`, searches each one
recursively for series-matrix files, reads `.txt` and `.txt.gz` directly with no manual
decompression, extracts the table between `!series_matrix_table_begin` and
`!series_matrix_table_end`, renames GEO's `ID_REF` column to `probe_id`, and writes the result to
`data/galaxy_file_input/`.

Tested on `GSE13485` (20,077 probes x 87 samples) and `GSE13699` (22,184 probes x 126 samples).

The matrix is one of the two files limma needs. The design file naming each sample's group is still
prepared per dataset by hand. This extractor does not currently write provenance.

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

Both calls return one structured record per accession with status, artifact paths, sizes, hashes,
errors, and timing. ImmPort records additionally include selection mode, selected release, and file
counts. GEO records include source URLs plus platform and sample counts. Pass `all_files=True` to
the ImmPort fetcher only when every remote file is needed.

## Outputs and options

Retrieval and parsing artifacts:

```text
data/
|-- immport_cache/
|   |-- <SDY_ID>/
|   |-- manifest.json
|   |-- parsed/                  # combined sample and structured GSE links
|   `-- provenance_log.jsonl
`-- geo_cache/
    |-- <GSE_ID>/
    |-- manifest.json
    |-- parsed/                  # bounded analysis units and run manifest
    |-- geo_fetch_selection.provenance.json
    `-- provenance_log.jsonl
```

Other generated artifacts:

```text
data/
|-- file_inventory/        # downloaded_files.tsv, parser_candidates.tsv, provenance
|-- galaxy_file_input/     # intensity CSVs; neutralizing-antibody TSV and provenance
|-- validator_input/       # the combined validator handoff bundle
`-- tests/                 # unittest suites for the modules in this directory
```

The two fetch commands write `manifest.json` and append to `provenance_log.jsonl`. Existing valid
files are reused. Use `--force` to download them again or `--destdir` to select another cache root.
Run a program with `--help` for its complete options.
