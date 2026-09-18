# Hypothesis2Omics: Agentic Multi-Omics Validation

**NIAID-BRCs AI Codeathon 2.0** · September 16–18, 2026 · Argonne National Laboratory

Translating a biological hypothesis into a reproducible search and analysis of public transcriptomic, proteomic, metabolomic, microbiome, or immune-system datasets.

Project page: https://niaid-brc-codeathons.github.io/projects/hypothesis2omics/

---

## What this does

A plain-English hypothesis goes in. A search over public repositories, an eligibility
judgment on every candidate dataset, a real analysis on the eligible ones, and an evidence
report come out. Every stage records its provenance, and the report says what would have
counted as an answer before any numbers existed.

The design goal is not to find a new result. It is to be honest about what the data can and
cannot show. A pipeline that says `supports` when handed a hypothesis its datasets cannot
test is worse than useless, so the components below are built to refuse rather than to
conclude when refusing is correct.

## Pipeline

```text
Hypothesis
    -> Test specification              mcp/parse-hypo
    -> ImmPort dataset discovery       mcp/keyword-to-immport
    -> ImmPort and GEO retrieval       data/, mcp/data_normalizer
    -> Normalized validator input      data/validator_input
    -> Scientific eligibility          scientific_validator
    -> Analysis                        Galaxy / limma, or a per-subject regression
    -> Evidence synthesis and report   evidence_rules
```

A human reads and corrects the parsed hypothesis before anything proceeds. Nothing in this
pipeline runs unattended from end to end, by design.

| Stage                                           | Where                          | Owner                      |
| ----------------------------------------------- | ------------------------------ | -------------------------- |
| Hypothesis -> test spec -> search spec          | `mcp/parse-hypo`               | Melanie Sadecki            |
| Keyword and spec search over ImmPort            | `mcp/keyword-to-immport`       | Yaphet Kebede              |
| MCP registration and installer                  | `mcp/`                         | Yaphet Kebede              |
| Retrieval, parsing, normalization               | `data/`, `mcp/data_normalizer` | Yijun Zhou, Amar Kumar     |
| Scientific eligibility and evidence fusion      | `scientific_validator/`        | Amar Kumar                 |
| Analysis execution                              | Galaxy, limma                  | Archit Vasan, Slim Fourati |
| Decision rules, independence, synthesis, report | `evidence_rules/`              | Rushikesh Lagad            |

### Current limitation, stated plainly

The step between "this dataset is eligible" and "the analysis has run" is only partly
automated. `data/extract_series_matrix.py` scans one or more `GSE*` dataset folders, reads
GEO series-matrix files directly from `.txt` or `.txt.gz` format, and generates Galaxy-ready
intensity CSV files in `data/galaxy_file_input/`. The design file naming each sample's group
is still prepared per dataset by hand, and limma needs both. The extractor does not currently
write provenance; retrieval and validator handoff do.

## Folder layout

```text
hypothesis2omics/
|-- mcp/
|   |-- parse-hypo/               # Hypothesis -> test specification
|   |-- keyword-to-immport/       # Search terms -> ImmPort studies
|   `-- data_normalizer/          # Fetch, parse, and inspect datasets
|
|-- data/
|   |-- immport_cache/            # Downloaded and normalized ImmPort data
|   |-- geo_cache/                # Downloaded and parsed GEO data
|   |-- validator_input/          # Validator-ready tables and provenance
|   |-- galaxy_file_input/        # Galaxy-ready intensity matrices
|   `-- extract_series_matrix.py  # GEO series-matrix -> intensity CSV
|
|-- scientific_validator/         # Dataset eligibility assessment
|
|-- evidence_rules/               # Decision rules, independence, synthesis, report
|
`-- docs/                         # Ground truth and benchmark documentation
```

Detailed data commands and MCP setup are in `data/README.md` and `mcp/README.md`.
The evidence layer is documented in `evidence_rules/README.md`.

## The benchmark

The MVP evaluates the hypothesis that GCN2/EIF2AK4 activity is associated with the magnitude
of the CD8+ T-cell response following YF-17D vaccination. ImmPort supplies study design,
sample linkage, and immune-response measurements; linked GEO series supply transcriptomic
measurements.

| ImmPort study | GEO series               | Platform             | Linked samples |
| ------------- | ------------------------ | -------------------- | -------------: |
| `SDY1264`     | `GSE13485`               | `GPL7567`            |             87 |
| `SDY1289`     | `GSE13699`               | `GPL6104`, `GPL6883` |       126 + 16 |
| `SDY1291`     | `GSE22768`               | `GPL10647`           |             50 |
| `SDY1294`     | `GSE82152`               | `GPL21975`           |            109 |
| `SDY1529`     | `GSE125921`, `GSE136163` | `GPL10558`           |       36 + 144 |

The committed parser run creates seven study/experiment/GSE/platform units containing 568 linked
samples. The `SDY1264` validator handoff extracts 87 expression values for the EIF2AK4 feature
`Hs.412102_at` and 25 quantitative `Act CD8 T Cell Response` measurements, and returns
`scientific_status: eligible`, `execution_status: ready`.

This checkout retains selected benchmark cache, parsed, validator, and Galaxy-input artifacts.
Treat new raw downloads as generated data and do not commit them. The regeneration commands and
known provenance gaps are documented in `data/README.md`.

Expected direction, predictor timepoint and cohort structure are documented in
`docs/ground_truth.md`, taken from Querec et al. 2009 and Ravindran et al. 2014 and checked
against the GEO records.

## What the benchmark found

```bash
python evidence_rules/run_yf17d.py --repo . --report docs/evidence_report_yf17d.md
```

That prints the run and writes the evidence report, plus `docs/evidence_report_yf17d.md.json` for
anything downstream.

**The headline verdict is `inconclusive`.** The pre-registered predictor is EIF2AK4 at day 7
as measured, because that is what Querec's signature was built from, and neither trial in
`GSE13485` reaches alpha on it.

Two other defensible definitions of the same predictor do better. Day 7 minus day 0 is
positive and significant in both trials, which reads as a replication. But baseline EIF2AK4
is itself correlated with the outcome, so a difference score inherits that; holding baseline
constant, the second trial is null. Three specifications, three overall verdicts, the same 25
subjects.

That is the pipeline working. The report leads with the specification that was declared
first, prints the other two next to it with their verdicts, and says in those words that the
result is specification-sensitive. Picking the flattering specification would have been one
line of code and no warning to the reader, which is the failure mode this stage exists to
catch.

## Design commitments

**The decision rule is frozen before results exist.** Alpha, expected direction, minimum
effect and the replication threshold are set during planning and printed in the report, so a
reader can check that what was promised is what was applied.

**The specification that counts is declared, not selected.** "Early expression" has several
defensible definitions and they do not agree. The report generator requires exactly one
specification marked primary, with a stated reason, and refuses to render without it. Every
other specification is printed beside it with its own verdict. This is the difference between
a result and a choice of result, and the report shows which one it is.

**Independence is counted in groups, not rows.** `GSE125921` and `GSE136163` are both
`SDY1529`. `GSE13485` is one accession holding two trials run a year apart. Treating related
accessions as independent replication is the easiest way to overstate a result.

**A result records what it is, not only what it equals.** A Galaxy limma run returns a log
fold change; a per-subject regression returns a correlation. `AnalysisResult` carries the
effect type, the estimand and the contrast direction so the two can sit in one evidence table
without being averaged.

**Controls are not evidence.** A sex-marker contrast with a known answer proves the machinery
works. It carries `role="positive_control"` and is refused entry to any evidence table about
the hypothesis.

**Inconclusive is a real verdict.** A hypothesis whose predictor was never measured, or whose
timepoint does not exist in the data, returns `inconclusive` rather than `refutes`. The
question went unasked, and saying otherwise would be a false claim about the biology.

## Credentials

Never commit or share an API key. Each person generates their own.

* **ImmPort**: create a key at the ImmPort API Keys page and pass the downloaded JSON with
  `--api-key-file`, or set `IMMPORT_API_KEY` in the environment you launch from for the MCP
  server. `.gitignore` blocks `**/immport-key-*.json`.
* **LLM gateway**: `OPENAI_API_KEY` for the hypothesis parser, defaulting to Argo. See
  `mcp/README.md`.

## Evaluation

Dataset-retrieval recall, workflow success, consistency with published findings, robustness
across datasets, and expert assessment of the final evidence report.

## Leads

* Slim Fourati
* Rushikesh Lagad

## Working here

This repository is the team's working space for the codeathon. Team members get access
through the [NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization.
