# Hypothesis2Omics: Agentic Multi-Omics Validation

**NIAID-BRCs AI Codeathon 2.0** · September 16–18, 2026 · Argonne National Laboratory

Translating a biological hypothesis into a reproducible search and analysis of public transcriptomic, proteomic, metabolomic, microbiome, or immune-system datasets.

Project page: https://niaid-brc-codeathons.github.io/projects/hypothesis2omics/

---

> **This is a draft pitch, not a plan.**
>
> What follows is a one-slide proposal from the organizing team. It exists
> to seed a team, not to constrain one. Scope, methods, target organism,
> and success criteria are all still open — expect them to change
> substantially. Turning this into a real plan is the team's first job, and
> it lands in the project charter due August 28, 2026.

---

## Goal (proposed)

Translate a biological hypothesis into a reproducible search and analysis of public transcriptomic, proteomic, metabolomic, microbiome, or immune-system datasets.

## Three-Day MVP (proposed)

Use one predefined hypothesis, for example, “A specified pathway or biomarker is associated with protection following vaccination.” The agent should discover appropriate datasets, assess eligibility, configure workflows, execute analyses through BRC Analytics or Galaxy, integrate results across datasets, and produce an evidence table.

A microbiome causal-inference branch could expose HUMAnN, differential abundance, sensitivity analysis, protein-language-model annotation, and SHAP interpretation as MCP tools.

## Example datasets

The MVP evaluates the hypothesis that GCN2/EIF2AK4 activity is associated with the magnitude of
the CD8+ T-cell response following YF-17D vaccination. ImmPort supplies study design, sample
linkage, and immune-response measurements; linked GEO series supply transcriptomic measurements.

| ImmPort study | GEO series | Platform | Linked samples |
|---|---|---|---:|
| `SDY1264` | `GSE13485` | `GPL7567` | 87 |
| `SDY1289` | `GSE13699` | `GPL6104`, `GPL6883` | 126 + 16 |
| `SDY1294` | `GSE82152` | `GPL21975` | 109 |
| `SDY1529` | `GSE125921`, `GSE136163` | `GPL10558` | 36 + 144 |

The current retrieval plan also resolves `GSE13486`, the SuperSeries containing `GSE13485`, but
marks it `skip_superseries` to avoid downloading duplicate expression data. `SDY1291` is fetched
and parsed successfully but does not currently produce a GEO analysis unit.

The reproducible parser run creates six study/experiment/GSE/platform analysis units containing
518 linked samples. The configured `SDY1264` validator handoff extracts 87 expression values for
the EIF2AK4 feature `Hs.412102_at` and 25 quantitative `Act CD8 T Cell Response` measurements.
Downloaded datasets and generated outputs are excluded from Git; the acquisition and parsing
steps below recreate them with provenance.


## Architecture (proposed)

<img width="1536" height="1024" alt="ChatGPT Image Sep 16, 2026, 11_36_29 AM" src="https://github.com/user-attachments/assets/87c40a8e-767e-42f6-a9ca-02c393795b10" />

## Pipeline

General pipeline flow
```text
Hypothesis
    → ImmPort dataset discovery
    → ImmPort and GEO normalization
    → Validator input
    → Scientific eligibility assessment
```

Folder layout
```text
hypothesis2omics/
├── mcp/
│   ├── parse-hypo/          # Hypothesis → test specification
│   ├── keyword-to-immport/  # Search terms → ImmPort studies
│   └── data_normalizer/     # Fetch, parse, and inspect datasets
│
├── data/
│   ├── immport_cache/       # Downloaded and normalized ImmPort data
│   ├── geo_cache/           # Downloaded and parsed GEO data
│   └── validator_input/     # Validator-ready tables and provenance
│
└── scientific_validator/    # Dataset eligibility assessment
```


Detailed data commands and MCP setup are documented in `data/README.md` and `mcp/README.md`.


## Evaluation (proposed)

Dataset-retrieval recall, workflow success, consistency with published findings, robustness across datasets, and expert assessment of the final evidence report.

## Leads

- Slim Fourati
- Rushikesh Lagad

Team assignments are still being finalized. Participants can review their project, and request a reassignment, in the participant spreadsheet circulated by the organizing team.

## Working here

This repository is the team's working space for the codeathon — code, notebooks, data pointers, and notes. Replace this README with the real thing once the charter is written. Team members get access through the [NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization; accept the invitation if you have not already.
