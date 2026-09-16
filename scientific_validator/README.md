# Hypothesis2Omics Scientific Validator

A proof-of-concept scientific eligibility layer for **Hypothesis2Omics**.

The goal is to determine whether a candidate omics dataset can actually test a biological 
hypothesis before sending that dataset into downstream analysis.

## Why this is needed

Finding a dataset in ImmPort, GEO, or another repository does not automatically mean that the 
dataset is scientifically suitable.

A candidate study may be missing:

* the required predictor
* the required biological outcome
* appropriate timepoints
* participant-level linkage
* accessible outcome data

This validator separates **scientific eligibility** from **technical executability**.

## Where it fits

```text
Natural-language hypothesis
        ↓
Test specification
        ↓
Dataset discovery / retrieval
        ↓
Dataset parsing
        ↓
Scientific eligibility validator
        ↓
Workflow execution
        ↓
Evidence synthesis / report
```

## Decision model

The validator returns one of three scientific states:

```text
eligible
uncertain
excluded
```

It also reports a separate execution status so that a scientifically relevant dataset is not 
automatically treated as analysis-ready.

## Example

The included proof of concept evaluates the hypothesis:

> GCN2/EIF2AK4 activity is associated with the magnitude of the CD8+ T-cell response following 
YF-17D vaccination.

### SDY1264

```text
scientific_status = eligible
execution_status  = needs_data_check
```

The study satisfies the scientific requirements, but participant-level CD8 outcome data still 
require an accessibility check.

### SDY1529

```text
scientific_status = uncertain
execution_status  = needs_scientific_resolution
```

Several required criteria remain unresolved, including the exact CD8 endpoint, timepoints, and 
participant-level linkage.

## Run

```bash
python3 eligibility_engine.py \
  --spec examples/yf17d_test_spec.json \
  --dataset examples/SDY1264.json
```

Second example:

```bash
python3 eligibility_engine.py \
  --spec examples/yf17d_test_spec.json \
  --dataset examples/SDY1529.json
```

## Run tests

```bash
python3 -m unittest discover -s tests -v
```

Current test suite:

```text
Ran 4 tests
OK
```

## Repository structure

```text
.
├── README.md
├── eligibility_engine.py
├── examples/
│   ├── yf17d_test_spec.json
│   ├── SDY1264.json
│   ├── SDY1529.json
│   └── synthetic_excluded.json
├── output/
│   ├── SDY1264_assessment.json
│   └── SDY1529_assessment.json
└── tests/
    └── test_eligibility_engine.py
```

## Design

The eligibility logic is config-driven.

The engine does not hard-code YF-17D, EIF2AK4, or CD8-specific logic. Instead, the required 
criteria are defined in the test specification.

Each dataset evidence record contains:

* criterion value
* confidence
* supporting evidence
* provenance/source

The final decision is deterministic and auditable.

## Current scope

This repository is a proof of concept.

The current dataset evidence records are manually curated examples used to test the eligibility 
logic independently from repository parsing.

The next integration step is to automatically convert metadata and files returned from ImmPort/GEO 
retrieval tools into the evidence format expected by the validator.

## Planned extensions

* automated ImmPort/GEO metadata mapping
* ontology/synonym-aware criterion matching
* automated analysis specification
* MCP tool exposure
* uncertainty scoring
* cross-dataset evidence synthesis
* provenance-aware reporting

## Status

Prototype under active development for the Hypothesis2Omics codeathon.
