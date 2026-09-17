# Hypothesis2Omics Scientific Validator

A scientific eligibility layer for **Hypothesis2Omics** that determines whether a retrieved omics dataset contains the biological information required to test a hypothesis before downstream statistical analysis.

## Why this is needed

Finding a dataset in ImmPort, GEO, or another repository does not automatically mean that the dataset can answer the biological question.

A candidate study may be missing:

* the required predictor or molecular feature
* the required biological outcome
* appropriate predictor/outcome timepoints
* participant-level linkage
* accessible quantitative measurements
* sufficient evidence to support analysis

The scientific validator separates **scientific eligibility** from dataset discovery, retrieval, parsing, and downstream workflow execution.

## Where it fits

```text
Natural-language hypothesis
        |
        v
Test specification
        |
        v
Dataset discovery / retrieval
        |
        v
Dataset parsing and normalization
        |
        v
Canonical validator input bundle
        |
        v
Scientific eligibility validator
        |
        v
eligible / uncertain / excluded
        |
        v
ready / needs scientific resolution
        |
        v
Harmonization and statistical workflow
        |
        v
Evidence synthesis / report
```

## Canonical validator input

The validator now consumes normalized parser outputs from:

```text
data/validator_input/
```

The primary input files are:

```text
sample_manifest.tsv
feature_expression.tsv
quantitative_outcome.tsv
```

Supporting provenance and manifest files are also retained in the input bundle.

### `sample_manifest.tsv`

Provides the ImmPort-derived sample structure, including information such as:

* study accession
* subject accession
* experimental sample
* assay / measurement technique
* study timepoint
* biosample
* repository
* GEO sample accession
* species

### `feature_expression.tsv`

Provides the requested GEO molecular feature measurements.

For the current YF-17D proof of concept:

```text
gene       = EIF2AK4
feature    = Hs.412102_at
dataset    = GSE13485
platform   = GPL7567
```

The current bundle contains 87 numeric EIF2AK4 expression measurements.

### `quantitative_outcome.tsv`

Provides participant-level quantitative biological outcomes.

For the current proof of concept, this contains 25 numeric day-15:

```text
Act CD8 T Cell Response
```

measurements from SDY1264.

## Current proof-of-concept hypothesis

> GCN2/EIF2AK4 activity is associated with the magnitude of the CD8+ T-cell response following YF-17D vaccination.

The current integrated test uses:

```text
ImmPort study: SDY1264
GEO dataset:   GSE13485
GEO platform:  GPL7567
Gene:          EIF2AK4
Feature ID:    Hs.412102_at
Outcome:       Act CD8 T Cell Response
```

## Scientific validation criteria

The current test specification evaluates nine required criteria:

```text
human_study
yf17d_vaccination
transcriptomics_available
eif2ak4_measurable
quantitative_cd8_response_available
appropriate_timepoints
participant_level_linkage
predictor_data_accessible
outcome_data_accessible
```

Each criterion retains:

* value
* confidence
* supporting evidence
* provenance / source

## Current validated result

Using the normalized files in `data/validator_input`, the validator currently identifies:

```text
112 SDY1264 manifest rows
87 transcriptomic samples
87 numeric EIF2AK4 expression measurements
25 quantitative CD8 outcomes
25 outcome subjects
25/25 outcome subjects linked to early transcriptomic measurements
25/25 outcome subjects with predictor measurements before the outcome
```

The resulting decision is:

```text
scientific_status = eligible
execution_status  = ready

failed_scientific_criteria = []
unresolved_criteria        = []
```

The recommended next action is:

```text
Proceed to harmonization and statistical workflow specification.
```

## Run the normalized handoff validator

From the repository root:

```bash
python scientific_validator/handoff_adapter.py \
  --input-dir data/validator_input \
  --study-accession SDY1264 \
  --gene EIF2AK4 \
  --output scientific_validator/output/SDY1264_validator_input_evidence.json
```

On Windows CMD, the same command can be run on one line:

```cmd
python scientific_validator\handoff_adapter.py --input-dir data\validator_input --study-accession SDY1264 --gene EIF2AK4 --output scientific_validator\output\SDY1264_validator_input_evidence.json
```

This produces:

```text
scientific_validator/output/SDY1264_validator_input_evidence.json
```

## Run the eligibility engine

```bash
python scientific_validator/eligibility_engine.py \
  --spec scientific_validator/examples/yf17d_test_spec.json \
  --dataset scientific_validator/output/SDY1264_validator_input_evidence.json
```

Windows CMD:

```cmd
python scientific_validator\eligibility_engine.py --spec scientific_validator\examples\yf17d_test_spec.json --dataset scientific_validator\output\SDY1264_validator_input_evidence.json
```

Expected result:

```text
scientific_status: eligible
execution_status: ready
failed_scientific_criteria: []
unresolved_criteria: []
```

## Decision model

The validator returns one of three scientific states:

```text
eligible
uncertain
excluded
```

It separately reports execution readiness.

This prevents a dataset from being treated as analysis-ready simply because it was successfully retrieved or parsed.

For example:

```text
eligible + ready
```

means the required scientific evidence is present and the dataset can proceed to downstream harmonization/statistical workflow specification.

```text
uncertain + needs_scientific_resolution
```

means one or more required scientific criteria remain unresolved.

```text
excluded
```

means at least one required criterion is explicitly not satisfied.

## Main validator components

```text
scientific_validator/
|
|-- README.md
|-- eligibility_engine.py
|-- handoff_adapter.py
|-- evidence_fusion.py
|-- geo_adapter.py
|-- geo_feature_verifier.py
|-- immport_adapter.py
|-- immport_manifest_adapter.py
|
|-- examples/
|   |-- yf17d_test_spec.json
|   |-- cross_repository_links.json
|   |-- SDY1264.json
|   |-- SDY1529.json
|   `-- synthetic_excluded.json
|
|-- output/
|   |-- SDY1264_handoff_evidence.json
|   |-- SDY1264_validator_input_evidence.json
|   `-- additional development / validation outputs
|
`-- tests/
    |-- test_eligibility_engine.py
    `-- test_handoff_adapter.py
```

## Current parser-to-validator interface

The current integrated workflow is:

```text
ImmPort / GEO retrieval
        |
        v
Repository-specific parsing
        |
        v
data/validator_input/
        |
        |-- sample_manifest.tsv
        |-- feature_expression.tsv
        `-- quantitative_outcome.tsv
        |
        v
handoff_adapter.py
        |
        v
criterion-level scientific evidence
        |
        v
eligibility_engine.py
        |
        v
eligible / uncertain / excluded
+
ready / needs scientific resolution
```

The validator therefore does not need to independently re-download or re-parse the original repository files once the canonical validator-input bundle has been produced.

## Cross-repository evidence

Earlier development components remain available for direct GEO/ImmPort evidence validation and cross-repository evidence fusion.

These include:

```text
geo_adapter.py
geo_feature_verifier.py
immport_adapter.py
immport_manifest_adapter.py
evidence_fusion.py
```

The cross-repository linkage registry is maintained in:

```text
scientific_validator/examples/cross_repository_links.json
```

For the current proof of concept, SDY1264 is linked with the corresponding GEO YF-17D expression dataset, allowing evidence from ImmPort and GEO to be combined while preserving provenance.

## Tests

Run:

```bash
python -m unittest discover -s scientific_validator/tests -v
```

Current test suite:

```text
Ran 7 tests
OK
```

The tests currently cover:

* known scientifically eligible datasets
* exclusion when a required criterion fails
* unresolved science not being promoted to eligible
* evidence provenance preservation
* canonical `data/validator_input` path resolution
* explicit input-path overrides
* backward compatibility with explicit input paths

The handoff path-resolution tests are cross-platform and work with Windows and POSIX-style path separators.

## Design principles

### Scientific evidence is explicit

A criterion is not promoted to `true` merely because related metadata exists.

Required scientific evidence must be established from the available normalized data.

### Unknown is different from false

Missing evidence is represented as unresolved rather than automatically treated as evidence against the dataset.

### Provenance is retained

Each criterion records the input source supporting the decision.

### Scientific eligibility is separate from retrieval

Successful dataset discovery or download does not imply that the dataset can test the hypothesis.

### Participant linkage matters

Predictor and outcome measurements must be linkable at the participant level when the hypothesis requires participant-level analysis.

### Cross-repository evidence can be complementary

ImmPort and GEO may contain different components of the same study. Evidence can therefore be combined when cross-repository linkage has been established.

## Current scope and limitations

The eligibility engine itself is driven by the supplied test specification.

The current normalized handoff adapter is still an MVP implementation demonstrated using the YF-17D / EIF2AK4 / CD8 test case. Some extraction rules and timing assumptions remain specific to this proof of concept.

The next step toward a general-purpose validator is to derive elements such as:

* target gene / molecular feature
* biological outcome
* predictor time window
* outcome time window
* assay requirements
* participant-linkage requirements

directly from the generated hypothesis test specification rather than encoding them in the handoff logic.

## Next integration priorities

* dynamically drive validation criteria from the generated test specification
* generalize gene / feature validation beyond EIF2AK4
* generalize quantitative outcome detection beyond CD8
* generalize predictor/outcome time windows
* connect validated evidence to harmonization and statistical workflow generation
* expand automated tests across multiple studies and hypotheses
* preserve structured provenance throughout downstream analysis and reporting

## Status

Active development for the **Hypothesis2Omics** codeathon.

Current milestone:

```text
Repository retrieval/parsing
        ->
Canonical validator-input bundle
        ->
Scientific validation
        ->
ELIGIBLE + READY
```

The parser-to-validator integration has been tested end-to-end using SDY1264 and GSE13485.
