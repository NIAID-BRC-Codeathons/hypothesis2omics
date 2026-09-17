import argparse
import json
from pathlib import Path

import pandas as pd

VALIDATOR_INPUT_FILENAMES = {
    "manifest": "sample_manifest.tsv",
    "feature_expression": "feature_expression.tsv",
    "outcome": "quantitative_outcome.tsv",
}


def resolve_input_paths(
    input_dir: str | None,
    manifest: str | None,
    feature_expression: str | None,
    outcome: str | None,
) -> tuple[str, str, str]:
    """Resolve explicit input paths, falling back to canonical bundle filenames."""
    base = Path(input_dir) if input_dir else None
    resolved = {
        "manifest": manifest,
        "feature_expression": feature_expression,
        "outcome": outcome,
    }
    if base is not None:
        for name, filename in VALIDATOR_INPUT_FILENAMES.items():
            if resolved[name] is None:
                resolved[name] = str(base / filename)
    missing = [name for name, path in resolved.items() if path is None]
    if missing:
        options = ", ".join(name.replace("_", "-") for name in missing)
        raise ValueError(f"Missing input paths for: {options}")
    return (
        str(resolved["manifest"]),
        str(resolved["feature_expression"]),
        str(resolved["outcome"]),
    )


def criterion(value, confidence, evidence, sources):
    return {
        "value": value,
        "confidence": confidence,
        "evidence": evidence,
        "sources": sources,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build scientific-validator evidence from normalized GEO/ImmPort "
            "handoff tables."
        )
    )
    parser.add_argument(
        "--input-dir",
        help="Directory containing the canonical validator-input bundle",
    )
    parser.add_argument("--manifest")
    parser.add_argument("--feature-expression")
    parser.add_argument("--outcome")
    parser.add_argument("--study-accession", required=True)
    parser.add_argument("--gene", default="EIF2AK4")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        args.manifest, args.feature_expression, args.outcome = resolve_input_paths(
            args.input_dir,
            args.manifest,
            args.feature_expression,
            args.outcome,
        )
    except ValueError as exc:
        parser.error(str(exc))

    study = args.study_accession
    gene = args.gene

    manifest = pd.read_csv(args.manifest, sep="\t")
    features = pd.read_csv(args.feature_expression, sep="\t")
    outcomes = pd.read_csv(args.outcome, sep="\t")

    manifest = manifest[manifest["study_accession"].eq(study)].copy()
    features = features[features["study_accession"].eq(study)].copy()
    outcomes = outcomes[outcomes["study_accession"].eq(study)].copy()

    # ---------- Human study ----------
    species = manifest["species"].fillna("").astype(str)
    human = species.str.contains("Homo sapiens", case=False, regex=False).any()

    # ---------- YF-17D ----------
    text_columns = [
        c for c in
        ["experiment_name", "study_time_t0_event", "arm_name"]
        if c in manifest.columns
    ]

    combined_text = (
        manifest[text_columns]
        .fillna("")
        .astype(str)
        .agg(" ".join, axis=1)
        if text_columns
        else pd.Series(dtype=str)
    )

    yf17d = combined_text.str.contains(
        r"YF[- ]?17D|YF17D|yellow fever",
        case=False,
        regex=True,
    ).any()

    # ---------- Transcriptomics ----------
    measurement = manifest["measurement_technique"].fillna("").astype(str)
    tx_mask = measurement.str.contains(
        r"transcription|transcriptomic|gene expression",
        case=False,
        regex=True,
    )

    tx = manifest[tx_mask].copy()

    tx["study_time_collected"] = pd.to_numeric(
        tx["study_time_collected"], errors="coerce"
    )

    early_tx = tx[
        tx["study_time_collected"].isin([0, 1, 3, 7])
    ].copy()

    # ---------- GEO feature ----------
    features["expression_value_numeric"] = pd.to_numeric(
        features["expression_value"], errors="coerce"
    )

    gene_features = features[
        features["gene"].fillna("").astype(str).str.upper().eq(gene.upper())
        & features["expression_value_numeric"].notna()
    ].copy()

    feature_samples = set(
        gene_features["sample_accession"].dropna().astype(str)
    )

    linked_predictor_rows = early_tx[
        early_tx["repository_accession"]
        .fillna("")
        .astype(str)
        .isin(feature_samples)
    ].copy()

    predictor_subjects = set(
        linked_predictor_rows["subject_accession"]
        .dropna()
        .astype(str)
    )

    # ---------- Quantitative CD8 outcome ----------
    outcomes["result_value_numeric"] = pd.to_numeric(
        outcomes["result_value"], errors="coerce"
    )

    outcomes["timepoint_day_numeric"] = pd.to_numeric(
        outcomes["timepoint_day"], errors="coerce"
    )

    cd8 = outcomes[
        outcomes["outcome_name"]
        .fillna("")
        .astype(str)
        .str.contains("CD8", case=False, regex=False)
        & outcomes["result_value_numeric"].notna()
    ].copy()

    outcome_subjects = set(
        cd8["subject_accession"].dropna().astype(str)
    )

    overlap_subjects = predictor_subjects & outcome_subjects

    all_outcomes_linked = (
        len(outcome_subjects) > 0
        and outcome_subjects.issubset(predictor_subjects)
    )

    # ---------- Timing check ----------
    timing_ok_subjects = set()

    for subject in overlap_subjects:
        predictor_days = linked_predictor_rows.loc[
            linked_predictor_rows["subject_accession"].astype(str).eq(subject),
            "study_time_collected",
        ].dropna()

        outcome_days = cd8.loc[
            cd8["subject_accession"].astype(str).eq(subject),
            "timepoint_day_numeric",
        ].dropna()

        if (
            len(predictor_days) > 0
            and len(outcome_days) > 0
            and predictor_days.min() < outcome_days.max()
        ):
            timing_ok_subjects.add(subject)

    all_timing_ok = (
        len(outcome_subjects) > 0
        and outcome_subjects.issubset(timing_ok_subjects)
    )

    criteria = {
        "human_study": criterion(
            "true" if human else "false",
            "high",
            f"{len(manifest)} study rows inspected; Homo sapiens detected={human}.",
            [args.manifest],
        ),

        "yf17d_vaccination": criterion(
            "true" if yf17d else "unknown",
            "high" if yf17d else "low",
            f"YF-17D/yellow-fever terminology detected={yf17d}.",
            [args.manifest],
        ),

        "transcriptomics_available": criterion(
            "true" if len(tx) > 0 else "false",
            "high",
            f"{len(tx)} transcriptomic experimental-sample rows detected.",
            [args.manifest],
        ),

        "eif2ak4_measurable": criterion(
            "true" if len(gene_features) > 0 else "false",
            "high",
            f"{len(gene_features)} numeric {gene} feature-expression measurements detected.",
            [args.feature_expression],
        ),

        "quantitative_cd8_response_available": criterion(
            "true" if len(cd8) > 0 else "false",
            "high",
            f"{len(cd8)} numeric CD8 response measurements detected.",
            [args.outcome],
        ),

        "appropriate_timepoints": criterion(
            "true" if all_timing_ok else (
                "unknown" if len(timing_ok_subjects) > 0 else "false"
            ),
            "high" if all_timing_ok else "medium",
            (
                f"{len(timing_ok_subjects)}/{len(outcome_subjects)} outcome subjects "
                f"have earlier predictor measurements."
            ),
            [args.manifest, args.feature_expression, args.outcome],
        ),

        "participant_level_linkage": criterion(
            "true" if all_outcomes_linked else (
                "unknown" if len(overlap_subjects) > 0 else "false"
            ),
            "high" if all_outcomes_linked else "medium",
            (
                f"{len(overlap_subjects)}/{len(outcome_subjects)} quantitative outcome "
                f"subjects link to early transcriptomic {gene} measurements."
            ),
            [args.manifest, args.feature_expression, args.outcome],
        ),

        "predictor_data_accessible": criterion(
            "true" if len(linked_predictor_rows) > 0 else "false",
            "high",
            (
                f"{len(linked_predictor_rows)} early transcriptomic samples have "
                f"numeric {gene} measurements."
            ),
            [args.manifest, args.feature_expression],
        ),

        "outcome_data_accessible": criterion(
            "true" if len(cd8) > 0 else "false",
            "high",
            f"{len(cd8)} directly accessible numeric CD8 outcome measurements.",
            [args.outcome],
        ),
    }

    result = {
        "dataset_id": study,
        "repository": "GEO+ImmPort normalized handoff",
        "criteria": criteria,
        "assay_context": {
            "gene": gene,
            "transcriptomic_rows": int(len(tx)),
            "feature_expression_rows": int(len(gene_features)),
            "quantitative_cd8_rows": int(len(cd8)),
            "predictor_subjects": int(len(predictor_subjects)),
            "outcome_subjects": int(len(outcome_subjects)),
            "participant_overlap": int(len(overlap_subjects)),
        },
        "data_summary": {
            "manifest_rows": int(len(manifest)),
            "early_transcriptomic_rows": int(len(early_tx)),
            "geo_feature_rows": int(len(gene_features)),
            "cd8_outcome_rows": int(len(cd8)),
            "all_outcome_subjects_linked": bool(all_outcomes_linked),
            "all_outcome_subjects_have_earlier_predictor": bool(all_timing_ok),
        },
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
