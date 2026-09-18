"""Resolve a hypothesis's predictor gene and outcome concept to real identifiers,
then produce feature_expression.tsv and quantitative_outcome.tsv.

The gap this closes: scientific_validator's two remaining validator-input tables
currently require a human to already know a study's exact GEO probe ID and exact
ImmPort outcome label, hand-typed into validator_handoff_config.json with no
recorded source. This script proposes and verifies both instead of asking someone
to already know them.

Three-stage resolution, deterministic first, LLM never involved:

    1. search_terms match   -- score platform/outcome candidates against
                                02_test_spec.json's own search_terms (already
                                generated, zero new cost). Prefers the longest
                                matching phrase, not flat word-union overlap, so a
                                specific synonym beats a generic family term.
    2. HGNC bridge           -- predictor only. search_terms carry the gene symbol
                                and pathway aliases but usually not the platform's
                                spelled-out full name (measured: none of three
                                EIF2AK4 trials generated it). Look the symbol up in
                                a vendored, pinned HGNC table (reference/) to get
                                the official full name, then match that.
    3. Needs human review    -- ties or zero candidates go to a review artifact,
                                same shape as the 01_parsed.yaml human gate. Never
                                guessed, per AGENTS.md: flag ambiguity, don't
                                silently resolve it.

Whatever resolves at stage 1 or 2 is then VERIFIED against the real data before
being accepted -- does the probe have real numeric values in the series matrix;
does the outcome panel have real numeric values with subject-level linkage. A
candidate that fails verification is rejected regardless of how it was matched.
This is the same discipline evidence_rules/run_h2.py already applies to its own
probe identity (FeatureIdentity.confirmed_by): a plausible match is not the same
as a confirmed one.

Run:
    python3 feature_outcome_resolver.py --help
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger("feature_outcome_resolver")

REFERENCE_DIR = Path(__file__).resolve().parent / "reference"
HGNC_TABLE_PATH = REFERENCE_DIR / "hgnc_protein_coding.tsv.gz"

FEATURE_EXPRESSION_COLUMNS = ["study_accession", "gene", "sample_accession", "expression_value"]
QUANTITATIVE_OUTCOME_COLUMNS = [
    "study_accession", "subject_accession", "timepoint_day", "outcome_name", "result_value",
]


# --------------------------------------------------------------------- resolution result


@dataclass
class ResolutionResult:
    """What was resolved, how, and whether it survived verification.

    status, in the order a result actually moves through them:
        "candidate"     one match found, not yet checked against real data
        "tied"          multiple equally-good matches found (candidates_considered
                        holds all of them), not yet checked against real data
        "verified"      exactly one candidate (from "candidate" or narrowed down
                        from "tied") was confirmed against real data
        "rejected"      a "candidate" was checked and did not hold up
        "needs_review"  zero candidates, or verification left zero or more than
                        one standing -- a human decides, nothing here guesses

    Only "verified" results feed the output tables.
    """

    concept: str
    status: str
    method: str
    evidence: str
    candidates_considered: list[dict[str, Any]] = field(default_factory=list)
    resolved_id: Optional[str] = None
    resolved_label: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _decide(scored: list[tuple[int, str, dict]]) -> tuple[str, Any]:
    """("single", row) | ("tied", [rows]) | ("none", None), from a scored list
    sorted descending by match length."""
    if not scored:
        return "none", None
    if len(scored) == 1 or scored[0][0] > scored[1][0]:
        return "single", scored[0][2]
    tied_len = scored[0][0]
    return "tied", [row for length, _term, row in scored if length == tied_len]


# --------------------------------------------------------------------- HGNC bridge


def load_hgnc_table(path: Path = HGNC_TABLE_PATH) -> dict[str, dict[str, Any]]:
    """symbol (and every known alias/previous symbol) -> {symbol, name, aliases}.

    Vendored, pinned snapshot (see reference/hgnc_protein_coding.provenance.json).
    No network call here or anywhere else in this module.
    """
    index: dict[str, dict[str, Any]] = {}
    with gzip.open(path, "rt", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            symbol = row["symbol"].strip()
            if not symbol:
                continue
            aliases = [a.strip() for a in (row.get("alias_symbol") or "").split("|") if a.strip()]
            prev = [a.strip() for a in (row.get("prev_symbol") or "").split("|") if a.strip()]
            prev_names = [a.strip() for a in (row.get("prev_name") or "").split("|") if a.strip()]
            entry = {
                "symbol": symbol, "name": row["name"].strip(),
                "aliases": aliases, "prev_names": prev_names,
            }
            for key in [symbol, *aliases, *prev]:
                index.setdefault(key.upper(), entry)
    return index


# --------------------------------------------------------------------- platform table


def parse_platform_table(family_soft_path: Path) -> list[dict[str, str]]:
    """Extract the !platform_table_begin/!platform_table_end section.

    Returns one dict per probe row, keyed by the table's own header (typically
    ID, Description, UniGene, GB_ACC, ... -- varies by platform, so nothing here
    assumes a fixed schema beyond ID and Description existing).
    """
    with gzip.open(family_soft_path, "rt", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    start = end = None
    for i, line in enumerate(lines):
        if line.startswith("!platform_table_begin"):
            start = i + 1
        elif line.startswith("!platform_table_end"):
            end = i
            break
    if start is None or end is None:
        raise ValueError(f"no platform table found in {family_soft_path}")

    reader = csv.DictReader(io.StringIO("".join(lines[start:end])), delimiter="\t")
    return list(reader)


# --------------------------------------------------------------------- predictor resolution


def _longest_match_score(text: str, terms: list[str]) -> tuple[int, str]:
    """Best (length, matched_term) for the longest term appearing in text.

    Longest-match-wins, not term-count, so a specific synonym ("interferon
    regulatory factor 7") outscores a generic one that also happens to match
    ("interferon regulatory factor") even though both are literal hits.
    """
    text_lower = text.lower()
    best_len, best_term = 0, ""
    for term in terms:
        t = term.strip()
        if t and t.lower() in text_lower and len(t) > best_len:
            best_len, best_term = len(t), t
    return best_len, best_term


def resolve_predictor(
    gene_symbol: str,
    search_terms: list[str],
    platform_rows: list[dict[str, str]],
    hgnc: dict[str, dict[str, Any]],
) -> ResolutionResult:
    if not platform_rows:
        return ResolutionResult(
            concept=gene_symbol, status="needs_review", method="none",
            evidence="platform table is empty",
        )
    columns = list(platform_rows[0])
    id_col = "ID" if "ID" in columns else columns[0]

    # Stage 0: an explicit symbol column (Illumina-style platforms: Symbol,
    # ILMN_Gene) makes this an exact lookup, not a text-matching problem at all.
    # Checked before search_terms, since an exact symbol match is strictly more
    # reliable than any phrase heuristic.
    symbol_col = next((c for c in ("Symbol", "ILMN_Gene", "Gene Symbol", "GENE_SYMBOL")
                       if c in columns), None)
    if symbol_col:
        hits = [r for r in platform_rows if (r.get(symbol_col) or "").strip().upper() == gene_symbol.upper()]
        if len(hits) == 1:
            row = hits[0]
            return ResolutionResult(
                concept=gene_symbol, status="candidate", method="symbol_column",
                evidence=f"exact match on platform column {symbol_col!r}",
                resolved_id=row.get(id_col), resolved_label=row.get(symbol_col),
            )
        elif len(hits) > 1:
            # Multiple probes for one gene (common: different exons/transcripts).
            # Not resolvable from this table alone -- verification against the
            # real expression data gets first crack at narrowing it; a human
            # only sees it if more than one candidate still stands afterward.
            return ResolutionResult(
                concept=gene_symbol, status="tied", method="symbol_column",
                evidence=f"{len(hits)} probes share Symbol={gene_symbol!r} on this platform",
                candidates_considered=[{"id": r.get(id_col), symbol_col: r.get(symbol_col)} for r in hits],
            )

    desc_col = next((c for c in ("Description", "Definition") if c in columns), None)
    if desc_col is None:
        return ResolutionResult(
            concept=gene_symbol, status="needs_review", method="none",
            evidence=f"platform table has no symbol or description-style column: {columns}",
        )

    def _rows_to_candidates(rows: list[dict]) -> list[dict[str, Any]]:
        return [{"id": r.get(id_col), "description": r.get(desc_col)} for r in rows]

    # Stage 1: search_terms, longest match wins.
    scored = []
    for row in platform_rows:
        desc = row.get(desc_col) or ""
        length, term = _longest_match_score(desc, search_terms)
        if length > 0:
            scored.append((length, term, row))
    scored.sort(key=lambda x: -x[0])

    outcome, payload = _decide(scored)
    if outcome == "single":
        row = payload
        return ResolutionResult(
            concept=gene_symbol, status="candidate", method="search_terms",
            evidence=f"matched {desc_col}={row.get(desc_col)!r}",
            resolved_id=row.get(id_col), resolved_label=row.get(desc_col),
        )

    # Stage 1 wasn't decisive (tied or empty) -- fall through to the HGNC bridge
    # rather than stopping here. A tie at stage 1 is often just a generic term
    # ("autophagy" matches 14 different ATG-family genes); the bridge's official
    # full name is usually more specific and can still resolve it outright.
    stage1_outcome, stage1_candidates = outcome, (payload if outcome == "tied" else [])

    hgnc_entry = hgnc.get(gene_symbol.upper())
    if hgnc_entry:
        bridged_terms = [
            hgnc_entry["name"], hgnc_entry["symbol"],
            *hgnc_entry["aliases"], *hgnc_entry["prev_names"],
        ]
        scored2 = []
        for row in platform_rows:
            desc = row.get(desc_col) or ""
            length, term = _longest_match_score(desc, bridged_terms)
            if length > 0:
                scored2.append((length, term, row))
        scored2.sort(key=lambda x: -x[0])

        outcome2, payload2 = _decide(scored2)
        if outcome2 == "single":
            row = payload2
            return ResolutionResult(
                concept=gene_symbol, status="candidate", method="hgnc_bridge",
                evidence=f"HGNC name/alias/former-name bridged to {desc_col}={row.get(desc_col)!r}",
                resolved_id=row.get(id_col), resolved_label=row.get(desc_col),
            )
        if outcome2 == "tied":
            # Prefer whichever stage produced the narrower, more informative tie.
            if stage1_outcome != "tied" or len(payload2) < len(stage1_candidates):
                stage1_outcome, stage1_candidates = "tied", payload2
                method = "hgnc_bridge"
            else:
                method = "search_terms"
            return ResolutionResult(
                concept=gene_symbol, status="tied", method=method,
                evidence=f"{len(stage1_candidates)} platform rows still tied after search_terms and HGNC bridge",
                candidates_considered=_rows_to_candidates(stage1_candidates),
            )

    if stage1_outcome == "tied":
        return ResolutionResult(
            concept=gene_symbol, status="tied", method="search_terms",
            evidence=f"{len(stage1_candidates)} platform rows tie on search_terms match length"
                    + ("" if hgnc_entry else f" (symbol {gene_symbol!r} not found in HGNC table, "
                                             "bridge unavailable)"),
            candidates_considered=_rows_to_candidates(stage1_candidates),
        )

    return ResolutionResult(
        concept=gene_symbol, status="needs_review", method="none",
        evidence="no match via search_terms or HGNC bridge"
                + ("" if hgnc_entry else f" (symbol {gene_symbol!r} not found in HGNC table)"),
    )


def verify_predictor(
    probe_id: str, series_matrix_gsms: dict[str, float]
) -> tuple[bool, str]:
    """Does this probe actually have real numeric expression values.

    series_matrix_gsms: {gsm_id: value} already extracted for this probe row.
    Mirrors geo_feature_verifier.py's "feature_in_matrix" check.
    """
    numeric = [v for v in series_matrix_gsms.values() if isinstance(v, (int, float))]
    if not numeric:
        return False, f"probe {probe_id} has no numeric values in the series matrix"
    return True, f"{len(numeric)} numeric values found for probe {probe_id}"


def finalize_predictor(
    result: ResolutionResult, matrix_lookup,
) -> ResolutionResult:
    """Run verification and produce the final status.

    matrix_lookup(probe_id) -> {gsm: value}, however the caller wants to source it
    (already-loaded dict, a function reading the series matrix lazily, etc.).

    "candidate" -> "verified" or "rejected", one probe either way.
    "tied" -> verify every tied candidate; "verified" only if exactly one survives,
    otherwise "needs_review" with each candidate's verification outcome attached,
    so a human sees which ones actually had data and which didn't.
    Anything already "needs_review" passes through unchanged.
    """
    if result.status == "candidate":
        ok, evidence = verify_predictor(result.resolved_id, matrix_lookup(result.resolved_id))
        return ResolutionResult(
            concept=result.concept, status="verified" if ok else "rejected",
            method=result.method, evidence=f"{result.evidence}; verification: {evidence}",
            resolved_id=result.resolved_id if ok else None,
            resolved_label=result.resolved_label if ok else None,
            candidates_considered=[] if ok else [{"id": result.resolved_id, "description": result.resolved_label}],
        )

    if result.status == "tied":
        verified = []
        checked = []
        for cand in result.candidates_considered:
            probe_id = cand.get("id")
            ok, evidence = verify_predictor(probe_id, matrix_lookup(probe_id))
            checked.append({**cand, "verified": ok, "verification_evidence": evidence})
            if ok:
                verified.append(cand)
        if len(verified) == 1:
            row = verified[0]
            return ResolutionResult(
                concept=result.concept, status="verified", method=f"{result.method}+verification",
                evidence=f"{result.evidence}; only {row.get('id')} had real expression data",
                resolved_id=row.get("id"), resolved_label=row.get("description") or row.get("Symbol"),
            )
        return ResolutionResult(
            concept=result.concept, status="needs_review", method=result.method,
            evidence=f"{result.evidence}; {len(verified)}/{len(checked)} candidates verified"
                    + (" -- more than one real candidate, a human call" if len(verified) > 1
                       else " -- none had real data"),
            candidates_considered=checked,
        )

    return result


# --------------------------------------------------------------------- outcome resolution


def _read_zip_table(zf: zipfile.ZipFile, filename: str) -> list[dict[str, str]] | None:
    matches = [n for n in zf.namelist() if n.rsplit("/", 1)[-1] == filename]
    if not matches:
        return None
    with zf.open(matches[0]) as f:
        content = f.read().decode("utf-8", errors="replace")
    return list(csv.DictReader(io.StringIO(content), delimiter="\t"))


def resolve_outcome(
    outcome_terms: list[str], tab_zip_path: Path
) -> ResolutionResult:
    with zipfile.ZipFile(tab_zip_path) as zf:
        panels = _read_zip_table(zf, "lab_test_panel.txt")
    if not panels:
        return ResolutionResult(
            concept="outcome", status="needs_review", method="none",
            evidence="no lab_test_panel.txt in Tab.zip",
        )

    def _row_candidate(r: dict) -> dict[str, Any]:
        return {"panel_accession": r.get("LAB_TEST_PANEL_ACCESSION"),
                "name_reported": r.get("NAME_REPORTED") or r.get("NAME_PREFERRED")}

    scored = []
    for row in panels:
        name = row.get("NAME_REPORTED") or row.get("NAME_PREFERRED") or ""
        length, term = _longest_match_score(name, outcome_terms)
        if length > 0:
            scored.append((length, term, row))
    scored.sort(key=lambda x: -x[0])

    outcome, payload = _decide(scored)
    if outcome == "single":
        row = payload
        name = row.get("NAME_REPORTED") or row.get("NAME_PREFERRED")
        return ResolutionResult(
            concept="outcome", status="candidate", method="search_terms",
            evidence=f"matched lab_test_panel NAME_REPORTED={name!r}",
            resolved_id=row.get("LAB_TEST_PANEL_ACCESSION"), resolved_label=name,
        )
    if outcome == "tied":
        return ResolutionResult(
            concept="outcome", status="tied", method="search_terms",
            evidence=f"{len(payload)} lab_test_panel rows tie on search_terms match length",
            candidates_considered=[_row_candidate(r) for r in payload],
        )

    return ResolutionResult(
        concept="outcome", status="needs_review", method="none",
        evidence="no match among lab_test_panel.txt rows via search_terms",
        candidates_considered=[_row_candidate(r) for r in panels[:10]],
    )


def verify_outcome(
    panel_accession: str, lab_test_rows: list[dict[str, str]], sample_manifest_subjects: set[str]
) -> tuple[bool, str]:
    """Does this panel have real numeric values, linkable to subjects we know about."""
    matching = [r for r in lab_test_rows if r.get("LAB_TEST_PANEL_ACCESSION") == panel_accession]
    numeric = []
    for r in matching:
        v = (r.get("RESULT_VALUE_PREFERRED") or r.get("RESULT_VALUE_REPORTED") or "").strip()
        try:
            float(v)
            numeric.append(r)
        except ValueError:
            continue
    if not numeric:
        return False, f"panel {panel_accession} has no numeric result values"
    return True, f"{len(numeric)} numeric outcome values found for panel {panel_accession}"


def finalize_outcome(
    result: ResolutionResult, lab_test_rows: list[dict[str, str]], sample_manifest_subjects: set[str],
) -> ResolutionResult:
    """Same shape as finalize_predictor, for the outcome side."""
    if result.status == "candidate":
        ok, evidence = verify_outcome(result.resolved_id, lab_test_rows, sample_manifest_subjects)
        return ResolutionResult(
            concept=result.concept, status="verified" if ok else "rejected",
            method=result.method, evidence=f"{result.evidence}; verification: {evidence}",
            resolved_id=result.resolved_id if ok else None,
            resolved_label=result.resolved_label if ok else None,
            candidates_considered=[] if ok else [
                {"panel_accession": result.resolved_id, "name_reported": result.resolved_label}
            ],
        )

    if result.status == "tied":
        verified = []
        checked = []
        for cand in result.candidates_considered:
            panel_accession = cand.get("panel_accession")
            ok, evidence = verify_outcome(panel_accession, lab_test_rows, sample_manifest_subjects)
            checked.append({**cand, "verified": ok, "verification_evidence": evidence})
            if ok:
                verified.append(cand)
        if len(verified) == 1:
            row = verified[0]
            return ResolutionResult(
                concept=result.concept, status="verified", method=f"{result.method}+verification",
                evidence=f"{result.evidence}; only {row.get('panel_accession')} had real outcome data",
                resolved_id=row.get("panel_accession"), resolved_label=row.get("name_reported"),
            )
        return ResolutionResult(
            concept=result.concept, status="needs_review", method=result.method,
            evidence=f"{result.evidence}; {len(verified)}/{len(checked)} candidates verified"
                    + (" -- more than one real candidate, a human call" if len(verified) > 1
                       else " -- none had real data"),
            candidates_considered=checked,
        )

    return result


# --------------------------------------------------------------------- human review artifact


REVIEW_FILENAME = "needs_review.json"


def load_confirmed_reviews(review_path: Path) -> dict[str, dict[str, Any]]:
    """concept -> {confirmed_id, confirmed_label}, for whatever a human already
    filled in. Missing file or nothing confirmed yet both return {} quietly --
    absence of review input is a normal state, not an error."""
    if not review_path.exists():
        return {}
    try:
        doc = json.loads(review_path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    out = {}
    for item in doc.get("pending", []):
        if item.get("confirmed_id"):
            out[item["concept"]] = {
                "confirmed_id": item["confirmed_id"],
                "confirmed_label": item.get("confirmed_label") or item["confirmed_id"],
                "confirmed_by": item.get("confirmed_by"),
            }
    return out


def write_review_artifact(
    study_accession: str, results: list[tuple[str, ResolutionResult]], review_path: Path,
) -> Optional[Path]:
    """Write needs_review.json for whatever didn't resolve to 'verified'.

    Same role as the 01_parsed.yaml human gate: a person reads this, fills in
    confirmed_id/confirmed_label/confirmed_by for each pending item, and a rerun
    picks the confirmed values up via load_confirmed_reviews -- still run through
    verification, not trusted blindly just because a human typed it.

    Returns None (and writes nothing) if everything already resolved.
    """
    pending = [
        {
            "concept": result.concept,
            "type": kind,
            "status": result.status,
            "method": result.method,
            "evidence": result.evidence,
            "candidates_considered": result.candidates_considered,
            "confirmed_id": None,
            "confirmed_label": None,
            "confirmed_by": None,
        }
        for kind, result in results
        if result.status != "verified"
    ]
    if not pending:
        if review_path.exists():
            review_path.unlink()
        return None

    review_path.parent.mkdir(parents=True, exist_ok=True)
    review_path.write_text(json.dumps(
        {"study_accession": study_accession, "pending": pending}, indent=2,
    ) + "\n")
    return review_path


# --------------------------------------------------------------------- table builders


def build_feature_expression_rows(
    study_accession: str, gene: str, gsm_values: dict[str, float], study_gsms: set[str],
) -> list[dict[str, Any]]:
    """One row per GEO sample that both has an expression value and is actually
    linked to this study in sample_manifest.tsv -- not every GSM in the matrix
    necessarily belongs to this ImmPort study."""
    return [
        {"study_accession": study_accession, "gene": gene, "sample_accession": gsm,
         "expression_value": value}
        for gsm, value in gsm_values.items()
        if gsm in study_gsms
    ]


def load_biosample_linkage(tab_zip_path: Path) -> dict[str, tuple[str, str]]:
    """biosample_accession -> (subject_accession, study_time_collected), from
    biosample.txt directly.

    Deliberately NOT sourced from sample_manifest.tsv: that manifest is scoped to
    experimental samples (things that went through an assay/experiment, which is
    what feeds GEO). A lab_test outcome measurement's biosample often never went
    through an experiment at all -- it's a separate ImmPort stream -- so it can be
    entirely absent from the experimental-sample manifest while still being a
    perfectly real biosample for this study. biosample.txt covers all of them.
    """
    with zipfile.ZipFile(tab_zip_path) as zf:
        rows = _read_zip_table(zf, "biosample.txt") or []
    return {
        r["BIOSAMPLE_ACCESSION"]: (r["SUBJECT_ACCESSION"], r["STUDY_TIME_COLLECTED"])
        for r in rows if r.get("BIOSAMPLE_ACCESSION") and r.get("SUBJECT_ACCESSION")
    }


def build_quantitative_outcome_rows(
    study_accession: str,
    outcome_name: str,
    panel_accession: str,
    lab_test_rows: list[dict[str, str]],
    biosample_to_subject_day: dict[str, tuple[str, str]],
) -> list[dict[str, Any]]:
    """One row per subject with a numeric result for this panel, linked via
    biosample_accession -> subject_accession/timepoint through sample_manifest.tsv.
    A lab_test row whose biosample isn't in the manifest is dropped, not guessed."""
    rows = []
    for r in lab_test_rows:
        if r.get("LAB_TEST_PANEL_ACCESSION") != panel_accession:
            continue
        v = (r.get("RESULT_VALUE_PREFERRED") or r.get("RESULT_VALUE_REPORTED") or "").strip()
        try:
            value = float(v)
        except ValueError:
            continue
        linkage = biosample_to_subject_day.get(r.get("BIOSAMPLE_ACCESSION", ""))
        if not linkage:
            continue
        subject, day = linkage
        rows.append({
            "study_accession": study_accession, "subject_accession": subject,
            "timepoint_day": day, "outcome_name": outcome_name, "result_value": value,
        })
    return rows


def _write_tsv(rows: list[dict[str, Any]], columns: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


# --------------------------------------------------------------------- orchestration


def _matrix_lookup_factory(series_matrix_path: Path):
    def lookup(probe_id: str) -> dict[str, float]:
        with gzip.open(series_matrix_path, "rt") as f:
            header = next(f).rstrip("\n").split("\t")
            for line in f:
                if line.startswith(probe_id + "\t"):
                    vals = line.rstrip("\n").split("\t")
                    return {
                        header[i]: float(vals[i])
                        for i in range(1, len(vals)) if vals[i] not in ("", "NA")
                    }
        return {}
    return lookup


def _load_study_gsms(sample_manifest_path: Path, study_accession: str) -> set[str]:
    with open(sample_manifest_path, newline="") as f:
        return {
            r["repository_accession"]
            for r in csv.DictReader(f, delimiter="\t")
            if r["study_accession"] == study_accession and r["repository_name"] == "GEO"
            and r["repository_accession"]
        }


def resolve_and_write(
    study_accession: str,
    gene_symbol: str,
    predictor_search_terms: list[str],
    outcome_search_terms: list[str],
    family_soft_path: Path,
    series_matrix_path: Path,
    tab_zip_path: Path,
    sample_manifest_path: Path,
    output_dir: Path,
    hgnc: Optional[dict[str, dict[str, Any]]] = None,
) -> dict[str, Any]:
    """Resolve both concepts, verify, write whatever is confirmed, and leave a
    review artifact for whatever isn't. Never writes a row for anything short of
    'verified' -- an unresolved gene or outcome produces no output row and a
    named reason, not a plausible-looking guess."""
    output_dir = Path(output_dir)
    review_path = output_dir / REVIEW_FILENAME
    confirmed = load_confirmed_reviews(review_path)
    hgnc = hgnc if hgnc is not None else load_hgnc_table()

    # --- predictor ---
    if gene_symbol in confirmed:
        c = confirmed[gene_symbol]
        predictor_result = ResolutionResult(
            concept=gene_symbol, status="candidate", method="human_review",
            evidence=f"confirmed by {c.get('confirmed_by') or 'a reviewer'}",
            resolved_id=c["confirmed_id"], resolved_label=c["confirmed_label"],
        )
    else:
        platform_rows = parse_platform_table(family_soft_path)
        predictor_result = resolve_predictor(gene_symbol, predictor_search_terms, platform_rows, hgnc)
    matrix_lookup = _matrix_lookup_factory(series_matrix_path)
    predictor_final = finalize_predictor(predictor_result, matrix_lookup)

    # --- outcome ---
    if "outcome" in confirmed:
        c = confirmed["outcome"]
        outcome_result = ResolutionResult(
            concept="outcome", status="candidate", method="human_review",
            evidence=f"confirmed by {c.get('confirmed_by') or 'a reviewer'}",
            resolved_id=c["confirmed_id"], resolved_label=c["confirmed_label"],
        )
    else:
        outcome_result = resolve_outcome(outcome_search_terms, tab_zip_path)
    with zipfile.ZipFile(tab_zip_path) as zf:
        lab_test_rows = _read_zip_table(zf, "lab_test.txt") or []
    study_subjects = {s for s, _ in load_biosample_linkage(tab_zip_path).values()}
    outcome_final = finalize_outcome(outcome_result, lab_test_rows, study_subjects)

    # --- write whatever verified ---
    written = {}
    if predictor_final.status == "verified":
        gsm_values = matrix_lookup(predictor_final.resolved_id)
        study_gsms = _load_study_gsms(sample_manifest_path, study_accession)
        rows = build_feature_expression_rows(study_accession, gene_symbol, gsm_values, study_gsms)
        path = output_dir / "feature_expression.tsv"
        _write_tsv(rows, FEATURE_EXPRESSION_COLUMNS, path)
        written["feature_expression_tsv"] = {"path": str(path), "rows": len(rows)}

    if outcome_final.status == "verified":
        linkage = load_biosample_linkage(tab_zip_path)
        rows = build_quantitative_outcome_rows(
            study_accession, outcome_final.resolved_label, outcome_final.resolved_id,
            lab_test_rows, linkage,
        )
        path = output_dir / "quantitative_outcome.tsv"
        _write_tsv(rows, QUANTITATIVE_OUTCOME_COLUMNS, path)
        written["quantitative_outcome_tsv"] = {"path": str(path), "rows": len(rows)}

    review_written = write_review_artifact(
        study_accession,
        [("predictor", predictor_final), ("outcome", outcome_final)],
        review_path,
    )

    provenance = {
        "study_accession": study_accession,
        "gene_symbol": gene_symbol,
        "predictor_resolution": predictor_final.to_dict(),
        "outcome_resolution": outcome_final.to_dict(),
        "written": written,
        "review_artifact": str(review_written) if review_written else None,
    }
    (output_dir / "resolution.provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


# --------------------------------------------------------------------- CLI


def _predictor_outcome_terms_from_spec(test_spec_path: Path) -> tuple[list[str], list[str]]:
    """Pull predictor.search_terms / outcome.search_terms straight out of an
    already-generated 02_test_spec.json -- zero new cost, reuses what step 2
    already produced rather than asking for terms a second time."""
    spec = json.loads(Path(test_spec_path).read_text())
    predictor_terms = (spec.get("predictor") or {}).get("search_terms") or []
    outcome_terms = (spec.get("outcome") or {}).get("search_terms") or []
    return predictor_terms, outcome_terms


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Resolve a hypothesis's predictor gene and outcome concept to real "
                    "identifiers, then produce feature_expression.tsv and "
                    "quantitative_outcome.tsv. Deterministic; escalates to a human-review "
                    "file rather than guessing on genuine ambiguity."
    )
    parser.add_argument("--study-accession", required=True, help="e.g. SDY1264")
    parser.add_argument("--gene-symbol", required=True,
                        help="the predictor's gene symbol, e.g. EIF2AK4. Explicit, not "
                             "inferred from the test spec's free-text predictor field.")
    parser.add_argument("--test-spec", required=True, type=Path,
                        help="path to 02_test_spec.json; predictor/outcome search_terms "
                             "are read from it directly")
    parser.add_argument("--family-soft", required=True, type=Path,
                        help="the GEO platform's family SOFT file (.soft.gz)")
    parser.add_argument("--series-matrix", required=True, type=Path,
                        help="the parsed expression matrix (expression.tsv.gz from "
                             "geo_matrix_parse_module.py)")
    parser.add_argument("--tab-zip", required=True, type=Path,
                        help="the study's ImmPort *_Tab.zip")
    parser.add_argument("--sample-manifest", required=True, type=Path,
                        help="sample_manifest.tsv covering this study")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = _build_parser().parse_args()

    predictor_terms, outcome_terms = _predictor_outcome_terms_from_spec(args.test_spec)
    if not predictor_terms:
        logger.warning("no predictor.search_terms found in %s", args.test_spec)
    if not outcome_terms:
        logger.warning("no outcome.search_terms found in %s", args.test_spec)

    result = resolve_and_write(
        study_accession=args.study_accession,
        gene_symbol=args.gene_symbol,
        predictor_search_terms=predictor_terms,
        outcome_search_terms=outcome_terms,
        family_soft_path=args.family_soft,
        series_matrix_path=args.series_matrix,
        tab_zip_path=args.tab_zip,
        sample_manifest_path=args.sample_manifest,
        output_dir=args.output_dir,
    )

    pred, out = result["predictor_resolution"], result["outcome_resolution"]
    logger.info("predictor %s: %s (%s) -> %s", args.gene_symbol, pred["status"], pred["method"], pred["resolved_id"])
    logger.info("outcome: %s (%s) -> %s", out["status"], out["method"], out["resolved_id"])
    for name, info in result["written"].items():
        logger.info("wrote %s rows to %s", info["rows"], info["path"])
    if result["review_artifact"]:
        logger.info("needs review: %s", result["review_artifact"])
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
