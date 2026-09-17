"""Recover cohort structure that lives only in GEO sample titles.

    python cohorts.py        # worked example on GSE13485 and self-tests

Depends on nothing but the standard library.


Why this exists
---------------

An accession is not always one study. GSE13485's own Overall design field
says:

    "This study is a combination of 2 separate YF-17D vaccination trials.
     Trial 1 has 15 subjects, and Trial 2 has 10 subjects. The two trials
     occurred 1 year apart and used different lots of the vaccine."

Querec's "independent, blinded trial", the one carrying the 90% accuracy
figure, is Trial 2. So the paper's own validation design is fit on Trial 1 and
test on Trial 2, and treating GSE13485 as a single unit throws that away.

Nothing in the pipeline currently recovers this, because the structured
metadata does not carry it. The sample titles do:

    Trial1 Subject ID 1901 Day 0
    Trial2 Subject ID 2205 Day 7

This module reads subject, timepoint and cohort out of those strings and
builds analysis units from them, which is what `synthesis.py` groups on.


What it does not do
-------------------

It does not invent cohorts. If no pattern matches, `cohort` is None and the
samples form one unit, which is the honest default. A split has to be visible
in the data to be made.

Parsing a cohort label does not, by itself, establish independent replication.
Callers must provide an explicit cohort-to-group mapping before separate cohorts
count as separate independence groups. Without that mapping, every cohort in an
accession remains in the same accession-level independence group.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Iterable, Optional, Pattern

__all__ = [
    "COHORT_PATTERNS",
    "DAY_PATTERNS",
    "HOUR_PATTERNS",
    "SUBJECT_PATTERNS",
    "TIMEPOINT_PATTERNS",
    "analysis_units_from_samples",
    "parse_sample_title",
    "summarize_cohorts",
]


# ---------------------------------------------------------------------------
# Patterns
# ---------------------------------------------------------------------------
#
# Ordered, first match wins. Each must expose exactly one capturing group.
# Add to these rather than editing them, so a new series cannot silently
# change how an existing one parses.

COHORT_PATTERNS: tuple[tuple[str, Pattern[str]], ...] = (
    # ("label shown in output", pattern). The label comes from whichever
    # pattern matched, so "Trial1" reads back as "Trial 1" and not "Cohort 1".
    ("Trial",  re.compile(r"\btrial[\s_-]*(\d+)\b", re.I)),
    ("Cohort", re.compile(r"\bcohort[\s_-]*([A-Za-z0-9]+)\b", re.I)),
    ("Batch",  re.compile(r"\bbatch[\s_-]*([A-Za-z0-9]+)\b", re.I)),
    ("Study",  re.compile(r"\bstudy[\s_-]*(\d+)\b", re.I)),
)

SUBJECT_PATTERNS: tuple[Pattern[str], ...] = (
    # "Subject ID 1901", "subject id: 1901"
    re.compile(r"\bsubject[\s_-]*id[\s_:-]*([A-Za-z0-9]+)\b", re.I),
    # "Subject 1901", "Donor 138", "Participant P07"
    re.compile(r"\b(?:subject|donor|participant)[\s_-]*([A-Za-z0-9]+)\b", re.I),
)

DAY_PATTERNS: tuple[Pattern[str], ...] = (
    # A minus counts as a sign only when something separates it from "day".
    # "day -14" is fourteen days before; "day-14" is day fourteen written with
    # a hyphen. Getting this backwards silently flips pre and post vaccination.
    re.compile(r"\bday[\s_]+(-\d+)\b", re.I),
    re.compile(r"\bday[\s_-]*(\d+)\b", re.I),
    re.compile(r"\bd(\d+)\b", re.I),
)

HOUR_PATTERNS: tuple[Pattern[str], ...] = (
    # "time -3hr" is negative; "time-3hr" uses a hyphen as a separator.
    # Keep signed patterns first so a later unsigned pattern cannot eat the sign.
    re.compile(r"\btime[\s_]+(-\d+)\s*(?:hr|hours?|h)\b", re.I),
    re.compile(r"(?<!\w)(-\d+)\s*(?:hr|hours?|h)\b", re.I),
    re.compile(r"\btime[\s_-]*(\d+)\s*(?:hr|hours?|h)\b", re.I),
    re.compile(r"\b(\d+)\s*(?:hr|hours?|h)\b", re.I),
)

TIMEPOINT_PATTERNS: tuple[Pattern[str], ...] = DAY_PATTERNS + HOUR_PATTERNS


def _first_match(text: str, patterns: Iterable[Pattern[str]]) -> Optional[str]:
    for pat in patterns:
        m = pat.search(text)
        if m:
            return m.group(1)
    return None


def parse_sample_title(title: str) -> dict[str, Optional[str]]:
    """Pull cohort, subject and timepoint out of one GEO sample title.

    Every field is optional. A field that does not match is None rather than
    a guess.

        >>> parsed = parse_sample_title("Trial1 Subject ID 1901 Day 0")
        >>> parsed["cohort"], parsed["subject"], parsed["timepoint"]
        ('Trial 1', '1901', 'Day 0')
    """
    text = (title or "").strip()

    cohort = None
    for label, pat in COHORT_PATTERNS:
        m = pat.search(text)
        if m:
            cohort = f"{label} {m.group(1)}"
            break

    subject = _first_match(text, SUBJECT_PATTERNS)
    day = None
    hour = None

    for pat in DAY_PATTERNS:
        m = pat.search(text)
        if m:
            day = m.group(1)
            break
    if day is None:
        for pat in HOUR_PATTERNS:
            m = pat.search(text)
            if m:
                hour = m.group(1)
                break

    if day is not None:
        timepoint = f"Day {day}"
    elif hour is not None:
        timepoint = f"Hour {hour}"
    else:
        timepoint = None

    return {
        "title": text or None,
        "cohort": cohort,
        "subject": subject,
        "timepoint": timepoint,
    }


# ---------------------------------------------------------------------------
# Building analysis units
# ---------------------------------------------------------------------------


def analysis_units_from_samples(
    samples: Iterable[dict[str, Any]],
    gse_id: str,
    platform: Optional[str] = None,
    min_subjects: int = 2,
    independence_groups: Optional[dict[str, str]] = None,
) -> list[dict[str, Any]]:
    """Group parsed samples into analysis units, one per cohort.

    `samples` are dicts carrying at least `title`, optionally `gsm`. Each unit
    comes back shaped for `synthesis.synthesize`, carrying `analysis_unit_id`,
    `gse_id`, `cohort`, `platform`, `independence_group`, plus the subjects and
    sample ids behind it.

    A cohort with fewer than `min_subjects` subjects is not treated as its own
    unit. A single person is not a replication cohort, and letting a stray
    label create one would manufacture structure rather than find it.

    `independence_groups` is an explicit cohort label -> independence group
    mapping. It must cover every usable cohort or be omitted. Without it,
    cohort-specific units remain in the same accession-level group. This keeps
    technical batches and treatment arms from silently becoming independent
    biological replications.

    When at least two usable cohorts exist, samples without a usable cohort
    label are excluded from cohort-specific units and the exclusion is recorded
    in each unit's note. They are never emitted as a third replication group.
    """
    independence_groups = dict(independence_groups or {})
    parsed = []
    for s in samples:
        rec = parse_sample_title(s.get("title", ""))
        rec["gsm"] = s.get("gsm")
        parsed.append(rec)

    by_cohort: dict[Optional[str], list[dict[str, Any]]] = defaultdict(list)
    for rec in parsed:
        by_cohort[rec["cohort"]].append(rec)

    named = {c: rows for c, rows in by_cohort.items() if c is not None}

    # Cohorts too small to stand alone get folded back in.
    too_small = {
        c for c, rows in named.items()
        if len({r["subject"] for r in rows if r["subject"]}) < min_subjects
    }
    usable = {c: rows for c, rows in named.items() if c not in too_small}

    if independence_groups:
        declared = set(independence_groups)
        observed = set(usable)
        if declared != observed:
            missing = sorted(observed - declared)
            unknown = sorted(declared - observed)
            raise ValueError(
                "independence_groups must cover every usable cohort exactly; "
                f"missing={missing}, unknown={unknown}"
            )
        if any(not str(group).strip() for group in independence_groups.values()):
            raise ValueError("independence group ids must be non-empty strings")

    folded: list[dict[str, Any]] = list(by_cohort.get(None, []))
    for c in too_small:
        folded.extend(named[c])

    units: list[dict[str, Any]] = []

    if len(usable) >= 2:
        unresolved_note = None
        if folded:
            unresolved_note = (
                f"{len(folded)} sample(s) had no usable cohort label and were "
                "excluded from cohort-specific units"
            )
        for cohort in sorted(usable):
            rows = usable[cohort]
            group = independence_groups.get(cohort, f"accession:{gse_id}")
            notes = []
            if not independence_groups:
                notes.append(
                    "cohort parsed, but independence was not explicitly declared; "
                    "counted at accession level"
                )
            if unresolved_note:
                notes.append(unresolved_note)
            units.append(
                _unit(
                    gse_id,
                    cohort,
                    rows,
                    platform,
                    independence_group=group,
                    folded_note="; ".join(notes) or None,
                )
            )
    else:
        # Fewer than two usable cohorts means no defensible split.
        reasons = []
        if len(usable) == 1:
            reasons.append(
                f"only one usable cohort label ({next(iter(usable))}), "
                "which is not a split"
            )
        if too_small:
            reasons.append(
                f"cohort label(s) {', '.join(sorted(too_small))} had fewer than "
                f"{min_subjects} subjects and were folded in"
            )
        if not named:
            reasons.append("no cohort labels found in the sample titles")
        note = "; ".join(reasons) or None
        units.append(
            _unit(
                gse_id,
                None,
                parsed,
                platform,
                independence_group=f"accession:{gse_id}",
                folded_note=note,
            )
        )

    return units


def _unit(
    gse_id: str,
    cohort: Optional[str],
    rows: list[dict[str, Any]],
    platform: Optional[str],
    independence_group: str,
    folded_note: Optional[str],
) -> dict[str, Any]:
    subjects = sorted({r["subject"] for r in rows if r["subject"]})
    timepoints = sorted(
        {r["timepoint"] for r in rows if r["timepoint"]},
        key=_timepoint_sort_key,
    )
    slug = re.sub(r"[^A-Za-z0-9]+", "", cohort).lower() if cohort else "all"
    unit = {
        "analysis_unit_id": f"{gse_id}_{slug}",
        "gse_id": gse_id,
        "cohort": cohort,
        "platform": platform,
        "independence_group": independence_group,
        "n_samples": len(rows),
        "n_subjects": len(subjects),
        "subjects": subjects,
        "timepoints": timepoints,
        "gsms": sorted(r["gsm"] for r in rows if r.get("gsm")),
    }
    if folded_note:
        unit["note"] = folded_note
    return unit


def _timepoint_sort_key(tp: str) -> tuple[int, float]:
    m = re.search(r"(-?\d+)", tp)
    n = float(m.group(1)) if m else 0.0
    return (0 if tp.lower().startswith("day") else 1, n)


def summarize_cohorts(units: list[dict[str, Any]]) -> str:
    """One readable block, for a report or a Slack message."""
    lines = []
    for u in units:
        label = u["cohort"] or "no cohort label"
        lines.append(
            f"{u['analysis_unit_id']}: {label}, {u['n_subjects']} subject(s), "
            f"{u['n_samples']} sample(s)"
        )
        if u["timepoints"]:
            lines.append(f"    timepoints: {', '.join(u['timepoints'])}")
        if u.get("note"):
            lines.append(f"    note: {u['note']}")
    lines.append(
        f"-> {len(units)} analysis unit(s) from "
        f"{sum(u['n_samples'] for u in units)} sample(s)"
    )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Worked example
# ---------------------------------------------------------------------------


def _gse13485_samples() -> list[dict[str, Any]]:
    """Titles in GSE13485's documented format.

    The format is taken from the GEO record. Subject numbers below are
    illustrative stand-ins for the real 25, which is enough to exercise the
    parser. Trial 1 has 15 subjects and Trial 2 has 10, as GEO states.
    """
    rows = []
    for i in range(15):
        for day in (0, 1, 3, 7, 21):
            rows.append({
                "gsm": f"GSM3400{len(rows):02d}",
                "title": f"Trial1 Subject ID {1901 + i} Day {day}",
            })
    for i in range(10):
        for day in (0, 1, 3, 7, 21):
            rows.append({
                "gsm": f"GSM3400{len(rows):02d}",
                "title": f"Trial2 Subject ID {2201 + i} Day {day}",
            })
    return rows


def _demo() -> None:
    print("=" * 74)
    print("Recovering GSE13485's two trials from sample titles.")
    print("=" * 74)

    samples = _gse13485_samples()
    print(f"\nInput: {len(samples)} sample titles, for example")
    for s in samples[:2] + samples[-2:]:
        print(f"  {s['gsm']}  {s['title']}")

    print("\nParsed (conservative default):")
    parsed_units = analysis_units_from_samples(
        samples, "GSE13485", platform="GPL7567"
    )
    print(summarize_cohorts(parsed_units))
    print("  Both cohorts remain one independence group until explicitly curated.")

    # GEO's Overall design documents two trials run a year apart with different
    # vaccine lots. That external design evidence, not the title parser alone,
    # authorizes two independence groups for this accession.
    independence_groups = {
        "Trial 1": "GSE13485|trial1",
        "Trial 2": "GSE13485|trial2",
    }
    units = analysis_units_from_samples(
        samples,
        "GSE13485",
        platform="GPL7567",
        independence_groups=independence_groups,
    )

    print("\nThis is what changes downstream:")
    print("  without the split -> 1 group, cannot meet a two group bar")
    print("  with the split    -> 2 groups, and it is the paper's own")
    print("                       train and test design, not one we invented")

    print("\n" + "-" * 74)
    print("A series with no cohort labels stays one unit, as it should.")
    print("-" * 74)
    plain = [
        {"gsm": f"GSM90{i:02d}", "title": f"Subject {i} day 7"}
        for i in range(8)
    ]
    print(summarize_cohorts(
        analysis_units_from_samples(plain, "GSE13699", platform="GPL10558")
    ))

    print("\nHanding these to synthesis:")
    print("  groups = {'Trial 1': 'GSE13485|trial1',")
    print("            'Trial 2': 'GSE13485|trial2'}")
    print("  units = analysis_units_from_samples(")
    print("      samples, 'GSE13485', 'GPL7567', independence_groups=groups)")
    print("  doc   = synthesize([{**u, **results[u['analysis_unit_id']]}")
    print("                      for u in units], rule)")


def _self_test() -> None:
    p = parse_sample_title("Trial1 Subject ID 1901 Day 0")
    assert p["cohort"] == "Trial 1", p
    assert p["subject"] == "1901", p
    assert p["timepoint"] == "Day 0", p

    p = parse_sample_title("Trial 2 Subject ID 2205 Day 21")
    assert p["cohort"] == "Trial 2"
    assert p["subject"] == "2205"
    assert p["timepoint"] == "Day 21"

    # GSE13484's format, hours rather than days, no trial label.
    p = parse_sample_title("Subject ID 138 time 3hr YF-17D stimulation")
    assert p["cohort"] is None
    assert p["subject"] == "138"
    assert p["timepoint"] == "Hour 3"

    # Other shapes that turn up in GEO.
    assert parse_sample_title("Cohort A donor 12 D7")["cohort"] == "Cohort A"
    assert parse_sample_title("Cohort A donor 12 D7")["subject"] == "12"
    assert parse_sample_title("batch_3 participant P07 day-14")["cohort"] == "Batch 3"
    # "day-14" is day 14 with a hyphen separator, not negative fourteen.
    # A genuine pre-vaccination timepoint is written with a space.
    assert parse_sample_title("batch_3 participant P07 day-14")["timepoint"] == "Day 14"
    assert parse_sample_title("Subject 4 day -14")["timepoint"] == "Day -14"
    assert parse_sample_title("Subject 4 time -3hr")["timepoint"] == "Hour -3"
    assert parse_sample_title("Subject 4 -3 hr")["timepoint"] == "Hour -3"
    assert parse_sample_title("Subject 4 time-3hr")["timepoint"] == "Hour 3"

    # Nothing to find is None, never a guess.
    blank = parse_sample_title("")
    assert blank["cohort"] is None and blank["subject"] is None
    assert blank["timepoint"] is None
    assert parse_sample_title("PBMC replicate 2")["cohort"] is None

    # Parsing finds two units but does not claim independent replication.
    parsed_units = analysis_units_from_samples(
        _gse13485_samples(), "GSE13485", "GPL7567"
    )
    assert len(parsed_units) == 2
    assert len({u["independence_group"] for u in parsed_units}) == 1
    assert all("not explicitly declared" in u["note"] for u in parsed_units)

    # GSE13485's external study design explicitly authorizes the trial split.
    declared_groups = {
        "Trial 1": "GSE13485|trial1",
        "Trial 2": "GSE13485|trial2",
    }
    units = analysis_units_from_samples(
        _gse13485_samples(),
        "GSE13485",
        "GPL7567",
        independence_groups=declared_groups,
    )
    assert len(units) == 2, units
    by_cohort = {u["cohort"]: u for u in units}
    assert by_cohort["Trial 1"]["n_subjects"] == 15
    assert by_cohort["Trial 2"]["n_subjects"] == 10
    assert by_cohort["Trial 1"]["n_samples"] == 75
    assert by_cohort["Trial 1"]["timepoints"] == [
        "Day 0", "Day 1", "Day 3", "Day 7", "Day 21"
    ]
    assert by_cohort["Trial 1"]["independence_group"] != \
        by_cohort["Trial 2"]["independence_group"]
    assert by_cohort["Trial 1"]["platform"] == "GPL7567"

    # A partial declaration is ambiguous and rejected rather than guessed.
    try:
        analysis_units_from_samples(
            _gse13485_samples(),
            "GSE13485",
            independence_groups={"Trial 1": "GSE13485|trial1"},
        )
        raise AssertionError("partial independence mapping should fail")
    except ValueError as exc:
        assert "cover every usable cohort" in str(exc)

    # No labels means one unit, not a manufactured split.
    plain = [{"gsm": f"G{i}", "title": f"Subject {i} day 7"} for i in range(8)]
    one = analysis_units_from_samples(plain, "GSE13699")
    assert len(one) == 1
    assert one[0]["cohort"] is None
    assert one[0]["n_subjects"] == 8

    # One label is not a split either.
    single = [{"gsm": f"G{i}", "title": f"Trial1 Subject {i} day 7"} for i in range(6)]
    assert len(analysis_units_from_samples(single, "GSE1")) == 1
    assert "not a split" in analysis_units_from_samples(single, "GSE1")[0]["note"]

    # A cohort of one subject cannot stand alone.
    lopsided = (
        [{"gsm": f"A{i}", "title": f"Trial1 Subject {i} day 7"} for i in range(6)]
        + [{"gsm": "B1", "title": "Trial2 Subject 99 day 7"}]
    )
    lop = analysis_units_from_samples(lopsided, "GSE2")
    assert len(lop) == 1, lop
    assert "fewer than 2 subjects" in lop[0]["note"]

    # Batch labels are parsed but never promoted to independent replication.
    batches = []
    for batch in (1, 2):
        for subject in (1, 2):
            batches.append(
                {
                    "gsm": f"B{batch}{subject}",
                    "title": f"Batch {batch} Subject {batch}{subject} Day 0",
                }
            )
    batch_units = analysis_units_from_samples(batches, "GSE3")
    assert len(batch_units) == 2
    assert len({u["independence_group"] for u in batch_units}) == 1

    # Unresolved samples are recorded but cannot become a third group.
    with_unlabelled = _gse13485_samples() + [
        {"gsm": "GSM_UNKNOWN", "title": "Subject X Day 0"}
    ]
    resolved = analysis_units_from_samples(
        with_unlabelled,
        "GSE13485",
        independence_groups=declared_groups,
    )
    assert len(resolved) == 2
    assert all("excluded from cohort-specific units" in u["note"] for u in resolved)

    # Units are shaped for synthesize().
    for u in units:
        assert set(("analysis_unit_id", "gse_id", "cohort", "platform",
                    "independence_group")) <= set(u)

    print("\nself-tests: all assertions passed")


if __name__ == "__main__":
    _demo()
    _self_test()

