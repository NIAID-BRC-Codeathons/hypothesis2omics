"""Turn a limma topTable into analysis units, without letting it pick winners.

    python limma_adapter.py <topTable.tsv>                  # checks only
    python limma_adapter.py <topTable.tsv> --feature ILMN_1802380
    python limma_adapter.py                                 # self-tests

Standard library only. Reads the TSV limma writes and hands rows to
`synthesis.synthesize`.


The trap this exists to close
-----------------------------

A limma table has one row per probe, often forty thousand of them. Picking the
probe after seeing the table is the cheapest p-hack available, and it does not
look like one in a log: every individual step is a legitimate lookup. So this
module makes you name the feature first. `select_feature` takes an identifier
and returns that row or raises. There is deliberately no "most significant
probe" helper. If a gene maps to several probes, declare the set up front and
`aggregate_features` combines them by a rule you chose in advance.


Two checks it runs on the way through
-------------------------------------

**Was the input log transformed.** limma's variance model assumes log-scale
input. On raw intensities the effect size stops being a fold change and starts
scaling with mean expression, which makes any `min_effect` threshold in the
decision rule mean something different for every gene. `scale_check` looks at
the AveExpr range and at whether |logFC| tracks AveExpr, and says what it sees.

**Was the p-value column already corrected.** `adj.P.Val` from limma is
already Benjamini-Hochberg across the whole table. Passing it as `p_value` and
correcting again is double correction; passing raw `P.Value` as `fdr` is no
correction. Both are silent. `load_limma` keeps the two separate and labelled.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path
from typing import Any, Optional, Sequence

__all__ = [
    "aggregate_features",
    "load_limma",
    "scale_check",
    "select_feature",
    "to_analysis_result",
    "to_analysis_unit",
]


# limma's column names, plus the variants that turn up in the wild.
_COLUMNS = {
    "feature": ("probe_id", "GeneID", "ID", "PROBEID", "probe", "gene",
                "Gene.symbol", "gene_id", "feature", "feature_id", "row_id"),
    "effect": ("logFC", "log2FoldChange", "coef", "estimate"),
    "mean": ("AveExpr", "baseMean", "aveexpr"),
    "t": ("t", "stat", "statistic"),
    "p_value": ("P.Value", "PValue", "pvalue", "p_value"),
    "fdr": ("adj.P.Val", "padj", "FDR", "qvalue", "adj_p_val"),
}


def _pick(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lower = {f.lower(): f for f in fieldnames}
    for c in candidates:
        if c in fieldnames:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def _as_float(v: Any) -> Optional[float]:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if f == f else None  # drop NaN


def load_limma(path: str | Path) -> list[dict[str, Any]]:
    """Read a limma topTable TSV into normalized rows.

    Each row carries `feature`, `effect`, `p_value`, `fdr`, `mean_expression`
    and `t`, with the original fields kept under `raw`. `p_value` and `fdr`
    stay separate and keep their meaning: `fdr` is limma's already-corrected
    column, `p_value` is the uncorrected one.
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as fh:
        sniff = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sniff.splitlines()[0] else ","
        reader = csv.DictReader(fh, delimiter=delim)
        fields = reader.fieldnames or []
        cols = {k: _pick(fields, v) for k, v in _COLUMNS.items()}
        if cols["effect"] is None:
            raise ValueError(
                f"no effect column in {path.name}; looked for {_COLUMNS['effect']}, "
                f"found {fields}"
            )
        # No recognised identifier column. limma writes row names into an
        # unnamed first column, so that case is fine. Anything else is not:
        # inventing row0, row1, row2 would silently discard every feature id
        # and turn `select_feature` into a lookup that can never succeed.
        if cols["feature"] is None:
            first = fields[0] if fields else None
            if first is not None and first.strip() == "":
                cols["feature"] = first
            else:
                raise ValueError(
                    f"no feature identifier column in {path.name}. Looked for "
                    f"{_COLUMNS['feature']} and for an unnamed first column; the "
                    f"columns present are {fields}. Add the right name to "
                    "_COLUMNS['feature'] rather than letting rows be numbered."
                )

        rows = []
        for i, r in enumerate(reader):
            feature = (r.get(cols["feature"]) or "").strip()
            if not feature:
                raise ValueError(
                    f"{path.name} row {i + 2}: empty value in the feature column "
                    f"{cols['feature']!r}. Refusing to number it instead."
                )
            rows.append({
                "feature": feature,
                "effect": _as_float(r.get(cols["effect"])),
                "mean_expression": _as_float(r.get(cols["mean"])) if cols["mean"] else None,
                "t": _as_float(r.get(cols["t"])) if cols["t"] else None,
                "p_value": _as_float(r.get(cols["p_value"])) if cols["p_value"] else None,
                "fdr": _as_float(r.get(cols["fdr"])) if cols["fdr"] else None,
                "raw": r,
            })
    if not rows:
        raise ValueError(f"{path.name} has no data rows")
    return rows


# ---------------------------------------------------------------------------
# Declare the feature first
# ---------------------------------------------------------------------------


def select_feature(rows: Sequence[dict[str, Any]], feature: str) -> dict[str, Any]:
    """The row for one named feature, or an error naming what was available.

    There is no `best_feature()` in this module and that absence is the point.
    """
    want = feature.strip().lower()
    hits = [r for r in rows if r["feature"].strip().lower() == want]
    if not hits:
        sample = ", ".join(r["feature"] for r in rows[:3])
        raise KeyError(
            f"{feature!r} is not in this table ({len(rows)} rows, e.g. {sample}). "
            "Declare a feature that exists rather than searching for one that "
            "gives the answer you want."
        )
    if len(hits) > 1:
        raise KeyError(
            f"{feature!r} appears {len(hits)} times. Use aggregate_features() "
            "with an explicit rule instead of taking the first."
        )
    return hits[0]


def aggregate_features(
    rows: Sequence[dict[str, Any]],
    features: Sequence[str],
    rule: str = "mean_effect",
) -> dict[str, Any]:
    """Combine several probes for one gene, by a rule chosen in advance.

    `mean_effect`   average effect, most conservative p among them
    `median_effect` median effect, most conservative p among them

    Neither rule looks at which probe is most significant, because that is the
    choice this module exists to prevent. "Most conservative p" means the
    largest, so combining probes can never improve significance.
    """
    if rule not in ("mean_effect", "median_effect"):
        raise ValueError(f"unknown rule {rule!r}")
    picked = [select_feature(rows, f) for f in features]
    effects = [p["effect"] for p in picked if p["effect"] is not None]
    if not effects:
        raise ValueError("none of the named features carry an effect")

    combine = statistics.fmean if rule == "mean_effect" else statistics.median
    ps = [p["p_value"] for p in picked if p["p_value"] is not None]
    fs = [p["fdr"] for p in picked if p["fdr"] is not None]
    means = [p["mean_expression"] for p in picked if p["mean_expression"] is not None]

    return {
        "feature": "+".join(features),
        "effect": combine(effects),
        "p_value": max(ps) if ps else None,
        "fdr": max(fs) if fs else None,
        "mean_expression": statistics.fmean(means) if means else None,
        "t": None,
        "n_features": len(picked),
        "aggregation": rule,
        "raw": {"features": list(features)},
    }


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def scale_check(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Does this table look like it came from log-transformed input.

    Two signals, neither conclusive alone:

      * log2 expression usually sits inside roughly 2 to 20. Raw microarray
        intensities run to tens of thousands.
      * on raw intensities the effect size scales with mean expression, so
        |effect| and mean correlate strongly. On log data that coupling is
        largely removed, which is the reason for logging in the first place.

    Returns what it measured and a plain verdict. It never edits the data.
    """
    means = [r["mean_expression"] for r in rows if r["mean_expression"] is not None]
    effs = [abs(r["effect"]) for r in rows if r["effect"] is not None]
    pairs = [(r["mean_expression"], abs(r["effect"])) for r in rows
             if r["mean_expression"] is not None and r["effect"] is not None]

    out: dict[str, Any] = {
        "n_rows": len(rows),
        "mean_expression_range": (min(means), max(means)) if means else None,
        "median_abs_effect": statistics.median(effs) if effs else None,
        "max_abs_effect": max(effs) if effs else None,
        "effect_mean_correlation": None,
        "looks_log_transformed": None,
        "note": "",
    }

    if len(pairs) > 10:
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        mx, my = statistics.fmean(xs), statistics.fmean(ys)
        num = sum((x - mx) * (y - my) for x, y in pairs)
        den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
        out["effect_mean_correlation"] = num / den if den else None

    if not means:
        out["note"] = "no mean-expression column, cannot check"
        return out

    big_range = max(means) > 100
    coupled = (out["effect_mean_correlation"] or 0) > 0.4
    out["looks_log_transformed"] = not (big_range or coupled)

    if out["looks_log_transformed"]:
        out["note"] = "Effect sizes look like they are on a log scale."
    else:
        bits = []
        if big_range:
            bits.append(f"mean expression reaches {max(means):,.0f}, far above a log2 range")
        if coupled:
            bits.append(
                f"|effect| tracks mean expression at r={out['effect_mean_correlation']:+.2f}, "
                "which is what raw intensities do"
            )
        out["note"] = (
            "This table does not look log transformed: " + "; and ".join(bits) + ". "
            "limma's variance model assumes log-scale input, and on this scale "
            "the effect is a raw difference rather than a fold change, so a "
            "single min_effect threshold means something different for every gene."
        )
    return out


# ---------------------------------------------------------------------------
# Handing off to synthesis
# ---------------------------------------------------------------------------


def to_analysis_unit(
    row: dict[str, Any],
    analysis_unit_id: str,
    gse_id: str,
    cohort: Optional[str] = None,
    platform: Optional[str] = None,
    independence_group: Optional[str] = None,
    n: Optional[int] = None,
    use_fdr: bool = True,
) -> dict[str, Any]:
    """Shape one selected limma row for `synthesis.synthesize`.

    `use_fdr=True` passes limma's already-corrected column, which is normally
    right. Set it False only when you will correct downstream yourself, and
    then do not correct twice.
    """
    unit = {
        "analysis_unit_id": analysis_unit_id,
        "gse_id": gse_id,
        "cohort": cohort,
        "platform": platform,
        "effect": row["effect"],
        "feature": row["feature"],
    }
    if independence_group:
        unit["independence_group"] = independence_group
    if n is not None:
        unit["n"] = n
    if use_fdr and row.get("fdr") is not None:
        unit["fdr"] = row["fdr"]
    elif row.get("p_value") is not None:
        unit["p_value"] = row["p_value"]
    else:
        raise ValueError(f"{row['feature']} carries neither an fdr nor a p-value")
    return unit


def to_analysis_result(
    row: dict[str, Any],
    analysis_unit_id: str,
    gse_id: str,
    n: int,
    outcome: dict[str, Any],
    test_group: str,
    reference_group: str,
    test_is_higher_outcome: bool,
    predictor_timepoint: Optional[str] = None,
    predictor_transform: Optional[str] = None,
    cohort: Optional[str] = None,
    platform: Optional[str] = None,
    independence_group: Optional[str] = None,
    limma_version: Optional[str] = None,
    scale: Optional[dict[str, Any]] = None,
) -> Any:
    """Emit the canonical `AnalysisResult` for the group-contrast branch.

    `to_analysis_unit` above is the older, flatter handoff. This is the one to
    use when a run has to sit in the same evidence table as a regression,
    because it records what the number is rather than only what it equals.

    `test_is_higher_outcome` has no default. A limma factors file can name its
    groups in either order, and a wrong answer here flips the conclusion while
    leaving every number identical.

    Pass `scale` from `scale_check()` and it is stored in provenance, so a
    table analysed off the log scale carries that fact with it.
    """
    # Works both as a script run from this directory and as
    # `evidence_rules.limma_adapter` imported as part of a package.
    try:
        from .analysis_result import from_group_contrast
    except ImportError:
        from analysis_result import from_group_contrast

    software: dict[str, Any] = {"tool": "limma"}
    if limma_version:
        software["version"] = limma_version

    provenance: dict[str, Any] = {}
    if scale is not None:
        provenance["scale_check"] = {
            "looks_log_transformed": scale.get("looks_log_transformed"),
            "mean_expression_range": scale.get("mean_expression_range"),
            "effect_mean_correlation": scale.get("effect_mean_correlation"),
            "note": scale.get("note"),
        }
    if row.get("aggregation"):
        provenance["probe_aggregation"] = {
            "rule": row["aggregation"],
            "n_features": row.get("n_features"),
        }

    predictor = {"name": row["feature"], "feature_id": row["feature"]}
    if predictor_timepoint:
        predictor["timepoint"] = predictor_timepoint
    if predictor_transform:
        predictor["transform"] = predictor_transform

    return from_group_contrast(
        analysis_unit_id=analysis_unit_id,
        gse_id=gse_id,
        effect=row["effect"],
        n=n,
        predictor=predictor,
        outcome=outcome,
        test_group=test_group,
        reference_group=reference_group,
        test_is_higher_outcome=test_is_higher_outcome,
        effect_type="logFC",
        p_value=row.get("p_value"),
        fdr=row.get("fdr"),
        cohort=cohort,
        platform=platform,
        independence_group=independence_group,
        software=software,
        provenance=provenance,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("table", help="limma topTable TSV or CSV")
    ap.add_argument("--feature", help="the feature declared in advance")
    ap.add_argument("--features", nargs="+", help="several probes for one gene")
    ap.add_argument("--rule", default="mean_effect",
                    choices=("mean_effect", "median_effect"))
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args(argv)

    rows = load_limma(args.table)
    print(f"{len(rows):,} features from {Path(args.table).name}")

    chk = scale_check(rows)
    # Every one of these is None on some real table: no mean column, too few
    # rows to correlate, or a constant column. Formatting None with :.3f
    # raises, so each is printed only when it exists.
    rng = chk["mean_expression_range"]
    print(f"\nmean expression   {rng[0]:,.1f} to {rng[1]:,.1f}" if rng
          else "\nmean expression   not reported")
    med = chk["median_abs_effect"]
    print(f"median |effect|   {med:,.3f}" if med is not None
          else "median |effect|   not reported")
    corr = chk["effect_mean_correlation"]
    print(f"|effect| vs mean  r = {corr:+.3f}" if corr is not None
          else "|effect| vs mean  too few rows to say")

    verdict = chk["looks_log_transformed"]
    label = "OK" if verdict else ("CHECK" if verdict is False else "UNKNOWN")
    print(f"\n{label}: {chk['note']}")

    # `r["fdr"] or 1` would turn an FDR of exactly 0.0 into 1 and drop the
    # most significant features from the count.
    sig = [r for r in rows if r["fdr"] is not None and r["fdr"] < args.alpha]
    print(f"\n{len(sig)} feature(s) at fdr < {args.alpha}")

    if args.features:
        row = aggregate_features(rows, args.features, args.rule)
    elif args.feature:
        row = select_feature(rows, args.feature)
    else:
        print("\nNo feature declared, so nothing was selected. That is the "
              "default on purpose: name the feature before you look.")
        return 0

    print(f"\ndeclared feature: {row['feature']}")
    print(f"  effect  {row['effect']:+.4g}")
    print(f"  p       {row['p_value']:.4g}" if row.get("p_value") is not None else "  p       n/a")
    print(f"  fdr     {row['fdr']:.4g}" if row.get("fdr") is not None else "  fdr     n/a")
    return 0


def _self_test() -> None:
    import tempfile

    content = (
        "probe_id\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB\n"
        "P1\t2.5\t8.0\t4.0\t0.001\t0.01\t3.0\n"
        "P2\t-1.0\t7.5\t-2.0\t0.04\t0.10\t1.0\n"
        "P3\t0.2\t9.0\t0.4\t0.70\t0.85\t-4.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(content)
        p = fh.name

    rows = load_limma(p)
    assert len(rows) == 3
    assert rows[0]["feature"] == "P1"
    assert rows[0]["effect"] == 2.5
    assert rows[0]["fdr"] == 0.01
    assert rows[0]["p_value"] == 0.001

    r = select_feature(rows, "P2")
    assert r["effect"] == -1.0
    assert select_feature(rows, "p2")["effect"] == -1.0   # case insensitive

    try:
        select_feature(rows, "NOPE")
        raise AssertionError("should refuse an absent feature")
    except KeyError as e:
        assert "not in this table" in str(e)

    # No helper exists for picking the winner.
    assert not any("best" in n.lower() or "top" in n.lower() for n in __all__)

    agg = aggregate_features(rows, ["P1", "P2"])
    assert agg["effect"] == 0.75                 # mean of 2.5 and -1.0
    assert agg["fdr"] == 0.10                    # most conservative, never the best
    assert agg["p_value"] == 0.04
    assert agg["n_features"] == 2
    assert aggregate_features(rows, ["P1", "P2", "P3"], "median_effect")["effect"] == 0.2

    # Log-scale data passes the check.
    chk = scale_check(rows)
    assert chk["looks_log_transformed"] is True, chk

    # Raw-intensity data fails it.
    raw = [{"feature": f"P{i}", "effect": float(i * 10), "mean_expression": float(i * 100),
            "p_value": 0.5, "fdr": 0.9, "t": None, "raw": {}} for i in range(1, 40)]
    chk2 = scale_check(raw)
    assert chk2["looks_log_transformed"] is False
    assert "does not look log transformed" in chk2["note"]
    assert chk2["effect_mean_correlation"] > 0.9

    unit = to_analysis_unit(rows[0], "GSE1_trial1", "GSE1", cohort="Trial 1",
                            platform="GPL1", independence_group="g1", n=15)
    assert unit["fdr"] == 0.01 and "p_value" not in unit
    assert unit["cohort"] == "Trial 1"
    unit2 = to_analysis_unit(rows[0], "u", "GSE1", use_fdr=False)
    assert unit2["p_value"] == 0.001 and "fdr" not in unit2

    try:
        to_analysis_unit({"feature": "X", "effect": 1.0}, "u", "GSE1")
        raise AssertionError("should refuse a row with no p-value at all")
    except ValueError:
        pass

    # The canonical contract, with the scale check carried in provenance.
    from analysis_result import oriented_effect, to_synthesis_unit, validate_result

    res = to_analysis_result(
        rows[0], "GSE1_all", "GSE1", n=36,
        outcome={"name": "day-84 nAb", "timepoint": "day 84",
                 "scale": "median split"},
        test_group="high", reference_group="low", test_is_higher_outcome=True,
        predictor_timepoint="baseline", predictor_transform="log2 quantile",
        platform="GPL10558", limma_version="3.58.1+galaxy0",
        scale=scale_check(rows),
    )
    assert validate_result(res) == [], validate_result(res)
    assert res.effect_type == "logFC"
    assert res.outcome["dichotomised"] is True
    assert res.provenance["scale_check"]["looks_log_transformed"] is True
    assert res.software["version"] == "3.58.1+galaxy0"
    assert oriented_effect(res) == 2.5
    assert to_synthesis_unit(res)["effect"] == 2.5

    # Reversed factors file, same numbers, opposite meaning.
    flipped = to_analysis_result(
        rows[0], "GSE1_all", "GSE1", n=36,
        outcome={"name": "day-84 nAb"},
        test_group="low", reference_group="high", test_is_higher_outcome=False,
    )
    assert oriented_effect(flipped) == -2.5

    # An aggregated probe set records how it was combined.
    agg_res = to_analysis_result(
        agg, "GSE1_all", "GSE1", n=36, outcome={"name": "out"},
        test_group="high", reference_group="low", test_is_higher_outcome=True,
    )
    assert agg_res.provenance["probe_aggregation"]["rule"] == "mean_effect"
    assert agg_res.provenance["probe_aggregation"]["n_features"] == 2

    # Galaxy's limma names the id column GeneID, not probe_id.
    gx = ("GeneID\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB\n"
          "ILMN_1\t3.6\t5.7\t49.0\t1.9e-44\t5.1e-40\t54.9\n")
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write(gx)
        gp = fh.name
    assert load_limma(gp)[0]["feature"] == "ILMN_1"

    # An unnamed first column is limma's row-name convention and is accepted.
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write("\tlogFC\tP.Value\nILMN_9\t1.0\t0.01\n")
        rp = fh.name
    assert load_limma(rp)[0]["feature"] == "ILMN_9"

    # An FDR of exactly 0.0 is the most significant value there is, and
    # `r["fdr"] or 1` would have turned it into 1 and dropped it.
    zero = [{"feature": "P1", "effect": 3.0, "fdr": 0.0, "p_value": 0.0,
             "mean_expression": 8.0, "t": None, "raw": {}}]
    assert [r for r in zero if r["fdr"] is not None and r["fdr"] < 0.05]
    assert not [r for r in zero if (r["fdr"] or 1) < 0.05]     # the old bug

    # scale_check returns None for what it cannot compute; nothing may assume
    # a number is there. Two rows is too few to correlate.
    thin = scale_check(zero)
    assert thin["effect_mean_correlation"] is None
    assert scale_check([{"feature": "x", "effect": 1.0, "mean_expression": None,
                         "t": None, "p_value": 0.1, "fdr": 0.1,
                         "raw": {}}])["looks_log_transformed"] is None

    # An unrecognised named id column is refused, not silently numbered.
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write("weird_name\tlogFC\tP.Value\nILMN_9\t1.0\t0.01\n")
        wp = fh.name
    try:
        load_limma(wp)
        raise AssertionError("should refuse a table with no known feature column")
    except ValueError as e:
        assert "no feature identifier column" in str(e)

    print("self-tests: 17 groups of assertions passed")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        _self_test()
    else:
        raise SystemExit(main())
