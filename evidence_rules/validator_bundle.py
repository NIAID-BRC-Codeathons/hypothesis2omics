"""Locate the validator input bundle, whichever directory name it is under.

    python validator_bundle.py            # worked example, self-tests

Melanie's end-to-end pipeline trace found this, and she was right. There are
three directory names in the repository for what is conceptually one
three-table bundle:

    data/validator_input/       what validator_handoff_parse_module.py writes,
                                what validator_input_manifest.json declares,
                                what scientific_validator/README.md calls
                                canonical, and what handoff_adapter.py reads
    data/validator_immport/     sample_manifest.tsv, quantitative_outcome.tsv
    data/validator_geo/         feature_expression.tsv

The split pair is organised by which repository the data came from. The
canonical directory is organised by which pipeline stage produced it. Those are
different axes, and the second one is the right one, because a consumer cares
what stage it is downstream of, not which upstream API the bytes arrived from.

`run_yf17d.py` was reading the split pair. So a freshly generated bundle from
step 10 could not feed it without someone hand-copying files into the other
naming, which is a real reason the pipeline did not join up end to end.

This module resolves either layout, prefers the canonical one, and refuses when
the two disagree.


Why a disagreement is fatal rather than a warning
-------------------------------------------------

While both layouts exist, the same table lives at two paths. Today they are
byte-identical. If someone regenerates step 10 and the split copies go stale,
then reading one while the other differs means the analysis is running on data
nobody chose, and the report would carry provenance for a file that is not the
one a reader would find if they went looking.

So `resolve` raises on divergence by default. Pass `allow_divergent=True` to
proceed deliberately, and the caveat it returns says which file was used and
which was ignored. This is the same discipline as the probe-identity gate in
`run_h2.py`: when the code cannot tell which of two answers is right, it should
stop rather than pick.

Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

__all__ = [
    "CANONICAL_DIR",
    "TABLES",
    "Resolved",
    "BundleDivergence",
    "resolve",
    "resolve_all",
    "bundle_caveats",
    "declared_rows",
]


CANONICAL_DIR = "data/validator_input"

# Where the same table also lives, from before the canonical directory existed.
FALLBACK_DIR = {
    "sample_manifest.tsv": "data/validator_immport",
    "quantitative_outcome.tsv": "data/validator_immport",
    "feature_expression.tsv": "data/validator_geo",
}

TABLES: tuple[str, ...] = (
    "sample_manifest.tsv",
    "quantitative_outcome.tsv",
    "feature_expression.tsv",
)

MANIFEST = "data/validator_input/validator_input_manifest.json"


class BundleDivergence(RuntimeError):
    """The canonical copy and the legacy copy of one table do not match."""


@dataclass(frozen=True)
class Resolved:
    """Which file was used for one table, and what else was lying around."""

    name: str
    path: Path
    directory: str
    is_canonical: bool
    duplicate_at: Optional[Path] = None
    duplicate_differs: bool = False

    def caveat(self) -> Optional[str]:
        if self.duplicate_differs:
            return (
                f"{self.name} exists at two paths that do not match. Used "
                f"{self.directory}/{self.name}; ignored {self.duplicate_at}. "
                "One of them is stale and the analysis below depends on which."
            )
        if not self.is_canonical:
            return (
                f"{self.name} was read from {self.directory}/, not from the "
                f"canonical {CANONICAL_DIR}/. That directory is what "
                "validator_handoff_parse_module.py writes, so this run is on "
                "pre-committed tables rather than a freshly generated bundle."
            )
        return None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(repo: Path, name: str, allow_divergent: bool = False) -> Resolved:
    """Find one table. Canonical location wins; legacy location is a fallback.

    Raises FileNotFoundError when neither exists, and BundleDivergence when
    both exist with different contents.
    """
    if name not in FALLBACK_DIR:
        raise KeyError(
            f"{name!r} is not part of the validator bundle. Known tables: "
            + ", ".join(TABLES)
        )
    repo = Path(repo)
    canonical = repo / CANONICAL_DIR / name
    legacy = repo / FALLBACK_DIR[name] / name

    if canonical.exists() and legacy.exists():
        differs = _sha256(canonical) != _sha256(legacy)
        if differs and not allow_divergent:
            raise BundleDivergence(
                f"{name} differs between {CANONICAL_DIR}/ and "
                f"{FALLBACK_DIR[name]}/. Both are present and they are not the "
                "same file, so which one this analysis should use is a "
                "question the code cannot answer. Regenerate the bundle with "
                "`python data/validator_handoff_parse_module.py --config "
                "data/validator_handoff_config.json` and delete the stale "
                "copy, or pass allow_divergent=True to proceed on the "
                "canonical one deliberately."
            )
        return Resolved(name, canonical, CANONICAL_DIR, True, legacy, differs)

    if canonical.exists():
        return Resolved(name, canonical, CANONICAL_DIR, True)
    if legacy.exists():
        return Resolved(name, legacy, FALLBACK_DIR[name], False)

    raise FileNotFoundError(
        f"{name} is in neither {CANONICAL_DIR}/ nor {FALLBACK_DIR[name]}/ "
        f"under {repo}. Generate it with `python "
        "data/validator_handoff_parse_module.py --config "
        "data/validator_handoff_config.json`."
    )


def resolve_all(
    repo: Path,
    names: Sequence[str] = TABLES,
    allow_divergent: bool = False,
) -> dict[str, Resolved]:
    return {n: resolve(repo, n, allow_divergent=allow_divergent) for n in names}


def bundle_caveats(resolved: dict[str, Resolved]) -> list[str]:
    """Anything a report should disclose about where its inputs came from."""
    out: list[str] = []
    for r in resolved.values():
        c = r.caveat()
        if c:
            out.append(c)
    return out


def declared_rows(repo: Path) -> dict[str, int]:
    """Row counts step 10 recorded, for cross-checking what we actually read.

    Returns an empty dict when the manifest is absent, because a missing
    manifest is a normal state for a checkout that has not run step 10.
    """
    path = Path(repo) / MANIFEST
    if not path.exists():
        return {}
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, int] = {}
    for key, art in (doc.get("artifacts") or {}).items():
        rows = art.get("rows")
        p = art.get("path", "")
        if isinstance(rows, int) and p:
            out[Path(p).name] = rows
    return out


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _write(p: Path, text: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def _self_test() -> None:
    import tempfile

    n = 0
    with tempfile.TemporaryDirectory() as td:
        repo = Path(td)

        # -- neither layout present ----------------------------------------
        try:
            resolve(repo, "sample_manifest.tsv")
            raise AssertionError("resolved a table that does not exist")
        except FileNotFoundError as e:
            assert "validator_handoff_parse_module" in str(e), "no recovery hint"
            n += 1

        # -- legacy only, which is the repo as it stands today -------------
        _write(repo / "data/validator_immport/sample_manifest.tsv", "a\n")
        _write(repo / "data/validator_immport/quantitative_outcome.tsv", "b\n")
        _write(repo / "data/validator_geo/feature_expression.tsv", "c\n")
        r = resolve_all(repo)
        assert all(not x.is_canonical for x in r.values())
        assert r["feature_expression.tsv"].directory == "data/validator_geo"
        caveats = bundle_caveats(r)
        assert len(caveats) == 3
        assert all("pre-committed tables" in c for c in caveats)
        n += 1

        # -- canonical present and identical: canonical wins, silently ------
        for name, text in (("sample_manifest.tsv", "a\n"),
                           ("quantitative_outcome.tsv", "b\n"),
                           ("feature_expression.tsv", "c\n")):
            _write(repo / CANONICAL_DIR / name, text)
        r = resolve_all(repo)
        assert all(x.is_canonical for x in r.values())
        assert all(x.duplicate_at is not None for x in r.values())
        assert not any(x.duplicate_differs for x in r.values())
        assert bundle_caveats(r) == [], "identical duplicates should not warn"
        n += 1

        # -- canonical present and DIFFERENT: refuse ------------------------
        _write(repo / CANONICAL_DIR / "feature_expression.tsv", "c CHANGED\n")
        try:
            resolve(repo, "feature_expression.tsv")
            raise AssertionError("resolved a table whose two copies disagree")
        except BundleDivergence as e:
            assert "which one this analysis should use" in str(e)
            assert "validator_handoff_parse_module" in str(e)
            n += 1

        # ...and proceed only when told to, with a caveat naming both
        r = resolve(repo, "feature_expression.tsv", allow_divergent=True)
        assert r.is_canonical and r.duplicate_differs
        c = r.caveat()
        assert c is not None and "do not match" in c and "stale" in c
        n += 1

        # the other two tables still resolve cleanly, so one bad table does
        # not take the whole bundle down when the override is used
        r = resolve_all(repo, allow_divergent=True)
        assert sum(1 for x in r.values() if x.duplicate_differs) == 1
        assert len(bundle_caveats(r)) == 1
        n += 1

        # -- canonical only, the state after a fresh step 10 ----------------
        for name in TABLES:
            (repo / FALLBACK_DIR[name] / name).unlink(missing_ok=True)
        _write(repo / CANONICAL_DIR / "feature_expression.tsv", "c\n")
        r = resolve_all(repo)
        assert all(x.is_canonical and x.duplicate_at is None for x in r.values())
        assert bundle_caveats(r) == []
        n += 1

        # -- the declared row counts ---------------------------------------
        assert declared_rows(repo) == {}, "absent manifest should be empty, not an error"
        _write(repo / MANIFEST, json.dumps({"artifacts": {
            "sample_manifest": {"path": "data/validator_input/sample_manifest.tsv",
                                "rows": 888},
            "feature_expression": {"path": "data/validator_input/feature_expression.tsv",
                                   "rows": 87},
            "quantitative_outcome": {"path": "data/validator_input/quantitative_outcome.tsv",
                                     "rows": 25},
        }}))
        assert declared_rows(repo) == {"sample_manifest.tsv": 888,
                                       "feature_expression.tsv": 87,
                                       "quantitative_outcome.tsv": 25}
        _write(repo / MANIFEST, "not json")
        assert declared_rows(repo) == {}, "unparseable manifest should not raise"
        n += 1

        # -- an unknown table is a mistake, not a fallback ------------------
        try:
            resolve(repo, "something_else.tsv")
            raise AssertionError("accepted a table outside the bundle")
        except KeyError:
            n += 1

    print(f"validator_bundle.py: {n} assertion groups passed")


def _demo(repo: Optional[Path]) -> None:
    if repo is None:
        print("Pass --repo to see how a real checkout resolves.")
        return
    print(f"Resolving the validator bundle under {repo}\n")
    try:
        r = resolve_all(repo)
    except (FileNotFoundError, BundleDivergence) as e:
        print(f"  {type(e).__name__}: {e}")
        return
    for name, res in r.items():
        tag = "canonical" if res.is_canonical else "LEGACY"
        dup = ""
        if res.duplicate_at is not None:
            dup = " (duplicate differs)" if res.duplicate_differs else " (duplicate identical)"
        print(f"  {name:<26} {tag:<10} {res.directory}{dup}")
    declared = declared_rows(repo)
    if declared:
        print("\n  row counts step 10 declared:")
        for k, v in sorted(declared.items()):
            print(f"    {k:<26} {v}")
    cav = bundle_caveats(r)
    print("\n  caveats for the report:" if cav else "\n  no caveats")
    for c in cav:
        print(f"    - {c}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", type=Path, help="checkout to resolve against")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return 0
    _demo(args.repo)
    print()
    _self_test()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
