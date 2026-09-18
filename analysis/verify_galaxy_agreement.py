#!/usr/bin/env python3
"""Re-check the Galaxy-versus-local agreement claim from the committed tables.

Reads each Galaxy per-feature table and its local reference table, and
recomputes the comparison rather than trusting the comparison_vs_local block
that run_galaxy_correlation.py wrote into the provenance JSON.

Also verifies each per-feature table against the per_feature_sha256 recorded in
its provenance JSON, so a substituted or truncated file is caught.

Usage:
    python analysis/verify_galaxy_agreement.py

Exits 0 if every check passes, 1 otherwise. Standard library only.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# analysis_key -> (provenance json, galaxy per-feature table, local reference table)
PAIRS = {
    "A__FC_d7-d0__rank": (
        "galaxy_out/galaxy_run_A__FC_d7md0__rank.json",
        "galaxy_out/perfeature_A__FC_d7md0__rank.tsv",
        "out/corr_A__FC_d7-d0__rank.tsv",
    ),
    "A_d60arm__FC_d7-d0__log2titer": (
        "galaxy_out/galaxy_run_A_d60arm__FC_d7md0__log2titer.json",
        "galaxy_out/perfeature_A_d60arm__FC_d7md0__log2titer.tsv",
        "out/corr_A_d60arm__FC_d7-d0__log2titer.tsv",
    ),
}

# The agreement asserted in analysis/galaxy_out/README.md and in the PR.
TOLERANCE = 1e-06


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path, key_col: str) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return {row[key_col]: row for row in csv.DictReader(fh, delimiter="\t")}


def check(analysis_key: str, prov_rel: str, galaxy_rel: str, local_rel: str) -> bool:
    prov_path, galaxy_path, local_path = (ROOT / p for p in (prov_rel, galaxy_rel, local_rel))

    missing = [p for p in (prov_path, galaxy_path, local_path) if not p.exists()]
    if missing:
        for p in missing:
            print(f"  MISSING {p.relative_to(ROOT)}")
        return False

    prov = json.loads(prov_path.read_text(encoding="utf-8"))

    recorded = prov["outputs"]["per_feature_sha256"]
    actual = sha256(galaxy_path)
    if recorded != actual:
        print(f"  FAIL sha256 mismatch on {galaxy_path.name}")
        print(f"       recorded {recorded}")
        print(f"       actual   {actual}")
        return False
    print(f"  ok   sha256 matches provenance ({recorded[:16]}...)")

    galaxy = read_tsv(galaxy_path, "feature")
    local = read_tsv(local_path, "probe_id")

    if set(galaxy) != set(local):
        only_g, only_l = len(set(galaxy) - set(local)), len(set(local) - set(galaxy))
        print(f"  FAIL feature sets differ: {only_g} only in Galaxy, {only_l} only in local")
        return False
    print(f"  ok   feature sets identical ({len(galaxy)} features)")

    d_rho = d_p = d_q = 0.0
    for probe, grow in galaxy.items():
        lrow = local[probe]
        d_rho = max(d_rho, abs(float(grow["statistic"]) - float(lrow["rho"])))
        d_p = max(d_p, abs(float(grow["p_value"]) - float(lrow["p"])))
        d_q = max(d_q, abs(float(grow["p_adjust_bh"]) - float(lrow["q_BH"])))

    print(f"  ok   max abs rho diff {d_rho:.6g}")
    print(f"  ok   max abs p   diff {d_p:.6g}")
    print(f"  ok   max abs q   diff {d_q:.6g}")

    if max(d_rho, d_p, d_q) > TOLERANCE:
        print(f"  FAIL exceeds stated tolerance of {TOLERANCE:g}")
        return False

    by_g = sorted(galaxy, key=lambda k: (float(galaxy[k]["p_value"]), k))[:1000]
    by_l = sorted(local, key=lambda k: (float(local[k]["p"]), k))[:1000]
    overlap = len(set(by_g) & set(by_l))
    print(f"  {'ok  ' if overlap == 1000 else 'FAIL'} top-1000 overlap {overlap}/1000")
    if overlap != 1000:
        return False

    n_g = sum(1 for r in galaxy.values() if float(r["p_adjust_bh"]) < 0.05)
    n_l = sum(1 for r in local.values() if float(r["q_BH"]) < 0.05)
    print(f"  ok   q<0.05: galaxy {n_g}, local {n_l}")
    if n_g or n_l:
        print("  note a feature now passes FDR; the committed READMEs say zero do")
        return False

    return True


def main() -> int:
    ok = True
    for analysis_key, rels in PAIRS.items():
        print(f"\n{analysis_key}")
        if not check(analysis_key, *rels):
            ok = False

    print("\n" + ("all checks passed" if ok else "CHECKS FAILED"))
    print("Agreement confirms the arithmetic, not a finding. A reproduced null is still null.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
