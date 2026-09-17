#!/usr/bin/env python3
"""
Run the YF-17D probe-vs-titer Spearman correlation on Galaxy, using the
Galaxy tool `Feature-wise Correlation Tests` instead of local scipy, and
verify the result against the local reference tables.

This is the scripted form of the cross-check recorded in notebook.md
(Session 5): Galaxy job bbd44e69cb8906b536b291633e51df67 in history
bbd44e69cb8906b54e30cdcb10bc7945 agreed with local scipy to ~5e-07.

WHAT RUNS WHERE
---------------
  Galaxy  : the per-feature Spearman + BH adjustment (the arithmetic)
  local   : subject-level collapse, d7-d0 fold change, within-arm rank
            de-confounding, the log2=7.0 floor filter

  The local half is NOT reproducible in Galaxy -- no tool implements it --
  and every consequential design decision lives there. This script drives
  the Galaxy half and checks it; it does not make the prep reproducible.

ENCODING
--------
The tool correlates feature f in matrix A against feature f in matrix B,
pairing features BY NAME when headers are present. To express "every probe
against ONE titer vector":

    matrix A = subjects x probes   (expression feature)
    matrix B = subjects x probes   (every column = the same titer vector)

Feature-wise correlation then reduces exactly to Spearman(probe_i, titer).
Both matrices are written by export_for_galaxy.py, which must be run first.

CREDENTIALS
-----------
Read from the environment only -- never inline a key in this file:
    export GALAXY_URL=https://usegalaxy.org
    export GALAXY_API_KEY=<your key>     # already set in this project

USAGE
-----
    python3 yf17d/export_for_galaxy.py --list
    python3 yf17d/export_for_galaxy.py "A__FC_d7-d0__rank"
    python3 yf17d/run_galaxy_correlation.py "A__FC_d7-d0__rank"

    # reuse an existing history instead of creating one:
    python3 yf17d/run_galaxy_correlation.py "A__FC_d7-d0__rank" \
        --history-id bbd44e69cb8906b54e30cdcb10bc7945
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time

from bioblend.galaxy import GalaxyInstance

TOOL_ID = ("toolshed.g2.bx.psu.edu/repos/goeckslab/featurewise_correlation/"
           "featurewise_correlation/0.1.0+galaxy3")

HERE = os.path.dirname(os.path.abspath(__file__))
GIN = os.path.join(HERE, "galaxy_in")
GOUT = os.path.join(HERE, "galaxy_out")
LOCAL_OUT = os.path.join(HERE, "out")

# Tool parameters. These mirror the verified run recorded in notebook.md;
# changing any of them invalidates the comparison against the local tables.
TOOL_PARAMS = {
    "delimiter_a": "tab",
    "delimiter_b": "tab",
    "has_header_a": "true",
    "has_header_b": "true",
    "id_column_a": 1,
    "id_column_b": 1,
    "feature_start_a": 2,
    "feature_start_b": 2,
    "transform_a": "none",
    "transform_b": "none",
    "method": "spearman",
    "alternative": "two-sided",
    "alpha": 0.05,
    "min_pairs": 3,
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def connect():
    url = os.environ.get("GALAXY_URL", "https://usegalaxy.org")
    key = os.environ.get("GALAXY_API_KEY")
    if not key:
        sys.exit("GALAXY_API_KEY is not set. Export it; do not paste it into "
                 "this file or into chat.")
    gi = GalaxyInstance(url=url, key=key)
    me = gi.users.get_current_user()
    print(f"[galaxy] {url} as {me['username']} "
          f"({me['nice_total_disk_usage']} used)")
    return gi


def upload(gi, history_id, path, name):
    """Upload a local TSV and block until Galaxy finishes ingesting it."""
    r = gi.tools.upload_file(path, history_id, file_name=name,
                             file_type="tabular")
    ds_id = r["outputs"][0]["id"]
    while True:
        ds = gi.datasets.show_dataset(ds_id)
        if ds["state"] in ("ok", "error"):
            break
        time.sleep(3)
    if ds["state"] != "ok":
        sys.exit(f"upload failed for {name}: state={ds['state']}")
    print(f"[upload] hid {ds['hid']:>2}  {name}  ({ds['file_size']} bytes)")
    return ds_id


def wait_for_job(gi, job_id, poll=5):
    """Poll a tool job to a terminal state. Returns the final job record."""
    while True:
        job = gi.jobs.show_job(job_id)
        state = job["state"]
        if state in ("ok", "error", "deleted", "paused"):
            return job
        print(f"[job] {state} ...")
        time.sleep(poll)


def compare_to_local(galaxy_tsv, local_tsv):
    """Probe-by-probe diff of the Galaxy table against the local scipy table.

    Returns a dict of check results. Deliberately compares rho, raw p AND the
    BH-adjusted q -- agreeing on rho alone would not establish that the two
    implementations adjust the same way.
    """
    gx = {}
    with open(galaxy_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            gx[r["feature"]] = (int(r["n_pairs"]), float(r["statistic"]),
                                float(r["p_value"]), float(r["p_adjust_bh"]),
                                r["significant"])
    lo = {}
    with open(local_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lo[r["probe_id"]] = (float(r["rho"]), float(r["p"]),
                                 float(r["q_BH"]))

    if set(gx) != set(lo):
        return {"feature_sets_identical": False,
                "n_galaxy": len(gx), "n_local": len(lo)}

    keys = sorted(lo)
    d_rho = max(abs(gx[k][1] - lo[k][0]) for k in keys)
    d_p = max(abs(gx[k][2] - lo[k][1]) for k in keys)
    d_q = max(abs(gx[k][3] - lo[k][2]) for k in keys)

    g_by_p = sorted(keys, key=lambda k: gx[k][2])
    l_by_p = sorted(keys, key=lambda k: lo[k][1])
    top = min(1000, len(keys))

    return {
        "feature_sets_identical": True,
        "n_features": len(keys),
        "max_abs_rho_diff": d_rho,
        "max_abs_p_diff": d_p,
        "max_abs_q_diff": d_q,
        "galaxy_p_lt_01": sum(1 for k in keys if gx[k][2] < 0.01),
        "local_p_lt_01": sum(1 for k in keys if lo[k][1] < 0.01),
        "galaxy_p_lt_001": sum(1 for k in keys if gx[k][2] < 0.001),
        "local_p_lt_001": sum(1 for k in keys if lo[k][1] < 0.001),
        "galaxy_significant_true": sum(1 for k in keys if gx[k][4] == "true"),
        "galaxy_q_lt_05": sum(1 for k in keys if gx[k][3] < 0.05),
        "local_q_lt_05": sum(1 for k in keys if lo[k][2] < 0.05),
        f"top{top}_overlap": len(set(g_by_p[:top]) & set(l_by_p[:top])),
        f"top{top}_of": top,
        "top_hit_galaxy": {"probe": g_by_p[0], "rho": gx[g_by_p[0]][1],
                           "p": gx[g_by_p[0]][2], "q": gx[g_by_p[0]][3]},
        "top_hit_local": {"probe": l_by_p[0], "rho": lo[l_by_p[0]][0],
                          "p": lo[l_by_p[0]][1], "q": lo[l_by_p[0]][2]},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("analysis_key",
                    help='e.g. "A__FC_d7-d0__rank" (see export_for_galaxy.py --list)')
    ap.add_argument("--history-id", default=None,
                    help="reuse an existing history instead of creating one")
    ap.add_argument("--history-name", default=None)
    args = ap.parse_args()

    key = args.analysis_key
    safe = key.replace("-", "m")
    mat_a = os.path.join(GIN, f"matA_{safe}.tsv")
    mat_b = os.path.join(GIN, f"matB_{safe}.tsv")
    for p in (mat_a, mat_b):
        if not os.path.exists(p):
            sys.exit(f"missing {p}\nRun first:  python3 yf17d/export_for_galaxy.py "
                     f'"{key}"')
    os.makedirs(GOUT, exist_ok=True)

    print(f"[inputs] matA sha256 {sha256(mat_a)}")
    print(f"[inputs] matB sha256 {sha256(mat_b)}")

    gi = connect()

    if args.history_id:
        hist_id = args.history_id
        print(f"[history] reusing {hist_id}")
    else:
        name = args.history_name or (
            f"YF-17D probe-vs-nAb-titer Spearman [{key}]")
        hist_id = gi.histories.create_history(name=name)["id"]
        print(f"[history] created {hist_id}  {name}")

    a_id = upload(gi, hist_id, mat_a, f"matA_expression [{key}]")
    b_id = upload(gi, hist_id, mat_b, f"matB_titer_replicated [{key}]")

    inputs = dict(TOOL_PARAMS)
    inputs["matrix_a"] = {"src": "hda", "id": a_id}
    inputs["matrix_b"] = {"src": "hda", "id": b_id}

    run = gi.tools.run_tool(hist_id, TOOL_ID, inputs)
    job_id = run["jobs"][0]["id"]
    per_feature_id = next(o["id"] for o in run["outputs"]
                          if o["output_name"] == "per_feature")
    summary_id = next(o["id"] for o in run["outputs"]
                      if o["output_name"] == "summary")
    print(f"[job] submitted {job_id}")

    job = wait_for_job(gi, job_id)
    print(f"[job] finished state={job['state']} exit_code={job.get('exit_code')}")
    if job["state"] != "ok":
        sys.exit(f"job did not succeed; inspect {job_id} in Galaxy")

    per_feature = os.path.join(GOUT, f"perfeature_{safe}.tsv")
    summary = os.path.join(GOUT, f"summary_{safe}.tsv")
    with open(per_feature, "wb") as f:
        f.write(gi.datasets.download_dataset(per_feature_id))
    with open(summary, "wb") as f:
        f.write(gi.datasets.download_dataset(summary_id))
    print(f"[download] {per_feature}")
    print(f"[download] {summary}")

    print("\n--- Galaxy summary ---")
    with open(summary) as f:
        print(f.read().strip())

    local_tsv = os.path.join(LOCAL_OUT, f"corr_{key}.tsv")
    result = {
        "analysis_key": key,
        "galaxy_url": os.environ.get("GALAXY_URL", "https://usegalaxy.org"),
        "tool_id": TOOL_ID,
        "tool_params": TOOL_PARAMS,
        "history_id": hist_id,
        "job_id": job_id,
        "job_state": job["state"],
        "inputs": {"matrix_a_sha256": sha256(mat_a),
                   "matrix_b_sha256": sha256(mat_b)},
        "outputs": {"per_feature": per_feature,
                    "per_feature_sha256": sha256(per_feature),
                    "summary": summary},
    }

    if os.path.exists(local_tsv):
        cmp = compare_to_local(per_feature, local_tsv)
        result["comparison_vs_local"] = cmp
        result["local_reference"] = local_tsv
        print("\n--- comparison vs local scipy ---")
        print(json.dumps(cmp, indent=2))
        agree = (cmp.get("feature_sets_identical")
                 and cmp.get("max_abs_rho_diff", 1) < 1e-5
                 and cmp.get("max_abs_p_diff", 1) < 1e-5
                 and cmp.get("max_abs_q_diff", 1) < 1e-5
                 and cmp.get("galaxy_p_lt_01") == cmp.get("local_p_lt_01")
                 and cmp.get("galaxy_q_lt_05") == cmp.get("local_q_lt_05"))
        result["agrees_with_local"] = bool(agree)
        print(f"\nAGREEMENT: {'PASS' if agree else 'FAIL'}")
        if cmp.get("galaxy_q_lt_05") == 0:
            print("NOTE: zero probes at FDR<0.05 -- agreement confirms the "
                  "arithmetic, not a finding. A reproduced null is still null.")
    else:
        result["comparison_vs_local"] = None
        print(f"\n[warn] no local reference at {local_tsv}; skipped comparison")

    rp = os.path.join(GOUT, f"galaxy_run_{safe}.json")
    with open(rp, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
    print(f"[provenance] {rp}")


if __name__ == "__main__":
    main()
