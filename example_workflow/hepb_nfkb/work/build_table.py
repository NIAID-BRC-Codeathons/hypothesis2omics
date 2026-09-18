#!/usr/bin/env python3
"""Join GSE65834 baseline NF-kB expression to anti-HBs outcome.

Outcome precedence: ImmPort elisa_result VALUE_PREFERRED (day 7 post-vax reading
is the study's post-vaccination titer) with GEO characteristics as cross-check.
"""
import csv, gzip, json, hashlib, re, sys
from pathlib import Path

D = Path("/nfs/lambda_stor_01/homes/avasan/.loom/mcp/hypothesis2omics/mcp/data")
W = Path("/nfs/lambda_stor_01/data/avasan/omics_agents/hepb_nfkb/work")
EXPR = D/"geo_cache/parsed/SDY1328/EXP16735/GSE65834_GPL15048/expression.tsv.gz"
META = D/"immport_cache/parsed/sample_manifest.tsv"

PROBES = {
 "NFKB1":["merck-NM_003998_at"],
 "RELA":["merck2-AA102377_at","merck2-AA102377_x_at","merck2-BC069248_at","merck-NM_021975_s_at"],
 "NFKB2":["merck-NM_002502_at"],
 "REL":["merck-BP294455_a_at","merck-NM_002908_at"],
 "RELB":["merck-NM_006509_at"],
}
allp = {p:g for g,ps in PROBES.items() for p in ps}

# ---- GSM -> subject (baseline, day 0, this unit only)
gsm2sub={}
for r in csv.DictReader(open(META),delimiter="\t"):
    if r["study_accession"]=="SDY1328" and r["study_time_collected"]=="0":
        if r["repository_accession"].startswith("GSM"):
            gsm2sub[r["repository_accession"]]=r["subject_accession"]
print(f"baseline GSM->subject: {len(gsm2sub)}")

# ---- outcome: ImmPort ELISA anti-HBs
import zipfile
z=zipfile.ZipFile(D/"immport_cache/SDY1328/SDY1328-DR58_Tab.zip")
name=[n for n in z.namelist() if n.endswith("elisa_result.txt")][0]
el=list(csv.DictReader(z.open(name).read().decode("utf8","replace").splitlines(),delimiter="\t"))
post={}; pre={}
for r in el:
    if r["ANALYTE_REPORTED"]!="Hepatitis B Virus Surface Antibody": continue
    v=r["VALUE_PREFERRED"].strip()
    if v in ("","NA"): continue
    (post if r["STUDY_TIME_COLLECTED"]=="7" else pre)[r["SUBJECT_ACCESSION"]]=float(v)
print(f"ImmPort anti-HBs: pre(d0)={len(pre)}  post(d7)={len(post)}")

# ---- expression rows for our probes
with gzip.open(EXPR,"rt") as fh:
    rd=csv.reader(fh,delimiter="\t"); hdr=next(rd)
    cols=hdr[1:]
    vals={}
    for row in rd:
        if row[0] in allp:
            vals[row[0]]={c:float(x) for c,x in zip(cols,row[1:]) if x not in ("","NA")}
print(f"probe rows found: {len(vals)}/{len(allp)}; arrays: {len(cols)}")

# ---- assemble one row per subject
recs=[]
for gsm,sub in sorted(gsm2sub.items()):
    if gsm not in cols: continue
    if sub not in post: continue
    rec={"gsm":gsm,"subject":sub,"antiHBs_post":post[sub],"antiHBs_pre":pre.get(sub,"")}
    ok=True
    for g,ps in PROBES.items():
        xs=[vals[p][gsm] for p in ps if p in vals and gsm in vals[p]]
        if not xs: ok=False; break
        rec[g]=sum(xs)/len(xs)      # mean across probes of a gene
    if ok: recs.append(rec)
print(f"analysis rows (baseline expr + post titer): {len(recs)}")

out=W/"analysis_table.tsv"
flds=["gsm","subject","antiHBs_pre","antiHBs_post"]+list(PROBES)
with open(out,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=flds,delimiter="\t"); w.writeheader()
    for r in recs: w.writerow(r)
print("wrote",out, hashlib.sha256(out.read_bytes()).hexdigest()[:16])
