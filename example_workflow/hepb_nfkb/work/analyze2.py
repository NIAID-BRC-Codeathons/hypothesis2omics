#!/usr/bin/env python3
import csv,json,math
import numpy as np
from scipy import stats
from pathlib import Path
D=Path("/nfs/lambda_stor_01/homes/avasan/.loom/mcp/hypothesis2omics/mcp/data")
W=Path("/nfs/lambda_stor_01/data/avasan/omics_agents/hepb_nfkb/work")
base=list(csv.DictReader(open(W/"analysis_table.tsv"),delimiter="\t"))
titer=np.array([float(r["antiHBs_post"]) for r in base])
meta={}
for r in csv.DictReader(open(D/"immport_cache/parsed/sample_manifest.tsv"),delimiter="\t"):
    if r["study_accession"]=="SDY1328" and r["repository_accession"].startswith("GSM"):
        meta[r["repository_accession"]]=(float(r["min_subject_age_in_years"]),r["gender"])
age=np.array([meta[r["gsm"]][0] for r in base])
rt=stats.rankdata(titer); ra=stats.rankdata(age)
def resid(a,b):
    b1=np.c_[np.ones(len(b)),b]; return a-b1@np.linalg.lstsq(b1,a,rcond=None)[0]
rt_r=resid(rt,ra)
def partial(x):
    r,p=stats.pearsonr(resid(stats.rankdata(x),ra),rt_r)
    return float(r),float(p)
def ci(r,n,k=1):
    z=np.arctanh(r); se=1/math.sqrt(n-3-k); return float(np.tanh(z-1.96*se)),float(np.tanh(z+1.96*se))

sc=json.load(open(W/"composite_scores.json"))
PRIMARY="TNF-alpha Signaling via NF-kB"
out={"n":len(base),"primary_set":PRIMARY,"sets":{}}
print(f"{'gene set':34}{'n':>5}{'rho raw':>10}{'p raw':>9}{'rho|age':>10}{'p|age':>9}{'95% CI (adj)':>22}")
for name,d in sc.items():
    x=np.array(d["score"]); rr,pr=stats.spearmanr(x,titer); pa,ppa=partial(x)
    lo,hi=ci(pa,len(base))
    out["sets"][name]={"n_genes":d["n_genes"],"rho_raw":float(rr),"p_raw":float(pr),
                       "rho_adj":pa,"p_adj":ppa,"ci_lo":lo,"ci_hi":hi}
    print(f"{name:34}{d['n_genes']:5}{rr:10.3f}{pr:9.4f}{pa:10.3f}{ppa:9.4f}   [{lo:+.3f},{hi:+.3f}]")
# BH across the 4 sets on the age-adjusted p
names=list(out["sets"]); ps=np.array([out["sets"][n]["p_adj"] for n in names])
o=np.argsort(ps); m=len(ps); q=np.empty(m); prev=1.0
for i in range(m-1,-1,-1):
    j=o[i]; v=ps[j]*m/(i+1); prev=min(prev,v); q[j]=min(prev,1.0)
for n,qq in zip(names,q): out["sets"][n]["q_BH"]=float(qq)
p=out["sets"][PRIMARY]
print(f"\nPRIMARY -> {PRIMARY}: age-adjusted rho={p['rho_adj']:+.3f}, p={p['p_adj']:.4f}, q={p['q_BH']:.4f}")
json.dump(out,open(W/"results_sets.json","w"),indent=2)
