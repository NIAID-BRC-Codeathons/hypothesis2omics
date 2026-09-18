#!/usr/bin/env python3
"""Exploratory: every gene vs anti-HBs, age-adjusted. BH across all genes."""
import csv,json
import numpy as np
from scipy import stats
from pathlib import Path
D=Path("/nfs/lambda_stor_01/homes/avasan/.loom/mcp/hypothesis2omics/mcp/data")
W=Path("/nfs/lambda_stor_01/data/avasan/omics_agents/hepb_nfkb/work")
X=np.load(W/"gene_matrix.npy"); ix=json.load(open(W/"gene_index.json")); genes=ix["genes"]
base=list(csv.DictReader(open(W/"analysis_table.tsv"),delimiter="\t"))
titer=np.array([float(r["antiHBs_post"]) for r in base])
meta={}
for r in csv.DictReader(open(D/"immport_cache/parsed/sample_manifest.tsv"),delimiter="\t"):
    if r["study_accession"]=="SDY1328" and r["repository_accession"].startswith("GSM"):
        meta[r["repository_accession"]]=float(r["min_subject_age_in_years"])
age=np.array([meta[r["gsm"]] for r in base])
n=len(base)
ra=stats.rankdata(age); rt=stats.rankdata(titer)
def resid_vec(a,b):
    b1=np.c_[np.ones(len(b)),b]; return a-b1@np.linalg.lstsq(b1,a,rcond=None)[0]
rt_r=resid_vec(rt,ra)
# rank each gene across subjects, residualise on age ranks, correlate
R=np.apply_along_axis(stats.rankdata,1,X)
B=np.c_[np.ones(n),ra]; H=B@np.linalg.pinv(B)
Rr=R-R@H.T
rt_c=rt_r-rt_r.mean()
num=Rr@rt_c
den=np.sqrt((Rr**2).sum(1)*(rt_c**2).sum())
rho=np.divide(num,den,out=np.zeros_like(num),where=den>0)
df=n-2-1
t=rho*np.sqrt(df/np.clip(1-rho**2,1e-300,None))
p=2*stats.t.sf(np.abs(t),df)
ok=np.isfinite(p)&(den>0)
pv=p[ok]; gs=[g for g,k in zip(genes,ok) if k]; rv=rho[ok]
m=len(pv); o=np.argsort(pv); q=np.empty(m); prev=1.0
for i in range(m-1,-1,-1):
    j=o[i]; v=pv[j]*m/(i+1); prev=min(prev,v); q[j]=min(prev,1.0)
print(f"genes tested: {m}")
print(f"p<0.05: {(pv<0.05).sum()}  (expected by chance ~{0.05*m:.0f})")
print(f"p<0.01: {(pv<0.01).sum()}  (expected ~{0.01*m:.0f})")
print(f"q<0.05: {(q<0.05).sum()}   q<0.10: {(q<0.10).sum()}")
print(f"\ntop 15 by p (age-adjusted):")
print(f"  {'gene':12}{'rho':>8}{'p':>11}{'q_BH':>8}")
for j in o[:15]: print(f"  {gs[j]:12}{rv[j]:8.3f}{pv[j]:11.2e}{q[j]:8.3f}")
# where do the hallmark sets' genes rank?
sets=json.load(open("genesets/sets_resolved.json"))
rank={g:i for i,j in enumerate(o) for g in [gs[j]]}
print("\nhallmark set members among the strongest signals:")
for name,d in sets.items():
    mem=[g for g in d["on_platform"] if g in rank]
    top=sum(1 for g in mem if rank[g]<m*0.05)
    print(f"  {name:34} {top:4}/{len(mem)} in top 5% (expected {0.05*len(mem):.1f})")
json.dump({"n_tested":m,"genes":gs,"rho":rv.tolist(),"p":pv.tolist(),"q":q.tolist(),
           "n_q05":int((q<0.05).sum()),"n_p05":int((pv<0.05).sum())},
          open(W/"sweep.json","w"))
