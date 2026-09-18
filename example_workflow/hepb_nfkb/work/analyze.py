#!/usr/bin/env python3
"""Pre-specified test: baseline NF-kB expression vs post-vaccination anti-HBs.

Primary: NFKB1. Secondary: RELA, NFKB2, REL, RELB. Two-sided Spearman,
BH-FDR across the 5 genes, alpha=0.05. Also a non-responder (<10 mIU/mL)
contrast as the pre-specified secondary outcome.
"""
import csv, json, hashlib, math
import numpy as np
from scipy import stats
from pathlib import Path

W=Path("/nfs/lambda_stor_01/data/avasan/omics_agents/hepb_nfkb/work")
rows=list(csv.DictReader(open(W/"analysis_table.tsv"),delimiter="\t"))
GENES=["NFKB1","RELA","NFKB2","REL","RELB"]
titer=np.array([float(r["antiHBs_post"]) for r in rows])
logt=np.log2(titer)
resp=(titer>=10).astype(int)     # seroprotection threshold
print(f"n={len(rows)}  responders(>=10)={resp.sum()}  non-responders={len(resp)-resp.sum()}")

res=[]
for g in GENES:
    x=np.array([float(r[g]) for r in rows])
    rho,p=stats.spearmanr(x,titer)             # rank-based: log is monotone, same rho
    # Fisher z CI for rho
    n=len(x); z=np.arctanh(rho); se=1/math.sqrt(n-3)
    lo,hi=np.tanh(z-1.96*se),np.tanh(z+1.96*se)
    u,pu=stats.mannwhitneyu(x[resp==1],x[resp==0],alternative="two-sided")
    # rank-biserial effect size
    n1,n0=int(resp.sum()),int((resp==0).sum())
    rb=1-2*u/(n1*n0)
    res.append(dict(gene=g,n=n,rho=rho,p=p,ci_lo=lo,ci_hi=hi,
                    mean_resp=float(x[resp==1].mean()),mean_nonresp=float(x[resp==0].mean()),
                    mwu_p=float(pu),rank_biserial=float(rb)))

# BH across the 5 genes (primary correlation test)
ps=np.array([r["p"] for r in res]); order=np.argsort(ps); m=len(ps)
q=np.empty(m); prev=1.0
for i in range(m-1,-1,-1):
    idx=order[i]; val=ps[idx]*m/(i+1); prev=min(prev,val); q[idx]=min(prev,1.0)
for r,qq in zip(res,q): r["q_BH"]=float(qq)

print(f"\n{'gene':7}{'rho':>8}{'95% CI':>20}{'p':>10}{'q_BH':>9}   {'MWU p':>8}  resp vs non")
for r in res:
    print(f"{r['gene']:7}{r['rho']:8.3f}  [{r['ci_lo']:6.3f},{r['ci_hi']:6.3f}]{r['p']:10.4f}{r['q_BH']:9.4f}   {r['mwu_p']:8.4f}  {r['mean_resp']:.3f}/{r['mean_nonresp']:.3f}")

out=W/"results.json"
json.dump({"n":len(rows),"n_responders":int(resp.sum()),"n_nonresponders":int(len(resp)-resp.sum()),
           "alpha":0.05,"method":"spearman two-sided, BH across 5 genes",
           "primary_gene":"NFKB1","results":res},open(out,"w"),indent=2)
print("\nwrote",out)
