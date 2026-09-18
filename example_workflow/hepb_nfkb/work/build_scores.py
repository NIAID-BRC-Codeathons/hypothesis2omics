#!/usr/bin/env python3
"""Composite inflammation scores + full gene-level matrix for GSE65834 baseline.

PRE-SPECIFIED (frozen before any result is inspected):
  Primary   : HALLMARK TNF-alpha Signaling via NF-kB composite z-score,
              AGE-ADJUSTED partial Spearman vs anti-HBs. ONE test.
  Secondary : 3 other hallmark sets (Inflammatory Response, IFN-a, IFN-g),
              BH across the 4 sets.
  Exploratory: every gene on the platform, BH across all tested genes.
Composite = mean of per-gene z-scores across subjects (genes standardised first,
so high-variance probes cannot dominate). Multi-probe genes averaged first.
"""
import csv, gzip, json, hashlib
import numpy as np
from pathlib import Path

D=Path("/nfs/lambda_stor_01/homes/avasan/.loom/mcp/hypothesis2omics/mcp/data")
W=Path("/nfs/lambda_stor_01/data/avasan/omics_agents/hepb_nfkb/work")
EXPR=D/"geo_cache/parsed/SDY1328/EXP16735/GSE65834_GPL15048/expression.tsv.gz"

base=list(csv.DictReader(open(W/"analysis_table.tsv"),delimiter="\t"))
keep=[r["gsm"] for r in base]                      # the validated n=165
print(f"subjects carried over: {len(keep)}")

# probe -> symbol
sym={}
for r in csv.DictReader(open(W/"GPL15048_annot.tsv"),delimiter="\t"):
    g=(r.get("GeneSymbol") or "").strip().upper()
    if g: sym[r["ID"]]=g

# read matrix, collapse probes -> gene means, restricted to our 165 columns
with gzip.open(EXPR,"rt") as fh:
    rd=csv.reader(fh,delimiter="\t"); hdr=next(rd)
    cols=hdr[1:]; idx=[cols.index(g) for g in keep]
    acc={}
    for row in rd:
        g=sym.get(row[0])
        if not g: continue
        try: v=np.array([float(row[1+i]) for i in idx])
        except ValueError: continue
        a=acc.get(g)
        acc[g]=v if a is None else a+v
        acc[g+"__n"]=acc.get(g+"__n",0)+1
genes=sorted(k for k in acc if not k.endswith("__n"))
X=np.vstack([acc[g]/acc[g+"__n"] for g in genes])   # genes x subjects
print(f"gene-level matrix: {X.shape[0]} genes x {X.shape[1]} subjects")

# z-score each gene across subjects
mu=X.mean(1,keepdims=True); sd=X.std(1,ddof=1,keepdims=True)
ok=(sd[:,0]>0)
Z=np.zeros_like(X); Z[ok]=(X[ok]-mu[ok])/sd[ok]
gi={g:i for i,g in enumerate(genes)}

sets=json.load(open("genesets/sets_resolved.json"))
scores={}
for name,d in sets.items():
    mem=[g for g in d["on_platform"] if g in gi and ok[gi[g]]]
    scores[name]={"n_genes":len(mem),"score":Z[[gi[g] for g in mem]].mean(0).tolist()}
    print(f"  {name}: composite over {len(mem)} genes")

np.save(W/"gene_matrix.npy",X)
json.dump({"genes":genes,"gsms":keep},open(W/"gene_index.json","w"))
json.dump(scores,open(W/"composite_scores.json","w"))
print("saved gene_matrix.npy, gene_index.json, composite_scores.json")
