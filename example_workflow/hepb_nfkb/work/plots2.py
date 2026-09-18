import csv,json,base64,io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
W=Path("work"); D=Path("/nfs/lambda_stor_01/homes/avasan/.loom/mcp/hypothesis2omics/mcp/data")
imgs=json.load(open(W/"plots.json"))
rs=json.load(open(W/"results_sets.json")); sw=json.load(open(W/"sweep.json"))
rep=json.load(open(W/"replication.json")); perm=json.load(open(W/"enrichment_perm.json"))
sc=json.load(open(W/"composite_scores.json"))
base=list(csv.DictReader(open(W/"analysis_table.tsv"),delimiter="\t"))
titer=np.array([float(r["antiHBs_post"]) for r in base]); lt=np.log2(titer)
def b64(f):
    b=io.BytesIO(); f.savefig(b,format="png",dpi=110,bbox_inches="tight"); plt.close(f)
    return base64.b64encode(b.getvalue()).decode()
SETS=list(rs["sets"])

# A: composite scatter (primary set)
P=rs["primary_set"]; x=np.array(sc[P]["score"])
f,a=plt.subplots(figsize=(6,4.3))
a.scatter(x,lt,c=np.where(titer>=10,"#2166ac","#b2182b"),s=26,alpha=.8,edgecolor="none")
m,b_=np.polyfit(x,lt,1); xs=np.linspace(x.min(),x.max(),50); a.plot(xs,m*xs+b_,"k--",lw=1.4)
a.axhline(np.log2(10),color="gray",ls=":",lw=1.2)
d=rs["sets"][P]
a.set_xlabel(f"Baseline {P} composite z-score ({d['n_genes']} genes)")
a.set_ylabel("log2 anti-HBs (mIU/mL)")
a.set_title(f"PRIMARY: inflammation composite vs anti-HBs\nage-adj rho={d['rho_adj']:+.3f}  p={d['p_adj']:.3f}  q={d['q_BH']:.3f}  n=165",fontsize=10)
imgs["comp_scatter"]=b64(f)

# B: forest, 4 sets raw vs adjusted
f,a=plt.subplots(figsize=(7,3.6))
for i,s in enumerate(SETS):
    d=rs["sets"][s]
    a.plot([d["ci_lo"],d["ci_hi"]],[i,i],color="#4393c3",lw=2.4)
    a.plot(d["rho_raw"],i,"o",color="#2166ac",ms=7,label="unadjusted" if i==0 else "")
    a.plot(d["rho_adj"],i,"D",color="#d6604d",ms=6,label="age-adjusted" if i==0 else "")
a.axvline(0,color="k",lw=1); a.set_yticks(range(len(SETS)))
a.set_yticklabels([s.replace(" Response","")+f"\n({rs['sets'][s]['n_genes']}g)" for s in SETS],fontsize=8)
a.set_xlabel("Spearman rho (negative = hypothesis direction)")
a.set_title("Four inflammation gene-set composites vs anti-HBs (95% CI, age-adjusted)",fontsize=10)
a.legend(fontsize=8,loc="lower right"); a.invert_yaxis()
imgs["comp_forest"]=b64(f)

# C: genome-wide sweep -- p histogram + QQ
pv=np.array(sw["p"]); m=len(pv)
f,axs=plt.subplots(1,2,figsize=(10,3.6))
axs[0].hist(pv,bins=50,color="#6b7280",edgecolor="none")
axs[0].axhline(m/50,color="#b2182b",ls="--",lw=1.5,label="uniform (null)")
axs[0].set_xlabel("p-value (age-adjusted)"); axs[0].set_ylabel("genes")
axs[0].set_title(f"All {m:,} genes: p-value distribution\n{sw['n_p05']} at p<0.05 vs {0.05*m:.0f} expected",fontsize=9.5)
axs[0].legend(fontsize=8)
o=np.sort(pv); exp=-np.log10((np.arange(1,m+1)-.5)/m)
axs[1].plot(exp,-np.log10(o),".",ms=2,color="#2166ac")
lim=max(exp.max(),(-np.log10(o)).max())
axs[1].plot([0,lim],[0,lim],"r--",lw=1.2)
axs[1].set_xlabel("expected -log10(p)"); axs[1].set_ylabel("observed")
axs[1].set_title(f"QQ plot: no departure from null\n0 genes at q<0.05",fontsize=9.5)
f.suptitle("Exploratory genome-wide sweep finds no transcriptional signal",fontsize=11)
imgs["sweep"]=b64(f)

# D: permutation enrichment + replication
f,axs=plt.subplots(1,2,figsize=(10,3.6))
obs=perm["observed"]; pp=perm["p"]
xs=np.arange(len(SETS)); expd=[0.05*rs["sets"][s]["n_genes"] for s in SETS]
axs[0].bar(xs-.2,[obs[s] for s in SETS],.4,label="observed",color="#2166ac")
axs[0].bar(xs+.2,expd,.4,label="expected",color="#9ca3af")
for i,s in enumerate(SETS): axs[0].text(i,max(obs[s],expd[i])+.8,f"p={pp[s]:.3f}",ha="center",fontsize=8)
axs[0].set_xticks(xs); axs[0].set_xticklabels([s.replace(" Response","").replace("TNF-alpha Signaling via ","")[:14] for s in SETS],fontsize=8,rotation=15)
axs[0].set_ylabel("genes in top 5%"); axs[0].legend(fontsize=8)
axs[0].set_title("Set enrichment vs 2,000 label permutations\n(preserves gene-gene correlation)",fontsize=9.5)
w=.35
axs[1].bar(xs-w/2,[rep[s]["set1"][0] for s in SETS],w,label="donor set 1 (n=121)",color="#2166ac")
axs[1].bar(xs+w/2,[rep[s]["set2"][0] for s in SETS],w,label="donor set 2 (n=44)",color="#d6604d")
axs[1].axhline(0,color="k",lw=1); axs[1].set_xticks(xs)
axs[1].set_xticklabels([s.replace(" Response","").replace("TNF-alpha Signaling via ","")[:14] for s in SETS],fontsize=8,rotation=15)
axs[1].set_ylabel("age-adjusted rho"); axs[1].legend(fontsize=8)
axs[1].set_title("Held-out replication: sign-consistent but null in both",fontsize=9.5)
imgs["perm_rep"]=b64(f)
json.dump(imgs,open(W/"plots.json","w"))
print("plots now:",list(imgs))
