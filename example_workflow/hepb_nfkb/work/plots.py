import csv,json,base64,io,numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
D=Path("/nfs/lambda_stor_01/homes/avasan/.loom/mcp/hypothesis2omics/mcp/data")
rows=list(csv.DictReader(open("work/analysis_table.tsv"),delimiter="\t"))
res=json.load(open("work/results.json"))
meta={}
for r in csv.DictReader(open(D/"immport_cache/parsed/sample_manifest.tsv"),delimiter="\t"):
    if r["study_accession"]=="SDY1328" and r["repository_accession"].startswith("GSM"):
        meta[r["repository_accession"]]=float(r["min_subject_age_in_years"])
age=np.array([meta[r["gsm"]] for r in rows])
titer=np.array([float(r["antiHBs_post"]) for r in rows]); lt=np.log2(titer)
GEN=["NFKB1","RELA","NFKB2","REL","RELB"]
def b64(fig):
    b=io.BytesIO(); fig.savefig(b,format="png",dpi=110,bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(b.getvalue()).decode()
imgs={}

# 1: NFKB1 scatter (primary)
f,a=plt.subplots(figsize=(6,4.4))
x=np.array([float(r["NFKB1"]) for r in rows])
a.scatter(x,lt,c=np.where(titer>=10,"#2166ac","#b2182b"),s=26,alpha=.8,edgecolor="none")
m,b_=np.polyfit(x,lt,1); xs=np.linspace(x.min(),x.max(),50); a.plot(xs,m*xs+b_,"k--",lw=1.4)
a.axhline(np.log2(10),color="gray",ls=":",lw=1.2)
a.text(x.min(),np.log2(10)+.25,"seroprotection 10 mIU/mL",fontsize=7.5,color="gray")
r0=res["results"][0]
a.set_xlabel("Baseline NFKB1 expression (log2, RMA)"); a.set_ylabel("log2 anti-HBs (mIU/mL)")
a.set_title(f"Primary test: NFKB1 vs anti-HBs\nrho={r0['rho']:.3f}  p={r0['p']:.3f}  q={r0['q_BH']:.3f}  n=165",fontsize=10)
imgs["scatter"]=b64(f)

# 2: forest of rho with CI, raw vs age-adjusted
f,a=plt.subplots(figsize=(6.4,4))
y=np.arange(len(GEN))
for i,r in enumerate(res["results"]):
    a.plot([r["ci_lo"],r["ci_hi"]],[i,i],color="#4393c3",lw=2.4)
    a.plot(r["rho"],i,"o",color="#2166ac",ms=7,label="unadjusted" if i==0 else "")
    a.plot(r["partial_rho_age"],i,"D",color="#d6604d",ms=6,label="age-adjusted" if i==0 else "")
a.axvline(0,color="k",lw=1); a.set_yticks(y); a.set_yticklabels(GEN)
a.set_xlabel("Spearman rho (negative = hypothesis direction)")
a.set_title("NF-\u03baB genes vs anti-HBs: effect sizes with 95% CI",fontsize=10)
a.legend(fontsize=8,loc="lower right"); a.invert_yaxis()
imgs["forest"]=b64(f)

# 3: responder vs non-responder boxplots
f,axs=plt.subplots(1,5,figsize=(12,3.4),sharey=False)
for ax,g in zip(axs,GEN):
    v=np.array([float(r[g]) for r in rows]); hi=v[titer>=10]; lo=v[titer<10]
    bp=ax.boxplot([lo,hi],labels=["non\n(<10)","resp\n(>=10)"],widths=.6,patch_artist=True)
    for p,c in zip(bp["boxes"],["#b2182b","#2166ac"]): p.set_facecolor(c); p.set_alpha(.45)
    rr=[q for q in res["results"] if q["gene"]==g][0]
    ax.set_title(f"{g}\nMWU p={rr['mwu_p']:.3f}",fontsize=9); ax.tick_params(labelsize=8)
axs[0].set_ylabel("baseline expression (log2)")
f.suptitle("Baseline NF-\u03baB by response status (n=107 non-responders / 58 responders)",fontsize=10)
imgs["box"]=b64(f)

# 4: the confounder
f,axs=plt.subplots(1,2,figsize=(9.5,3.8))
axs[0].scatter(age,lt,s=24,c="#762a83",alpha=.75,edgecolor="none")
m,b_=np.polyfit(age,lt,1); xs=np.linspace(age.min(),age.max(),50); axs[0].plot(xs,m*xs+b_,"k--",lw=1.4)
av=res["age_vs_titer"]; axs[0].set_title(f"Age vs anti-HBs\nrho={av['rho']:.3f}  p={av['p']:.4f}",fontsize=10)
axs[0].set_xlabel("Age (years)"); axs[0].set_ylabel("log2 anti-HBs")
axs[1].scatter(age,x,s=24,c="#1b7837",alpha=.75,edgecolor="none")
m,b_=np.polyfit(age,x,1); axs[1].plot(xs,m*xs+b_,"k--",lw=1.4)
axs[1].set_title(f"Age vs baseline NFKB1\nrho={r0['age_vs_expr_rho']:.3f}  p={r0['age_vs_expr_p']:.4f}",fontsize=10)
axs[1].set_xlabel("Age (years)"); axs[1].set_ylabel("NFKB1 (log2)")
f.suptitle("Age confounds both sides of the hypothesis",fontsize=11)
imgs["confound"]=b64(f)

json.dump(imgs,open("work/plots.json","w"))
print("plots:",{k:len(v)//1024 for k,v in imgs.items()},"KB")
