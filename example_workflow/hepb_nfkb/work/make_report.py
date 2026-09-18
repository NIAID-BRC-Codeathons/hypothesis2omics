import json,datetime
from pathlib import Path
res=json.load(open("work/results.json")); imgs=json.load(open("work/plots.json"))
syn=json.load(open("work/synthesis.json")); co=res["cohort"]
G=["NFKB1","RELA","NFKB2","REL","RELB"]
by={r["gene"]:r for r in res["results"]}; p=by["NFKB1"]
def row(r):
    sig="yes" if r["q_BH"]<0.05 else "no"
    return (f"<tr><td class=g>{r['gene']}{' <span class=pri>primary</span>' if r['gene']=='NFKB1' else ''}</td>"
            f"<td>{r['rho']:+.3f}</td><td>[{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}]</td>"
            f"<td>{r['p']:.4f}</td><td>{r['q_BH']:.4f}</td>"
            f"<td>{r['partial_rho_age']:+.3f}</td><td>{r['partial_p_age']:.4f}</td>"
            f"<td class='{'y' if sig=='yes' else 'n'}'>{sig}</td></tr>")
H=f"""<!doctype html><meta charset=utf-8>
<title>NF-\u03baB baseline vs Hepatitis B vaccine response</title>
<style>
body{{font:15px/1.65 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1000px;margin:2rem auto;padding:0 1.2rem;color:#1a1a1a}}
h1{{font-size:1.65rem;margin-bottom:.2rem}} h2{{margin-top:2.2rem;border-bottom:2px solid #eee;padding-bottom:.3rem;font-size:1.2rem}}
.sub{{color:#666;margin-top:0}} table{{border-collapse:collapse;width:100%;margin:1rem 0;font-size:14px}}
th,td{{border:1px solid #ddd;padding:.45rem .6rem;text-align:right}} th{{background:#f5f7fa;text-align:right}}
td.g,th:first-child{{text-align:left;font-weight:600}} .y{{color:#1b7837;font-weight:700}} .n{{color:#b2182b}}
.pri{{background:#2166ac;color:#fff;font-size:10px;padding:1px 5px;border-radius:3px;margin-left:5px;font-weight:600}}
.verdict{{background:#fff8e1;border-left:5px solid #f0a202;padding:1rem 1.2rem;margin:1.4rem 0;border-radius:0 5px 5px 0}}
.verdict b{{font-size:1.1rem}} img{{max-width:100%;border:1px solid #e3e3e3;border-radius:5px;margin:.7rem 0}}
.kv{{display:grid;grid-template-columns:200px 1fr;gap:.3rem .9rem;font-size:14px}} .kv div:nth-child(odd){{color:#666}}
code{{background:#f3f4f6;padding:1px 5px;border-radius:3px;font-size:13px}}
.warn{{background:#fdecea;border-left:5px solid #b2182b;padding:.9rem 1.1rem;margin:1.2rem 0;border-radius:0 5px 5px 0}}
.ok{{background:#edf7ed;border-left:5px solid #1b7837;padding:.9rem 1.1rem;margin:1.2rem 0;border-radius:0 5px 5px 0}}
</style>
<h1>Baseline NF-\u03baB expression and Hepatitis B vaccine response</h1>
<p class=sub>Hypothesis test report &middot; generated {datetime.date.today().isoformat()} &middot; ImmPort SDY1328 / GEO GSE65834</p>

<div class=verdict><b>Verdict: INCONCLUSIVE &mdash; hypothesis not supported by this cohort.</b><br>
All five NF-\u03baB genes trend in the hypothesised direction (higher baseline &rarr; lower antibody),
but none reaches significance. Primary gene <b>NFKB1: &rho; = {p['rho']:+.3f}, p = {p['p']:.3f},
q<sub>BH</sub> = {p['q_BH']:.3f}</b> (n = {res['n']}). After adjusting for age the effect roughly halves
to &rho; = {p['partial_rho_age']:+.3f} (p = {p['partial_p_age']:.3f}).</div>

<h2>1. Hypothesis and pre-specified test</h2>
<p><i>"High level of NFKB, a transcription factor mediating inflammation, before vaccination is
associated with poor response to Hepatitis B vaccine."</i></p>
<div class=kv>
<div>Predictor</div><div>Baseline (day 0, pre-vaccination) whole-blood expression of NF-\u03baB subunits</div>
<div>Primary gene</div><div><code>NFKB1</code> &mdash; pre-specified to avoid multiple-testing laundering</div>
<div>Secondary genes</div><div><code>RELA</code>, <code>NFKB2</code>, <code>REL</code>, <code>RELB</code></div>
<div>Outcome</div><div>Post-vaccination anti-HBs (mIU/mL), ImmPort <code>VALUE_PREFERRED</code></div>
<div>Test</div><div>Spearman rank correlation, <b>two-sided</b> (a positive association stays detectable)</div>
<div>Multiplicity</div><div>Benjamini-Hochberg across the 5 genes, &alpha; = 0.05</div>
<div>Secondary outcome</div><div>Non-responder contrast at the 10 mIU/mL seroprotection threshold (Mann-Whitney)</div>
</div>
<p><b>The decision rule was frozen before any result existed</b> (<code>evidence_rules.DecisionRule</code>,
direction <code>down</code>, &alpha; 0.05, &ge;1 independent group), so alpha could not be widened
nor the direction flipped after seeing the numbers.</p>

<h2>2. Cohort</h2>
<div class=kv>
<div>ImmPort study</div><div>SDY1328 &mdash; arm {co['arm']} (single arm, no treatment confound)</div>
<div>GEO series</div><div>{co['gse_baseline']} "Transcriptional profiling of HBV-na&iuml;ve subjects <b>before</b> vaccination"</div>
<div>Platform</div><div>{co['platform']} (Affymetrix), {co['normalization']}-normalised, 60,607 probes</div>
<div>Tissue / timepoint</div><div>{co['tissue']}, {co['timepoint']}</div>
<div>Vaccine</div><div>{co['vaccine']}</div>
<div>Subjects analysed</div><div><b>{res['n']}</b> with baseline expression + outcome
({res['n_responders']} responders &ge;10, {res['n_nonresponders']} non-responders &lt;10 mIU/mL)</div>
<div>Age</div><div>{res['age_summary']['min']:.0f}&ndash;{res['age_summary']['max']:.0f} y, mean {res['age_summary']['mean']}</div>
</div>
<div class=ok><b>Outcome cross-validated.</b> ImmPort ELISA <code>VALUE_PREFERRED</code> and the
independent GEO sample characteristics agree <b>exactly on all 165 subjects</b> (0 mismatches),
so the response variable is not an artefact of one metadata source.</div>

<h2>3. Primary result</h2>
<img src="data:image/png;base64,{imgs['scatter']}">
<p>Each point is one subject. Blue = seroprotected (&ge;10 mIU/mL), red = non-responder.
The fitted slope is slightly negative, consistent with the hypothesis, but the scatter is dominated
by subjects spanning three orders of magnitude of antibody at essentially identical NFKB1 levels.</p>

<h2>4. All five NF-\u03baB genes</h2>
<table><tr><th>Gene</th><th>&rho;</th><th>95% CI</th><th>p</th><th>q<sub>BH</sub></th>
<th>&rho; | age</th><th>p | age</th><th>q&lt;0.05</th></tr>
{''.join(row(by[g]) for g in G)}</table>
<img src="data:image/png;base64,{imgs['forest']}">
<p>Every confidence interval crosses zero. The consistent negative sign across five genes is
mild corroborating evidence of direction, but the five are co-regulated members of one pathway,
so they are <b>not five independent tests</b> &mdash; the synthesis module counts them as one group.</p>

<h2>5. Responder vs non-responder</h2>
<img src="data:image/png;base64,{imgs['box']}">
<p>Dichotomising at the 10 mIU/mL clinical seroprotection threshold reproduces the same
null: the largest separation is RELB (Mann-Whitney p = {by['RELB']['mwu_p']:.3f}), and group means
differ in the third decimal place.</p>

<h2>6. Age confounds both sides of the hypothesis</h2>
<img src="data:image/png;base64,{imgs['confound']}">
<div class=warn><b>This is the substantive finding.</b> Age predicts antibody response
(&rho; = {res['age_vs_titer']['rho']:+.3f}, p = {res['age_vs_titer']['p']:.4f}) <i>and</i> baseline NFKB1
(&rho; = {p['age_vs_expr_rho']:+.3f}, p = {p['age_vs_expr_p']:.4f}) &mdash; older subjects have higher
baseline NF-\u03baB and worse responses. Controlling for age shrinks the NFKB1 association from
{p['rho']:+.3f} to {p['partial_rho_age']:+.3f}. <b>About half of the already-non-significant raw
association is attributable to age</b>, consistent with inflammaging rather than a direct NF-\u03baB effect.</div>

<h2>7. Verdict</h2>
<div class=verdict><b>{syn['overall_verdict'].upper()}</b><br>
<pre style="white-space:pre-wrap;font:13px/1.55 ui-monospace,monospace;margin:.6rem 0 0">{syn['narrative']}</pre></div>

<h2>8. Limitations &mdash; what this does and does not show</h2>
<ul>
<li><b>Not a refutation.</b> Failing to reject the null at n=165 is not evidence of no effect.
The CI for NFKB1 ([{p['ci_lo']:+.3f}, {p['ci_hi']:+.3f}]) still admits a true &rho; near &minus;0.27.</li>
<li><b>One cohort, one platform.</b> Five genes collapse to <b>one independent group</b>. No replication.</li>
<li><b>mRNA is a weak proxy for NF-\u03baB activity</b>, which is regulated by nuclear translocation,
not transcript abundance. A negative transcript result does not exclude a pathway-activity effect.</li>
<li><b>Whole blood</b> mixes cell types; a shift in composition can mask or mimic a within-cell effect.</li>
<li><b>Elderly, Twinrix (HepA/B) cohort</b> &mdash; may not generalise to young adults or monovalent HepB.</li>
<li><b>Other candidate studies were not testable:</b> of 5 HepB studies passing the readiness gate,
only SDY1328 had a resolvable GEO series; SDY299/1538/2312/2584 had no linked expression matrices.</li>
</ul>

<h2>9. Provenance</h2>
<div class=kv>
<div>Expression matrix</div><div><code>GSE65834_GPL15048/expression.tsv.gz</code><br>sha256 <code>80df7f81ec66bb36&hellip;</code> &middot; 169 samples &times; 60,607 probes</div>
<div>ImmPort manifest</div><div><code>sample_manifest.tsv</code> sha256 <code>9d2f3262c9dbd1e5&hellip;</code></div>
<div>Analysis table</div><div><code>work/analysis_table.tsv</code> sha256 <code>3e0d68bd78e0a55b&hellip;</code></div>
<div>Results / verdict</div><div><code>work/results.json</code>, <code>work/synthesis.json</code>
(schema-valid, <code>validate_synthesis == []</code>)</div>
<div>Probe mapping</div><div>GPL15048 platform table from <code>GSE65834_family.soft.gz</code>;
multi-probe genes averaged (RELA 4, REL 2, others 1)</div>
</div>
"""
Path("report.html").write_text(H)
print("report.html",len(H),"bytes")
