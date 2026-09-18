import json, datetime
from pathlib import Path
steps=json.load(open("work/steps.json"))
PH=["Setup","Parse","Search","Gate","Fetch","Prepare","Verify","Analyse","Synthesise","Report"]
COL={"Setup":"#6b7280","Parse":"#7c3aed","Search":"#0891b2","Gate":"#f0a202","Fetch":"#0e7490",
     "Prepare":"#2563eb","Verify":"#059669","Analyse":"#dc2626","Synthesise":"#9333ea","Report":"#1b7837"}
def esc(s): return (s.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;"))
nodes=[]
for s in steps:
    c=COL[s["phase"]]; dec = s["kind"]=="decision"
    nodes.append(f'''<div class="node {'dec' if dec else 'tool'}" data-id="{s['id']}" style="--c:{c}">
      <div class="idx">{s['id']}</div>
      <div class="body"><div class="ph" style="background:{c}">{s['phase']}</div>
      <div class="ti">{esc(s['title'])}</div>
      {'<div class="badge">DECISION POINT</div>' if dec else '<div class="badge tl">tool call</div>'}</div></div>''')
det={}
for s in steps:
    det[s["id"]]=f'''<h3><span class="ph" style="background:{COL[s['phase']]}">{s['phase']}</span> Step {s['id']}: {esc(s['title'])}</h3>
<div class="sec"><h4>Chat at this point</h4><p class="chat">{esc(s['chat'])}</p></div>
<div class="sec"><h4>Tools planned</h4><ul>{''.join(f'<li>{esc(x)}</li>' for x in s['planned'])}</ul></div>
<div class="sec"><h4>Tools actually called</h4><ul>{''.join(f'<li><code>{esc(x)}</code></li>' for x in s['called'])}</ul></div>
<div class="sec dc"><h4>Decision</h4><p><b>{esc(s['decision'])}</b></p>
<h4>Why</h4><p>{esc(s['why'])}</p></div>
<div class="sec oc"><h4>Outcome</h4><p>{esc(s['outcome'])}</p></div>'''
H=f'''<!doctype html><meta charset=utf-8><title>NF-\u03baB / HepB workflow orchestration</title>
<style>
*{{box-sizing:border-box}}
body{{font:15px/1.6 -apple-system,Segoe UI,Roboto,sans-serif;margin:0;background:#f7f8fa;color:#18181b}}
header{{background:#111827;color:#fff;padding:1.3rem 2rem}}
header h1{{margin:0;font-size:1.35rem}} header p{{margin:.35rem 0 0;color:#9ca3af;font-size:13px}}
.wrap{{display:grid;grid-template-columns:minmax(330px,1fr) 1.25fr;gap:1.4rem;padding:1.4rem 2rem;align-items:start}}
.flow{{display:flex;flex-direction:column;gap:0}}
.node{{display:flex;gap:.8rem;background:#fff;border:2px solid #e5e7eb;border-left:6px solid var(--c);
 border-radius:8px;padding:.7rem .9rem;cursor:pointer;transition:.13s;position:relative}}
.node:hover{{transform:translateX(5px);box-shadow:0 4px 14px rgba(0,0,0,.13);border-color:var(--c)}}
.node.active{{border-color:var(--c);box-shadow:0 0 0 3px color-mix(in srgb,var(--c) 25%,transparent)}}
.node.dec{{background:#fffdf5}}
.idx{{width:26px;height:26px;border-radius:50%;background:var(--c);color:#fff;font-weight:700;
 display:flex;align-items:center;justify-content:center;font-size:13px;flex:0 0 auto}}
.ti{{font-weight:600;font-size:14px;margin:.15rem 0}}
.ph{{display:inline-block;color:#fff;font-size:9.5px;font-weight:700;padding:1px 6px;border-radius:3px;
 letter-spacing:.04em;text-transform:uppercase}}
.badge{{font-size:10px;font-weight:700;color:#b45309;background:#fef3c7;display:inline-block;
 padding:1px 6px;border-radius:3px;margin-top:.25rem}}
.badge.tl{{color:#3730a3;background:#e0e7ff}}
.conn{{height:14px;width:2px;background:#d1d5db;margin-left:calc(.9rem + 13px)}}
.panel{{background:#fff;border:1px solid #e5e7eb;border-radius:10px;padding:1.3rem 1.5rem;
 position:sticky;top:1.4rem;max-height:calc(100vh - 3rem);overflow-y:auto}}
.panel h3{{margin:0 0 1rem;font-size:1.08rem;display:flex;align-items:center;gap:.5rem;flex-wrap:wrap}}
.sec{{margin:.9rem 0;padding:.75rem .9rem;background:#f9fafb;border-radius:7px}}
.sec h4{{margin:0 0 .4rem;font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:#6b7280}}
.sec p,.sec li{{margin:.25rem 0;font-size:13.5px}} .sec ul{{margin:.2rem 0;padding-left:1.15rem}}
.sec.dc{{background:#fffbeb;border-left:4px solid #f0a202}}
.sec.oc{{background:#f0fdf4;border-left:4px solid #1b7837}}
.chat{{font-style:italic;color:#374151}}
code{{background:#eef2ff;padding:1px 5px;border-radius:3px;font-size:12.5px;word-break:break-word}}
.hint{{color:#6b7280;font-size:13px;text-align:center;padding:2rem 1rem}}
.legend{{display:flex;gap:1rem;flex-wrap:wrap;padding:0 2rem 1rem;font-size:12px;color:#4b5563}}
.legend span{{display:flex;align-items:center;gap:.35rem}}
.sw{{width:11px;height:11px;border-radius:2px;display:inline-block}}
.verdict{{margin:0 2rem 1.2rem;padding:.9rem 1.2rem;background:#fff8e1;border-left:5px solid #f0a202;border-radius:0 6px 6px 0;font-size:14px}}
</style>
<header><h1>Hypothesis validation workflow &mdash; baseline NF-\u03baB vs Hepatitis B vaccine response</h1>
<p>Living document &middot; {datetime.date.today().isoformat()} &middot; hover or click any step for the chat, the tools planned vs called, and the conclusion</p></header>
<div class=verdict><b>Final verdict: INCONCLUSIVE.</b> All 5 NF-\u03baB genes trend as hypothesised
(NFKB1 &rho;=&minus;0.118) but none is significant (q=0.237, n=165); ~half the raw association is
explained by age. <b>{sum(1 for s in steps if s['kind']=='decision')} decision points</b> are highlighted below.</div>
<div class=legend>{''.join(f'<span><i class=sw style="background:{COL[p]}"></i>{p}</span>' for p in PH)}
<span><i class=sw style="background:#fffdf5;border:1px solid #f0a202"></i>cream = decision point</span></div>
<div class=wrap>
 <div class=flow>{'<div class=conn></div>'.join(nodes)}</div>
 <div class=panel id=panel><div class=hint>&#8592; Hover or click a step to inspect it.</div></div>
</div>
<script>
const D={json.dumps(det)};
const panel=document.getElementById('panel');
let locked=null;
function show(id){{ panel.innerHTML=D[id];
  document.querySelectorAll('.node').forEach(n=>n.classList.toggle('active',n.dataset.id==id)); }}
document.querySelectorAll('.node').forEach(n=>{{
  n.addEventListener('mouseenter',()=>{{ if(!locked) show(n.dataset.id); }});
  n.addEventListener('click',()=>{{ locked = locked===n.dataset.id ? null : n.dataset.id; show(n.dataset.id); }});
}});
</script>'''
Path("workflow.html").write_text(H)
print("workflow.html",len(H),"bytes")
