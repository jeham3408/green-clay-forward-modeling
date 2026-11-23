#!/usr/bin/env python3
"""15/9-13: log profile through the Hordaland clay (20 ms bins) and what the Green Clay model predicts there.
Read-only."""
import json
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026"
W = json.load(open(R + "/sleipner/cubes/sleipner94/wells.json"))
w = W["wells"][0]
print("tops:", json.dumps(w["tops"])[:900])
print("gc_sites:", json.dumps(W.get("gc_sites"))[:900])
L = w["logs"]
n = L["n"]
t = L["t0"] + L["dt"] * np.arange(n)


def f(k):
    return np.array([np.nan if v is None else v for v in L[k]], float)


md, gr, vp, rho, rhob, ai = f("md"), f("gr"), f("vp"), f("rho_used"), f("rhob"), f("ai")
print("log samples", n, "TWT", t[0], t[-1], "md", np.nanmin(md), np.nanmax(md))
ok = np.isfinite(rhob)
print("density measured from TWT", t[ok][0] if ok.any() else None, "to", t[ok][-1] if ok.any() else None,
      "; share measured in 1100-1990 ms:", round(float(np.mean(ok[(t >= 1100) & (t <= 1990)])), 3))
print(" TWT ms |  MD m |  Vp m/s | rho used | rho measured |   GR  |   AI")
for a in range(780, 1990, 20):
    m = (t >= a) & (t < a + 20)
    if not m.any() or not np.isfinite(vp[m]).any():
        continue
    rm = f"{np.nanmedian(rhob[m]):.3f}" if np.isfinite(rhob[m]).any() else "gardner"
    print(f"{a:6d}   {np.nanmedian(md[m]):6.0f}   {np.nanmedian(vp[m]):6.0f}    {np.nanmedian(rho[m]):5.3f}    {rm:>8s}     {np.nanmedian(gr[m]):6.1f}   {np.nanmedian(ai[m]):6.0f}")

G = json.load(open(R + "/natures-plug/data/greenclay.json"))
ex = G["ex"]
names = ex["n"]
for i, nm in enumerate(names):
    if nm in ("15/9-13", "15/9-A-16", "15/9-11", "15/9-9", "15/9-19 A", "15/9-19 SR"):
        print("greenclay model:", {k: ex[k][i] for k in ex if isinstance(ex[k], list) and len(ex[k]) == len(names)})
print("rule:", {k: v for k, v in G["rule"].items() if not isinstance(v, (list, dict))})
L1 = json.load(open(R + "/natures-plug/data/layer1.json"))
items = L1 if isinstance(L1, list) else L1.get("wells", L1.get("ex", []))
if isinstance(items, list):
    for it in items:
        if isinstance(it, dict) and str(it.get("n", it.get("name", ""))).startswith("15/9-13"):
            print("layer1:", {k: v for k, v in it.items() if not isinstance(v, (list, dict)) or len(str(v)) < 80})
else:
    print("layer1 structure:", type(L1).__name__, list(L1.keys())[:12] if isinstance(L1, dict) else "")
