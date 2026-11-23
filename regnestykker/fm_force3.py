"""Scratch control for fm_force2.py: is 'no step at the predicted top, a step at the base' an artefact of testing a predicted depth
against an observed one?  The observed base is tested with the same depth error deliberately added, and the 'largest step nearby'
statistic is compared with control depths in the mudstone above."""
import sys, json, glob, numpy as np
R = sys.argv[1]
g = json.load(open(R + "/natures-plug/data/greenclay.json")); ex = g["ex"]; idx = {n: k for k, n in enumerate(ex["n"])}
print("model: average miss %s m, 9 in 10 within %s m (leave-one-field-out, from the model file)" % (g["rule"].get("mae"), g["rule"].get("p90")))
def read_las(p):
    lines = open(p, errors="replace").read().split("\n")
    c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
    end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
    names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
    well = next(l for l in lines if l.startswith("WELL.")).split(".", 1)[1].split(":")[0].strip()
    arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
    return well, {n: arr[:, k] for k, n in enumerate(names)}
wells = []
for p in sorted(glob.glob(R + "/force2020/las/*/*.las")):
    well, c = read_las(p)
    if well not in idx or "DTC" not in c or "RHOB" not in c: continue
    k = idx[well]; top, th, ho, ro = ex["top"][k], ex["th"][k], ex["ho"][k], ex["ro"][k]
    if top is None or ro is None or ho is None: continue
    md, dtc, rho = c["DEPT"], c["DTC"], c["RHOB"]
    ok = np.isfinite(dtc) & np.isfinite(rho) & (dtc > 40) & (dtc < 250) & (rho > 1.5) & (rho < 3.0)
    if "CALI" in c and "BS" in c:
        wash = c["CALI"] - c["BS"]; ok &= ~(np.isfinite(wash) & (wash > 3.0))
    wells.append(dict(well=well, md=md, ai=304800.0 / dtc * rho, ok=ok, top=top, base=top + th, ho=ho, ro=ro, kn=ex["kn"][k]))
def step(w, z, W):
    def win(a, b):
        m = w["ok"] & (w["md"] >= a) & (w["md"] < b)
        return float(np.median(w["ai"][m])) if m.sum() * 0.152 >= 0.6 * (b - a) else None
    u, d = win(z - W, z), win(z, z + W)
    return None if u is None or d is None else (d - u) / (d + u)
def best(w, z, W, reach=100.0):
    b = None
    for q in np.arange(z - reach, z + reach + 0.1, 5.0):
        r = step(w, q, W)
        if r is not None and (b is None or abs(r) > abs(b)): b = r
    return b
def show(label, vals):
    a = np.array([v for v in vals if v is not None])
    if a.size < 5: print("  %-66s n=%d (too few)" % (label, a.size)); return
    print("  %-66s n=%2d  median %+.3f  16-84 %% %+.3f..%+.3f  |R|>0.03: %2d  negative: %2d" % (label, a.size, np.median(a), *np.percentile(a, [16, 84]), (np.abs(a) > 0.03).sum(), (a < 0).sum()))
for W in (50.0, 30.0):
    print("windows of %.0f m above and below" % W)
    show("observed base (top Rogaland Gp), exact depth", [step(w, w["ro"], W) for w in wells])
    show("observed base moved 55 m up", [step(w, w["ro"] - 55, W) for w in wells])
    show("observed base moved 55 m down", [step(w, w["ro"] + 55, W) for w in wells])
    show("observed base moved 100 m up", [step(w, w["ro"] - 100, W) for w in wells])
    show("predicted top of the clay", [step(w, w["top"], W) for w in wells])
    show("control: middle of the Hordaland Gp above the predicted top", [step(w, 0.5 * (w["ho"] + w["top"]), W) for w in wells if w["top"] - w["ho"] > 2.5 * W])
    show("control: middle of the predicted clay", [step(w, 0.5 * (w["top"] + w["base"]), W) for w in wells])
print("largest step (either sign) within 100 m, 50 m windows")
show("around the predicted top", [best(w, w["top"], 50.0) for w in wells])
show("around the observed base", [best(w, w["ro"], 50.0) for w in wells])
show("control: around the middle of the Hordaland Gp above", [best(w, 0.5 * (w["ho"] + w["top"]), 50.0) for w in wells if w["top"] - w["ho"] > 300])
show("control: around the middle of the predicted clay", [best(w, 0.5 * (w["top"] + w["base"]), 50.0) for w in wells if w["base"] - w["top"] > 250])
print("the same for the", sum(1 for w in wells if w["kn"] <= 10), "wells within 10 km of a report site")
near = [w for w in wells if w["kn"] <= 10]
show("predicted top, 50 m windows", [step(w, w["top"], 50.0) for w in near]); show("largest step within 100 m of the predicted top", [best(w, w["top"], 50.0) for w in near]); show("observed base", [step(w, w["ro"], 50.0) for w in near])
