"""Scratch: is there an impedance step at the predicted top and base of the Green Clay in the open FORCE wells? (read-only)"""
import sys, json, glob, numpy as np
R = sys.argv[1]
g = json.load(open(R + "/natures-plug/data/greenclay.json"))["ex"]; ex = {n: k for k, n in enumerate(g["n"])}
def read_las(p):
    lines = open(p, errors="replace").read().split("\n")
    c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
    end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
    names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
    well = next(l for l in lines if l.startswith("WELL.")).split(".", 1)[1].split(":")[0].strip()
    arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
    return well, {n: arr[:, k] for k, n in enumerate(names)}
W = 50.0
out = []
for p in sorted(glob.glob(R + "/force2020/las/*/*.las")):
    well, c = read_las(p)
    if well not in ex: continue
    k = ex[well]; top, th, ho, ro, kn = g["top"][k], g["th"][k], g["ho"][k], g["ro"][k], g["kn"][k]
    if top is None or ro is None or "DTC" not in c or "RHOB" not in c: continue
    md, dtc, rho = c["DEPT"], c["DTC"], c["RHOB"]
    ok = np.isfinite(dtc) & np.isfinite(rho) & (dtc > 40) & (dtc < 250) & (rho > 1.5) & (rho < 3.0)
    if "CALI" in c and "BS" in c:
        wash = c["CALI"] - c["BS"]; ok &= ~(np.isfinite(wash) & (wash > 3.0))          # badly washed-out hole: density not trusted
    ai = 304800.0 / dtc * rho
    def win(a, b):
        m = ok & (md >= a) & (md < b)
        return (float(np.median(ai[m])), float(np.median(304800.0 / dtc[m])), float(np.median(rho[m]))) if m.sum() * 0.152 >= 0.6 * (b - a) else None
    rec = dict(well=well, kn=kn, top=top, base=top + th)
    a, b = win(top - W, top), win(top, top + W)
    if a and b: rec["r_top"] = (b[0] - a[0]) / (b[0] + a[0])
    a2, b2 = win(top + th - W, top + th), win(top + th, top + th + W)
    if a2 and b2: rec["r_base"] = (b2[0] - a2[0]) / (b2[0] + a2[0])
    # the largest 50 m-window impedance step within 100 m of the predicted top (either sign)
    best = None
    for z in np.arange(top - 100, top + 100.1, 5.0):
        u, d = win(z - W, z), win(z, z + W)
        if u and d:
            r = (d[0] - u[0]) / (d[0] + u[0])
            if best is None or abs(r) > abs(best[1]): best = (z, r)
    if best: rec["r_best"] = best
    inside, above = win(top, top + th), win(max(ho, top - th), top)
    if inside and above: rec["inside"], rec["above"] = inside, above
    out.append(rec)
rt = np.array([r["r_top"] for r in out if "r_top" in r]); rb = np.array([r["r_base"] for r in out if "r_base" in r]); rbest = np.array([r["r_best"][1] for r in out if "r_best" in r])
print("FORCE wells with sonic and density around the predicted Green Clay: top %d, base %d" % (rt.size, rb.size))
q = lambda a: "median %+.3f, 16-84 %% %+.3f..%+.3f, |R|>0.03 in %d of %d, negative in %d of %d" % (np.median(a), *np.percentile(a, [16, 84]), (np.abs(a) > 0.03).sum(), a.size, (a < 0).sum(), a.size)
print("R0 across the PREDICTED TOP (50 m above vs 50 m below):   " + q(rt))
print("R0 across the PREDICTED BASE = top Rogaland Gp:           " + q(rb))
print("largest 50 m step within 100 m of the predicted top:      " + q(rbest))
near = [r for r in out if r["kn"] <= 10 and "r_top" in r]
print("wells within 10 km of a report site (n=%d): R0 at predicted top: %s" % (len(near), " ".join("%s %+.3f" % (r["well"], r["r_top"]) for r in near)))
io = [(r["inside"], r["above"]) for r in out if "inside" in r]
d_ai = np.array([(i[0] - a[0]) / a[0] for i, a in io]); d_vp = np.array([(i[1] - a[1]) / a[1] for i, a in io]); d_rho = np.array([(i[2] - a[2]) / a[2] for i, a in io])
print("inside the predicted clay vs the same thickness of Hordaland above (n=%d): AI %+.1f %% (16-84: %+.1f..%+.1f), Vp %+.1f %%, density %+.1f %%; AI lower inside in %d wells" % (len(io), 100 * np.median(d_ai), *(100 * np.percentile(d_ai, [16, 84])), 100 * np.median(d_vp), 100 * np.median(d_rho), (d_ai < 0).sum()))
print("typical inside values: Vp median %.0f m/s (16-84 %% %.0f-%.0f), density %.2f, AI %.0f" % (np.median([i[1] for i, a in io]), *np.percentile([i[1] for i, a in io], [16, 84]), np.median([i[2] for i, a in io]), np.median([i[0] for i, a in io])))
