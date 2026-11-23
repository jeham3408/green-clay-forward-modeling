"""Scratch: is the predicted clay interval softer than the depth trend of the Hordaland mudstone above it? (open FORCE wells, read-only)
For each well a straight line is fitted to 10 m medians of Vp, density and impedance in the Hordaland Group above the predicted top,
and carried down to the middle of the predicted clay.  Rough: a straight line carried 100-200 m beyond its data."""
import sys, json, glob, numpy as np
R = sys.argv[1]
ex = json.load(open(R + "/natures-plug/data/greenclay.json"))["ex"]; idx = {n: k for k, n in enumerate(ex["n"])}
def read_las(p):
    lines = open(p, errors="replace").read().split("\n")
    c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
    end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
    names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
    well = next(l for l in lines if l.startswith("WELL.")).split(".", 1)[1].split(":")[0].strip()
    arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
    return well, {n: arr[:, k] for k, n in enumerate(names)}
rows = []
for p in sorted(glob.glob(R + "/force2020/las/*/*.las")):
    well, c = read_las(p)
    if well not in idx or "DTC" not in c or "RHOB" not in c: continue
    k = idx[well]; top, th, ho, ro, kn = ex["top"][k], ex["th"][k], ex["ho"][k], ex["ro"][k], ex["kn"][k]
    if top is None or ho is None: continue
    md, dtc, rho = c["DEPT"], c["DTC"], c["RHOB"]
    ok = np.isfinite(dtc) & np.isfinite(rho) & (dtc > 40) & (dtc < 250) & (rho > 1.5) & (rho < 3.0)
    if "CALI" in c and "BS" in c:
        wash = c["CALI"] - c["BS"]; ok &= ~(np.isfinite(wash) & (wash > 3.0))
    vp = 304800.0 / dtc
    def binned(a, b):
        zs, v, r = [], [], []
        for z in np.arange(a, b - 10 + 0.1, 10.0):
            m = ok & (md >= z) & (md < z + 10)
            if m.sum() >= 30: zs.append(z + 5); v.append(np.median(vp[m])); r.append(np.median(rho[m]))
        return np.array(zs), np.array(v), np.array(r)
    za, va, ra = binned(max(ho, top - 400), top)
    zb, vb, rb = binned(top, top + th)
    if za.size < 12 or zb.size < 0.6 * th / 10: continue
    zm = top + th / 2
    pv = np.polyfit(za, va, 1); pr = np.polyfit(za, ra, 1)
    v_tr, r_tr = np.polyval(pv, zm), np.polyval(pr, zm)
    rows.append(dict(well=well, kn=kn, grad=pv[0], dv=np.median(vb) / v_tr - 1, dr=np.median(rb) / r_tr - 1, dai=np.median(vb) * np.median(rb) / (v_tr * r_tr) - 1,
                     vin=np.median(vb), vab=np.median(va), gap=zm - za.mean(), slope_in=np.polyfit(zb, vb, 1)[0]))
a = lambda k, sel=None: np.array([r[k] for r in rows if sel is None or sel(r)])
print("wells: %d (of these within 10 km of a report site: %d)" % (len(rows), sum(1 for r in rows if r["kn"] <= 10)))
for lab, sel in (("all", None), ("within 10 km of a report site", lambda r: r["kn"] <= 10), ("farther away", lambda r: r["kn"] > 10)):
    dv, dr, dai = a("dv", sel), a("dr", sel), a("dai", sel)
    print("%-30s n=%2d | Vp vs trend: median %+.1f %% (16-84 %%: %+.1f..%+.1f), lower in %d | density: %+.1f %% | impedance: %+.1f %% (16-84 %%: %+.1f..%+.1f), lower in %d"
          % (lab, dv.size, 100 * np.median(dv), *(100 * np.percentile(dv, [16, 84])), (dv < 0).sum(), 100 * np.median(dr), 100 * np.median(dai), *(100 * np.percentile(dai, [16, 84])), (dai < 0).sum()))
print("velocity gradient in the Hordaland above the predicted top: median %.2f m/s per m (16-84 %%: %.2f..%.2f); inside the predicted clay: median %.2f m/s per m" % (np.median(a("grad")), *np.percentile(a("grad"), [16, 84]), np.median(a("slope_in"))))
print("distance the line is carried: median %.0f m; Vp above %.0f m/s, inside %.0f m/s (medians of wells)" % (np.median(a("gap")), np.median(a("vab")), np.median(a("vin"))))
