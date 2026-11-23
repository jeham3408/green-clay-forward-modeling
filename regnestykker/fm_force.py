"""Scratch: shear-log coverage of the predicted Green Clay interval in the open FORCE 2020 wells (read-only)."""
import sys, os, json, glob, numpy as np
R = sys.argv[1]
g = json.load(open(R + "/natures-plug/data/greenclay.json"))["ex"]
ex = {n: k for k, n in enumerate(g["n"])}
def read_las(p):
    lines = open(p, errors="replace").read().split("\n")
    c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
    end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
    names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
    well = next(l for l in lines if l.startswith("WELL.")).split(".", 1)[1].split(":")[0].strip()
    arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
    return well, {n: arr[:, k] for k, n in enumerate(names)}
rows = []; n_dts_any = 0; n_in_model = 0; n_dts_h = 0
for p in sorted(glob.glob(R + "/force2020/las/*/*.las")):
    well, c = read_las(p)
    if well not in ex: continue
    k = ex[well]; top, th, ho, ro, kn = g["top"][k], g["th"][k], g["ho"][k], g["ro"][k], g["kn"][k]
    if top is None or ho is None or ro is None: continue
    n_in_model += 1
    md = c["DEPT"]; dtc = c.get("DTC"); dts = c.get("DTS"); rho = c.get("RHOB"); gr = c.get("GR"); cal = c.get("CALI"); bs = c.get("BS")
    if dts is None or dtc is None or not np.isfinite(dts).any(): continue
    n_dts_any += 1
    ok = np.isfinite(dts) & np.isfinite(dtc) & (dtc > 40) & (dtc < 250) & (dts > 80) & (dts < 800) & (dts / dtc > 1.4) & (dts / dtc < 5)
    inH = (md >= ho) & (md < ro)
    if (ok & inH).sum() > 30: n_dts_h += 1
    A = ok & (md >= ho) & (md < top); B = ok & (md >= top) & (md < top + th)
    step = float(np.nanmedian(np.diff(md)))
    la, lb = A.sum() * step, B.sum() * step
    if lb < 30: 
        rows.append((well, top, th, kn, la, lb, None)); continue
    def med(m):
        vp, vs = 304800 / dtc[m], 304800 / dts[m]; r = rho[m] if rho is not None else np.full(m.sum(), np.nan)
        res = np.polyfit(dtc[m], dts[m], 1); scat = float(np.std(dts[m] - np.polyval(res, dtc[m])))
        return dict(vp=float(np.median(vp)), vs=float(np.median(vs)), vv=float(np.median(vp / vs)), rho=float(np.nanmedian(r)) if np.isfinite(r).any() else float("nan"),
                    gr=float(np.nanmedian(gr[m])) if gr is not None else float("nan"), scat=scat,
                    wash=float(np.nanmedian((cal[m] - bs[m]))) if cal is not None and bs is not None and np.isfinite(cal[m] - bs[m]).any() else float("nan"))
    rows.append((well, top, th, kn, la, lb, (med(A) if la >= 30 else None, med(B))))
print("FORCE wells in the Green Clay model: %d; with a shear log anywhere: %d; with >30 samples of shear in the Hordaland Group: %d" % (n_in_model, n_dts_any, n_dts_h))
print("%-11s %5s %4s %5s | shear m above/in | above: Vp   Vs  Vp/Vs rho  GR  | in clay: Vp   Vs  Vp/Vs rho  GR  wash  DTS-scatter | R0     dVp/Vs" % ("well", "top", "th", "km"))
full = []
for well, top, th, kn, la, lb, m in rows:
    if m is None: print("%-11s %5.0f %4.0f %5.1f | %5.0f /%5.0f | too little shear in the predicted interval" % (well, top, th, kn, la, lb)); continue
    a, b = m
    if a is None:
        print("%-11s %5.0f %4.0f %5.1f | %5.0f /%5.0f | (no shear above)                 | %4.0f %4.0f %5.2f %4.2f %3.0f %5.1f %6.1f" % (well, top, th, kn, la, lb, b["vp"], b["vs"], b["vv"], b["rho"], b["gr"], b["wash"], b["scat"])); continue
    r0 = (b["vp"] * b["rho"] - a["vp"] * a["rho"]) / (b["vp"] * b["rho"] + a["vp"] * a["rho"]) if np.isfinite(a["rho"]) and np.isfinite(b["rho"]) else float("nan")
    print("%-11s %5.0f %4.0f %5.1f | %5.0f /%5.0f | %4.0f %4.0f %5.2f %4.2f %3.0f | %4.0f %4.0f %5.2f %4.2f %3.0f %5.1f %6.1f | %+.3f %+.2f" % (well, top, th, kn, la, lb, a["vp"], a["vs"], a["vv"], a["rho"], a["gr"], b["vp"], b["vs"], b["vv"], b["rho"], b["gr"], b["wash"], b["scat"], r0, b["vv"] - a["vv"]))
    full.append((well, kn, a, b, r0))
if full:
    dv = np.array([b["vv"] - a["vv"] for _, _, a, b, _ in full]); r0 = np.array([r for *_, r in full]); vv = np.array([b["vv"] for _, _, a, b, _ in full])
    print("\nwells with shear both above and inside the predicted interval: %d" % len(full))
    print("Vp/Vs inside: median %.2f (min %.2f, max %.2f); difference inside minus above: median %+.2f, positive in %d of %d" % (np.median(vv), vv.min(), vv.max(), np.median(dv), (dv > 0).sum(), len(dv)))
    print("normal-incidence reflection coefficient, Hordaland above -> predicted clay (medians): median %+.3f, range %+.3f..%+.3f, negative in %d of %d" % (np.nanmedian(r0), np.nanmin(r0), np.nanmax(r0), (r0 < 0).sum(), np.isfinite(r0).sum()))
    near = [x for x in full if x[1] <= 10]
    print("of these, within 10 km of a report that names the Green Clay: %d (%s)" % (len(near), ", ".join(x[0] for x in near)))
