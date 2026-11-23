"""Point 1, the crux: how much does the result for ONE low-Vs layer depend on where exactly its top and base are put?

The experiment lowered Vs in 1358-1390 m. Here: boxes with the top between 1340 and 1366 m and the base between 1382 and
1410 m (2 m steps), Vs 7 / 15 / 20 % below trend, and the layer bounded by the impedance log itself (cells softer than
the 61 m running mean of the acoustic impedance, inside the zone). Two set-ups: 25 Hz Ricker without NMO stretch, and
the ghosted 25 Hz Ricker with NMO (closest to the experiment). Writes p7_results.json.
"""
import json, os
import numpy as np
import slp, vz
HERE = os.path.dirname(os.path.abspath(__file__))
m = vz.load_model(); z = m["z"]; n = len(z)
F = vz.Forward(m)
ai = m["vp"] * m["rho"]
def runmean(x, L):
    k = np.ones(L) / L
    return np.convolve(np.pad(x, L // 2, mode="edge"), k, mode="valid")[:len(x)]
def ricker(f, nn=401):
    t = (np.arange(nn) - nn // 2) * vz.DT
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)
def ghosted(w, zs=6.0, zr=8.0, c=1480.0):
    ts, tr = 2 * zs / c / vz.DT, 2 * zr / c / vz.DT
    f = np.fft.rfftfreq(len(w)); sh = lambda s: np.exp(-2j * np.pi * f * s)
    return np.fft.irfft(np.fft.rfft(w) * (1 - sh(ts)) * (1 - sh(tr)) * sh(-(ts + tr) / 2), len(w))
zt, zb = [float(np.interp(t, m["twt"], m["ze"])) for t in vz.ZONE]
OUT = dict(zone_depth_m=[zt, zb])
# which cells does "softer than the running mean" pick in the zone?
iv = {}
for L in (31, 61, 101):
    d = runmean(ai, L) - ai
    c = (z >= zt) & (z < zb) & (d > 0)
    idx = np.nonzero(c)[0]
    runs, s = [], None
    for i in idx:
        if s is None: s = p = i
        elif i == p + 1: p = i
        else: runs.append([float(m["ze"][s]), float(m["ze"][p + 1])]); s = p = i
    if s is not None: runs.append([float(m["ze"][s]), float(m["ze"][p + 1])])
    iv["L%d" % L] = dict(intervals_m=runs, total_m=int(c.sum()))
OUT["soft_bed_intervals_in_zone"] = iv
prof = [dict(z=float(a), ai=round(float(ai[(z >= a) & (z < a + 2)].mean())), ai_mean61=round(float(runmean(ai, 61)[(z >= a) & (z < a + 2)].mean()))) for a in np.arange(1340.0, 1412.0, 2.0)]
OUT["ai_profile_2m"] = prof
def box(a, b, X):
    f = np.ones(n); f[(z >= a) & (z < b)] = 1 - X; return f
tops = np.arange(1340.0, 1366.1, 2.0); bases = np.arange(1382.0, 1410.1, 2.0)
RES = {}
for wn, w, st in (("ricker25|nmofree", ricker(25.0), False), ("ricker25_ghosted|nmo", ghosted(ricker(25.0)), True)):
    F.wavelet = lambda f=None, n=301, w=w: w
    S0 = F.stacks(F.R0, stretch=st)
    ev0 = slp.pick_events(S0["ful"])
    base = dict(rms=vz.rms_measure(S0)["fn"], ev=vz.event_measure(S0, ev0)["fn"])
    row = dict(baseline=base)
    for X in (0.07, 0.15, 0.20):
        g_r = np.zeros((len(tops), len(bases))); g_e = np.zeros_like(g_r)
        for i, a in enumerate(tops):
            for j, b in enumerate(bases):
                S = F.stacks(F.refl(box(a, b, X)), stretch=st)
                g_r[i, j] = vz.rms_measure(S)["fn"]
                e = vz.event_measure(S, ev0)["fn"]
                g_e[i, j] = np.nan if e is None else e
        ie = (int(np.argmin(np.abs(tops - 1358))), int(np.argmin(np.abs(bases - 1390))))
        im = np.unravel_index(np.nanargmin(g_r), g_r.shape)
        row["%02d" % round(100 * X)] = dict(experiment_box_1358_1390=dict(rms=float(g_r[ie]), ev=float(g_e[ie])),
            min_rms=dict(top=float(tops[im[0]]), base=float(bases[im[1]]), rms=float(g_r[im]), ev=float(g_e[im])),
            rms_min_max=[float(np.nanmin(g_r)), float(np.nanmax(g_r))], ev_min_max=[float(np.nanmin(g_e)), float(np.nanmax(g_e))],
            share_rms_below_0p75=float(np.mean(g_r < 0.75)), share_rms_below_0p65=float(np.mean(g_r < 0.65)),
            rms_grid=[[round(float(v), 2) for v in r] for r in g_r], ev_grid=[[None if not np.isfinite(v) else round(float(v), 2) for v in r] for r in g_e])
    RES[wn] = row
OUT["tops_m"] = [float(v) for v in tops]; OUT["bases_m"] = [float(v) for v in bases]
OUT["boxes"] = RES
json.dump(OUT, open(os.path.join(HERE, "p7_results.json"), "w"), indent=1, default=float)
print("zone depth", OUT["zone_depth_m"]); print("soft-bed intervals in the zone:", json.dumps(iv))
print("AI (2 m) and its 61 m running mean:", " ".join("%d:%d/%d" % (p["z"], p["ai"], p["ai_mean61"]) for p in prof))
for wn, row in RES.items():
    print("\n", wn, "baseline", {k: round(v, 2) for k, v in row["baseline"].items()})
    for X in ("07", "15", "20"):
        r = row[X]
        print("  Vs -%s %%: experiment's box 1358-1390: RMS %.2f ev %.2f | min RMS %.2f (ev %.2f) at %d-%d m | RMS range %.2f-%.2f | ev range %.2f-%.2f | share RMS<0.75: %.2f, <0.65: %.2f" % (
            X, r["experiment_box_1358_1390"]["rms"], r["experiment_box_1358_1390"]["ev"], r["min_rms"]["rms"], r["min_rms"]["ev"], r["min_rms"]["top"], r["min_rms"]["base"],
            r["rms_min_max"][0], r["rms_min_max"][1], r["ev_min_max"][0], r["ev_min_max"][1], r["share_rms_below_0p75"], r["share_rms_below_0p65"]))
    print("  RMS grid at 15 %% (rows: top %s; columns: base %s)" % (OUT["tops_m"], OUT["bases_m"]))
    for t_, r_ in zip(tops, row["15"]["rms_grid"]):
        print("   %5.0f  " % t_ + " ".join("%.2f" % v for v in r_))
