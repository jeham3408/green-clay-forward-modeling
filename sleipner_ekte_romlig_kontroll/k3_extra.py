#!/usr/bin/env python3
"""Independent check, step 3: (B) far/near after taking out a scaled, delayed copy of the data 110/220 ms above
(per stack lags), zone and control windows; (C) fine time profile for the onset of the deficit; (D) datum check at the
well block and the time of the soft clay from the time-depth curve; (E) copy details in control windows.
Uses the functions of k2_analyse.py (everything above its main part). Adds its numbers to results.json under 'extra'."""
import json
import numpy as np

OUTD = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig_kontroll"
src = open(OUTD + "/k2_analyse.py").read().split('RES = {"setup"')[0]
exec(src)

LAG = {"nea": (55, 110), "mid": (56, 112), "far": (60, 120)}      # samples: first and second copy


def resid_energy(X, w, lags):
    """Energy of X in window w per block, raw and after least-squares removal of X(t-L) for L in lags."""
    s = sl(w)
    y = X[:, s]
    A = np.stack([X[:, s.start - L:s.stop - L] for L in lags], 2)            # nb, nt, nl
    G = np.einsum("btl,btm->blm", A, A)
    r = np.einsum("btl,bt->bl", A, y)
    c = np.linalg.solve(G + 1e-12 * np.eye(len(lags))[None] * np.trace(G, axis1=1, axis2=2)[:, None, None], r[..., None])[..., 0]
    raw = (y * y).sum(1)
    return raw, raw - np.einsum("bl,bl->b", c, r), c


EX = {}
for proc in ("p07", "p10"):
    Fl, valid, sh, NBI, NBX, wb = build(proc, 2, 3)
    ok = valid
    iv = np.where(ok)[0]
    groups = superblocks(NBI, NBX, ok)
    ex = {}
    # ---- B
    wins = {"zone_1400_1500": (1400, 1500), "copywin_1366_1466": (1366, 1466), "c_1580_1680": (1580, 1680), "c_1680_1780": (1680, 1780), "c_1780_1880": (1780, 1880)}
    for nl, tag in ((1, "one_lag"), (2, "two_lags")):
        Er, Es, Cf = {}, {}, {}
        for k in ("nea", "mid", "far"):
            for nm, w in wins.items():
                Er[k, nm], Es[k, nm], Cf[k, nm] = resid_energy(Fl[k][0], w, LAG[k][:nl])
        refs = ("c_1580_1680", "c_1680_1780", "c_1780_1880")
        o = {}
        for nm in wins:
            rr = [r for r in refs if r != nm]

            def D(idx, a, b, top, bot, nm=nm, rr=rr):
                num = top[a, nm][idx].sum() / top[b, nm][idx].sum()
                den = sum(bot[a, r][idx].sum() for r in rr) / sum(bot[b, r][idx].sum() for r in rr)
                return float(np.sqrt(num / den)) if num > 0 and den > 0 else np.nan
            d = {}
            for a, b in (("far", "nea"), ("mid", "nea")):
                d[f"{a}_{b}_raw"] = dict(value=r4(D(iv, a, b, Er, Er)), ci95=boot(lambda i: D(i, a, b, Er, Er), groups, 500))
                d[f"{a}_{b}_copy_removed_in_window_only"] = dict(value=r4(D(iv, a, b, Es, Er)), ci95=boot(lambda i: D(i, a, b, Es, Er), groups, 500))
                d[f"{a}_{b}_copy_removed_in_window_and_reference"] = dict(value=r4(D(iv, a, b, Es, Es)), ci95=boot(lambda i: D(i, a, b, Es, Es), groups, 500))
            for k in ("nea", "mid", "far"):
                d[f"energy_fraction_removed_{k}"] = r4(1 - Es[k, nm][iv].sum() / Er[k, nm][iv].sum())
                d[f"median_coef_first_lag_{k}"] = r4(np.median(Cf[k, nm][iv, 0]))
            o[nm] = d
        ex["copy_removed_" + tag] = o
    # ---- C fine profile, 30 ms windows every 6 ms, normalised by 1150-1360 + 1580-1880
    Rw = [(1150, 1360), (1580, 1880)]
    ER = {k: energy(Fl[k], Rw, "rms")[iv].sum() for k in ("nea", "mid", "far")}
    nref = sum(si(w[1]) - si(w[0]) + 1 for w in Rw)
    prof = []
    for tc in range(1330, 1622, 6):
        w = [(tc - 14, tc + 14)]
        e = {k: energy(Fl[k], w, "rms")[iv].sum() for k in ("nea", "mid", "far")}
        prof.append(dict(t=tc, far_near=r4(np.sqrt(e["far"] / e["nea"] / (ER["far"] / ER["nea"]))), mid_near=r4(np.sqrt(e["mid"] / e["nea"] / (ER["mid"] / ER["nea"]))),
                         **{f"level_{k}": r4(np.sqrt(e[k] / 15.0 / (ER[k] / nref))) for k in ("nea", "mid", "far")}))
    ex["fine_profile_30ms"] = prof
    ex["fine_profile_first_centre_far_near_below_0p8"] = next(p["t"] for p in prof if p["far_near"] < 0.8)
    ex["fine_profile_first_centre_far_near_below_0p7"] = next(p["t"] for p in prof if p["far_near"] < 0.7)
    # ---- D datum: extrema of the well-block near and full traces 1250-1420 ms
    for k in ("nea", "ful"):
        tr = Fl[k][0][wb]
        j0, j1 = si(1250), si(1420)
        seg = tr[j0:j1]
        ext = [(T0 + DT * (j0 + i), r4(seg[i] / np.abs(seg).max())) for i in range(1, len(seg) - 1)
               if (seg[i] - seg[i - 1]) * (seg[i + 1] - seg[i]) < 0 and abs(seg[i]) > 0.3 * np.abs(seg).max()]
        ex[f"well_block_{k}_extrema_1250_1420"] = ext
    mn = Fl["nea"][0][iv].mean(0)
    j0, j1 = si(1250), si(1500)
    seg = mn[j0:j1]
    ex["area_mean_near_extrema_1250_1500"] = [(T0 + DT * (j0 + i), r4(seg[i] / np.abs(seg).max())) for i in range(1, len(seg) - 1)
                                              if (seg[i] - seg[i - 1]) * (seg[i + 1] - seg[i]) < 0 and abs(seg[i]) > 0.25 * np.abs(seg).max()]
    # ---- E copy in control windows, per stack (global coefficient against lag)
    cc = {}
    for nm, srcw in (("zone_source_1256_1356", (1256, 1356)), ("source_1436_1536", (1436, 1536)), ("source_1536_1636", (1536, 1636)), ("source_1636_1736", (1636, 1736))):
        row = {}
        for k in ("nea", "mid", "far"):
            Ls, num, den, tt = copy_fit(Fl[k][0], srcw, (90, 134))
            bg = num[iv].sum(0) / den[iv].sum()
            rg = num[iv].sum(0) / np.sqrt(den[iv].sum() * tt[iv].sum(0))
            j = int(np.argmin(rg))
            row[k] = dict(coef=r4(bg[j]), corr=r4(rg[j]), lag_ms=float(Ls[j]))
        cc[nm] = row
    ex["copy_by_source_window"] = cc
    EX[proc] = ex
    print(proc, "done", flush=True)

# time of the soft clay (1383-1415 m MD) from the time-depth curve of 15/9-13
try:
    W = json.load(open("/mnt/ssd/Data/force-hackathon-2026/sleipner/cubes/sleipner94/wells.json"))["wells"][0]
    md = np.array(W["td"]["md"], float)
    tw = np.array(W["td"]["twt"], float)
    m = np.isfinite(md) & np.isfinite(tw)
    EX["soft_clay_twt_ms_from_td_curve"] = dict(md_m=[1383, 1415], twt_ms=[r4(np.interp(1383, md[m], tw[m])), r4(np.interp(1415, md[m], tw[m]))], well=W.get("name"), td_md_range=[float(md[m].min()), float(md[m].max())])
except Exception as e:
    EX["soft_clay_twt_ms_from_td_curve"] = "not available: " + repr(e)
R = json.load(open(OUTD + "/results.json"))
R["extra"] = EX
json.dump(R, open(OUTD + "/results.json", "w"), indent=1)
for p in ("p07", "p10"):
    e = EX[p]
    for tag in ("copy_removed_one_lag", "copy_removed_two_lags"):
        for nm, d in e[tag].items():
            print(p, tag, nm, json.dumps(d))
    print(p, "FP", [(q["t"], q["far_near"], q["mid_near"], q["level_nea"], q["level_mid"], q["level_far"]) for q in e["fine_profile_30ms"]])
    print(p, e["fine_profile_first_centre_far_near_below_0p8"], e["fine_profile_first_centre_far_near_below_0p7"])
    print(p, "well nea", e["well_block_nea_extrema_1250_1420"]); print(p, "well ful", e["well_block_ful_extrema_1250_1420"]); print(p, "mean", e["area_mean_near_extrema_1250_1500"])
    print(p, "CC", json.dumps(e["copy_by_source_window"]))
print(EX["soft_clay_twt_ms_from_td_curve"])
