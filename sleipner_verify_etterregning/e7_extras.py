"""Extras on the saved own stacks (e2_stacks.npz):
(a) estimators restricted to samples with far fold >= 5 (needed for the stricter mutes),
(b) how special is the zone? the same far/near RMS ratio in every 84 ms window between 1.1 and 1.9 s,
(c) copies on the small validation model with own estimator and own theory,
(d) which offset does a copy lag of 110 / 112 / 119 ms (real data) correspond to,
(e) own stacks against the experiment's stacks.npz, and the experiment's far fold,
(f) own three-point multiple filter on the stacks (theory coefficients and coefficients measured in the stacks)."""
import numpy as np

import e1_lib as L
import e3_zoep as Z

E = np.load(L.OUT + "/e2_stacks.npz")
T = E["t"]
def S_of(cfg, case, wn, mn):
    return {k: E["%s_%s_%s_%s__%s" % (cfg, case, wn, mn, k)] for k in ("near", "mid", "far", "n_near", "n_mid", "n_far")}
R = {}

# ---- (a)
print("=== (a) zone/reference with far fold >= 5 only ===")
a = {}
for mn in ("a48", "a40", "s50"):
    P0 = S_of("G", "trend", "w25", mn)
    ok = P0["n_far"] >= 5
    rw = L.RW & ok
    R.setdefault("ref_window_samples", {})[mn] = dict(ref_ms_total=int(L.RW.sum()), ref_ms_with_far_fold_ge5=int(rw.sum()), first_time_with_far_fold_ge5_ms=float(1000 * T[np.argmax(ok & (T > 1.0))]))
    pk = [i for i in L.picks(P0["near"]) if ok[i - 12:i + 13].all()]
    for cfg in ("G", "FS"):
        for case in ("trend", "m07", "m15", "m20"):
            S = S_of(cfg, case, "w25", mn)
            rr = (L.rms(S["far"], L.ZW) / L.rms(S["near"], L.ZW)) / (L.rms(S["far"], rw) / L.rms(S["near"], rw))
            rows = L.per_reflector(S, pk, "far", "near")
            zp, ze = L.zone_ref(rows, "pk"), L.zone_ref(rows, "env")
            a["%s_%s_%s" % (cfg, case, mn)] = dict(rms=rr, pk=zp["ratio"], env=ze["ratio"], n=zp["n"], zone_pk=zp["zone"], ref_pk_minmax=zp["ref_minmax"])
            print("%-16s rms %.2f  peak %s  env %s  n %s  zone %s  ref min/max %s" % ("%s_%s_%s" % (cfg, case, mn), rr, "%.2f" % zp["ratio"] if zp["ratio"] else "nan", "%.2f" % ze["ratio"] if ze["ratio"] else "nan", zp["n"], zp["zone"], zp["ref_minmax"]))
R["fold_ge5"] = a

# ---- (b)
print("\n=== (b) every 84 ms window, 1.1-1.9 s: far/near (and mid/near) RMS in the window over the same outside (window -28 ms .. +68 ms excluded) ===")
b = {}
for cfg in ("G", "FS"):
    for case in ("trend", "m15"):
        S = S_of(cfg, case, "w25", "a48")
        res = []
        for c0 in np.arange(1.100, 1.900 - 0.084 + 1e-9, 0.004):
            zw = (T >= c0) & (T <= c0 + 0.084)
            rw = (T >= 1.1) & (T <= 1.9) & ~((T >= c0 - 0.028) & (T <= c0 + 0.084 + 0.068))
            fn = (L.rms(S["far"], zw) / L.rms(S["near"], zw)) / (L.rms(S["far"], rw) / L.rms(S["near"], rw))
            mnr = (L.rms(S["mid"], zw) / L.rms(S["near"], zw)) / (L.rms(S["mid"], rw) / L.rms(S["near"], rw))
            res.append((c0, fn, mnr))
        res = np.array(res)
        iz = int(np.argmin(np.abs(res[:, 0] - 1.408)))
        q = lambda x: [round(float(v), 2) for v in np.percentile(x, [0, 10, 25, 50, 75, 90, 100])]
        b["%s_%s" % (cfg, case)] = dict(fn_percentiles_0_10_25_50_75_90_100=q(res[:, 1]), mn_percentiles=q(res[:, 2]), zone_fn=float(res[iz, 1]), zone_mn=float(res[iz, 2]),
                                        share_of_windows_with_fn_at_or_below_zone=float(np.mean(res[:, 1] <= res[iz, 1])),
                                        share_of_windows_with_fn_below_0p65=float(np.mean(res[:, 1] < 0.65)), share_fn_above_1p5=float(np.mean(res[:, 1] > 1.5)),
                                        windows_ms_fn_below_0p65=[int(round(1000 * x)) for x in res[res[:, 1] < 0.65, 0]], n_windows=int(len(res)))
        print(cfg, case, {k: v for k, v in b["%s_%s" % (cfg, case)].items()})
R["sliding_windows"] = b

# ---- (c)
print("\n=== (c) validation model: copy / primary on raw traces, own estimator, and own theory ===")
c = {}
vcol = np.where(np.arange(800) + 0.5 < 81.0, 1480.0, 1800.0)
offv = np.arange(0.0, 900.1, 25.0)
offv_t = np.maximum(offv, 1.0)
tx0, px0, dx0 = L.raytab(vcol, offv_t)
tx1, px1, dx1 = L.raytab(vcol, offv_t, n_water=1)
for cfg in ("FS", "G"):
    z1, z0 = np.load(L.OUT + "/theirs/gather_%s_val1.npz" % cfg), np.load(L.OUT + "/theirs/gather_%s_val0.npz" % cfg)
    d = z1["d"].astype(float) - z0["d"].astype(float)
    tv, ov = z1["t"], z1["offsets"]
    rows = []
    for wname, dd in (("35Hz_as_modelled", d), ("25Hz_converted", L.to_wavelet(d, 35.0, 25.0, point=True))):
        for x in (0.0, 200.0, 400.0, 600.0, 800.0):
            ix = int(np.argmin(np.abs(ov - x)))
            tp = tx0[-1, ix]
            w = (tv > tp - 0.03) & (tv < tp + 0.05)
            ic = int(np.where(w)[0][np.argmax(np.abs(dd[ix][w]))])
            m_ = L.copy_fit(dd[ix], dd[ix], ic, 0.085, 0.130, -1)
            pm, pp = px1[-1, ix], px0[-1, ix]
            rsb = Z.fluid_solid_rpp(1480.0, 1.03, 1800.0, 600.0, 1.95, pm)
            cm, cp = np.sqrt(1 - (pm * 1480.0) ** 2), np.sqrt(1 - (pp * 1480.0) ** 2)
            th = -2 * rsb * np.sqrt((dx0[-1, ix] * cp * cp) / (dx1[-1, ix] * cm * cm))
            rows.append(dict(wavelet=wname, offset=x, measured=m_["coef"], lag_ms=m_["lag_ms"], ncc=m_["ncc"], theory=round(float(th), 3), theory_lag_ms=round(1000 * float(tx1[-1, ix] - tx0[-1, ix]), 1), R_seabed=round(float(rsb), 4)))
            print(cfg, rows[-1])
    c[cfg] = rows
R["validation_copies"] = c

# ---- (d)
print("\n=== (d) lag of the first / second peg-leg after NMO against offset (parents at 1276 and 1316 ms) ===")
m = L.load_model()
off = np.arange(150.0, 3150.1, 25.0)
kin = L.Kin(m, off, smooth=11)
dd_ = {}
for nw in (1, 2):
    txm = L.raytab(kin.v, off, n_water=nw)[0]
    for tpar in (1.276, 1.316):
        k0 = int(np.argmin(np.abs(kin.t0 - tpar)))
        lag = np.full(len(off), np.nan)
        for ix in range(len(off)):
            col = kin.tx[:, ix]
            ok = np.isfinite(col)
            if np.isfinite(txm[k0, ix]) and ok.sum() > 10:
                lag[ix] = 1000 * (np.interp(txm[k0, ix], col[ok], kin.t0[ok]) - kin.t0[k0])
        g = np.isfinite(lag)
        tabl = {str(int(x)): round(float(lag[int(np.argmin(np.abs(off - x)))]), 1) for x in (150, 450, 760, 1000, 1250, 1500, 1760, 2000, 2500, 3000)}
        inv = {str(v): round(float(np.interp(v, lag[g], off[g])), 0) for v in ((110, 112, 119) if nw == 1 else (220, 224, 234))}
        means = {sk: round(float(np.nanmean(lag[(off >= a - 0.5) & (off <= b + 0.5)])), 1) for sk, (a, b) in L.RANGES.items()}
        dd_["order%d_parent_%d" % (nw, int(1000 * tpar))] = dict(lag_ms_by_offset=tabl, offset_m_for_real_lags=inv, mean_lag_per_offset_range=means)
        print("order", nw, "parent", tpar, tabl, "| offset for the real lags:", inv, "| mean per range:", means)
R["lag_vs_offset"] = dd_

# ---- (e)
print("\n=== (e) own stacks against the experiment's stacks.npz (1.1-1.9 s) ===")
TH = np.load(L.OUT + "/theirs/stacks.npz")
tt = TH["t"]
e = {}
w = (T >= 1.1) & (T <= 1.9)
for cfg in ("G", "FS"):
    for case in ("trend", "m15"):
        S = S_of(cfg, case, "w25", "a48")
        row = {}
        for mine, th in (("near", "nea"), ("mid", "mid"), ("far", "far")):
            y = np.interp(T, tt, TH["%s_%s_main__%s" % (cfg, case, th)])
            row[mine] = dict(corr=round(float(np.corrcoef(S[mine][w], y[w])[0, 1]), 4), rms_ratio_mine_over_theirs=round(float(L.rms(S[mine], w) / np.sqrt(np.mean(y[w] ** 2))), 4))
        e["%s_%s" % (cfg, case)] = row
        print(cfg, case, row)
nf = np.interp(T, tt, TH["G_trend_main__n_far"].astype(float))
dn = np.diff(TH["G_trend_main__n_far"].astype(float)[(tt >= 1.1) & (tt <= 1.9)])
e["their_far_fold"] = dict(at_ms={str(int(1000 * x)): float(nf[int(round(x / L.DT))]) for x in (1.10, 1.20, 1.30, 1.38, 1.408, 1.45, 1.492, 1.56, 1.70, 1.90)},
                           n_decreases_1100_1900=int((dn < 0).sum()), largest_decrease=float(dn.min()), n_increases=int((dn > 0).sum()),
                           zone_min_max=[float(nf[L.ZW].min()), float(nf[L.ZW].max())])
print("their far fold:", e["their_far_fold"])
R["vs_their_stacks"] = e

# ---- (f)
print("\n=== (f) own three-point multiple filter on the with-surface stacks ===")
def three(x, a1, l1, a2, l2):
    y = x.copy()
    n1, n2 = int(round(l1)), int(round(l2))
    y[n1:] -= a1 * x[:-n1]
    y[n2:] += (a1 * a1 - a2) * x[:-n2]
    return y
f = {}
G0 = S_of("G", "trend", "w25", "a48")
pk_true = L.picks(G0["near"])
theory = {"near": (-0.75, 111, 0.43, 221), "mid": (-0.60, 117, 0.20, 232), "far": (-0.62, 128, 0.19, 254)}     # e2 theory, stacked isolated parent (1245 ms)
for case in ("trend", "m15"):
    S = S_of("FS", case, "w25", "a48")
    G = S_of("G", case, "w25", "a48")
    meas = {}
    par = sorted([i for i in L.picks(0.20 * S["near"] + 0.32 * S["mid"] + 0.45 * S["far"]) if 1.10 <= T[i] <= 1.36], key=lambda i: -abs(S["near"][i]))[:2]
    for k in ("near", "mid", "far"):
        c1 = [L.copy_fit(S[k], S[k], i, 0.095, 0.140, -1) for i in par]
        c2 = [L.copy_fit(S[k], S[k], i, 0.200, 0.270, +1) for i in par]
        meas[k] = (float(np.mean([x["coef"] for x in c1])), float(np.mean([x["lag_ms"] for x in c1])), float(np.mean([x["coef"] for x in c2])), float(np.mean([x["lag_ms"] for x in c2])))
    for name, co in (("theory_coefficients", theory), ("coefficients_measured_in_stacks", meas)):
        Fs = {k: three(S[k], *co[k]) for k in ("near", "mid", "far")}
        rr = L.est_rms(Fs)
        rows = L.per_reflector(Fs, pk_true, "far", "near")
        own = L.per_reflector(Fs, L.picks(Fs["near"]), "far", "near")
        f["%s_%s" % (case, name)] = dict(coefficients={k: [round(float(v), 3) for v in co[k]] for k in co}, rms=rr, pk_true=L.zone_ref(rows, "pk")["ratio"], env_true=L.zone_ref(rows, "env")["ratio"],
                                         pk_own=L.zone_ref(own, "pk")["ratio"], env_own=L.zone_ref(own, "env")["ratio"], n_own=L.zone_ref(own, "pk")["n"],
                                         residual_rms_vs_primaries={k: round(L.rms(Fs[k] - G[k], w) / L.rms(G[k], w), 3) for k in Fs},
                                         unfiltered_residual_rms_vs_primaries={k: round(L.rms(S[k] - G[k], w) / L.rms(G[k], w), 3) for k in Fs})
        print(case, name, f["%s_%s" % (case, name)])
R["three_point_filter"] = f
L.jdump(R, "e7_results.json")
print("saved")
