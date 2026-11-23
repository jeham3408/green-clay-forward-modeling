"""Etterregning step 1-3: own stacks from the raw FD gathers, zone/reference with three own estimators, copies,
NMO flatness, mute / fold check. Writes e2_results.json, e2_stacks.npz, prints a summary."""
import sys

import numpy as np

import e1_lib as L

m = L.load_model()
CASES = ("trend", "m07", "m15", "m20")
water = {c: L.load_gather(c, "water") for c in ("G", "FS")}
off = water["G"][2]
t = water["G"][1]
kin = L.Kin(m, off, smooth=11)
R = {"model": {}}
v = m["vp"]
R["model"] = dict(vp_max_above_1500m=float(v[:1500].max()), z_of_that=float(m["z"][:1500][np.argmax(v[:1500])]),
                  n_cells_vp_gt_2700_above_1900m=int((v[:1900] > 2700).sum()),
                  twt_own_minus_model_ms_at_1358=float(1000 * (kin.t0[1357] - m["twt"][1358])),
                  twt_soft_own_ms=[float(1000 * kin.t0[1357]), float(1000 * kin.t0[1389])])
print("MODEL", R["model"], flush=True)

MUTES = {"a48": ("angle", 48.0), "a40": ("angle", 40.0), "s50": ("stretch", 0.5)}
WAVS = {"w25": 25.0, "w35": 35.0}
ST = {}
GATH = {}
for cfg in ("G", "FS"):
    for case in CASES:
        d, tt, o = L.load_gather(cfg, case)
        d = d - water[cfg][0]
        for wn, fw in WAVS.items():
            dd = L.to_wavelet(d, 35.0, fw, point=True)
            for mn, mu in MUTES.items():
                if wn == "w35" and mn != "a48":
                    continue
                c, live = L.nmo(dd, tt, kin, mu)
                ST[(cfg, case, wn, mn)] = L.stack3(c, live, off)
                if case == "trend" and wn == "w25" and mn == "a48":
                    GATH[cfg] = (c, live, dd)
print("stacks done", flush=True)
np.savez_compressed(L.OUT + "/e2_stacks.npz", t=L.T0, **{"%s_%s_%s_%s__%s" % (k + (kk,)): vv for k, S in ST.items() for kk, vv in S.items()})

# ---------------------------------------------------------------- 1. zone / reference
P0 = ST[("G", "trend", "w25", "a48")]
pk_true = L.picks(P0["near"])                         # reflectors picked on the primaries-only near stack
tab = {}
for key, S in ST.items():
    cfg, case, wn, mn = key
    name = "_".join(key)
    row = dict(rms_fn=L.est_rms(S, "far", "near"), rms_mn=L.est_rms(S, "mid", "near"))
    plists = {"true": pk_true}
    if cfg == "FS":
        plists["own"] = L.picks(S["near"])           # reflectors as one would pick them on the contaminated near stack
        full = 0.20 * S["near"] + 0.32 * S["mid"] + 0.45 * S["far"]
        plists["full"] = L.picks(full)
    for pn, pk in plists.items():
        rows = L.per_reflector(S, pk, "far", "near")
        rowsm = L.per_reflector(S, pk, "mid", "near", hw_num=8)
        row["fn_pk_" + pn] = L.zone_ref(rows, "pk")
        row["fn_env_" + pn] = L.zone_ref(rows, "env")
        row["mn_pk_" + pn] = L.zone_ref(rowsm, "pk")
        row["mn_env_" + pn] = L.zone_ref(rowsm, "env")
        if wn == "w25" and mn == "a48" and case in ("trend", "m15"):
            row["rows_" + pn] = rows
    # absolute RMS levels in the zone and reference windows
    row["rms_abs"] = {k: dict(zone=L.rms(S[k], L.ZW), ref=L.rms(S[k], L.RW)) for k in ("near", "mid", "far")}
    tab[name] = row
R["zone_ref"] = tab
R["picks_true_ms"] = [round(1000 * L.T0[i], 1) for i in pk_true]

print("\n=== zone/reference far/near (and mid/near) ===")
print("%-22s %7s %7s %7s | %7s %7s %7s | n(z/r)" % ("case", "RMS", "peak", "env", "mnRMS", "mnpeak", "mnenv"))
for name, row in tab.items():
    for pn in ("true", "own", "full"):
        if "fn_pk_" + pn not in row:
            continue
        f = lambda x: "%7.2f" % x if x is not None else "    nan"
        print("%-22s %7.2f %s %s | %7.2f %s %s | %s %s  zone pk=%s" % (name + ":" + pn, row["rms_fn"], f(row["fn_pk_" + pn]["ratio"]), f(row["fn_env_" + pn]["ratio"]),
              row["rms_mn"], f(row["mn_pk_" + pn]["ratio"]), f(row["mn_env_" + pn]["ratio"]), row["fn_pk_" + pn]["n"], "", row["fn_pk_" + pn]["zone"]))

# decomposition: primaries, multiples only (FS - G), total; trend, 25 Hz, 48 deg
G, F = ST[("G", "trend", "w25", "a48")], ST[("FS", "trend", "w25", "a48")]
dec = {}
for k in ("near", "mid", "far"):
    M = F[k] - G[k]
    dec[k] = dict(prim_zone=L.rms(G[k], L.ZW), prim_ref=L.rms(G[k], L.RW), mult_zone=L.rms(M, L.ZW), mult_ref=L.rms(M, L.RW),
                  tot_zone=L.rms(F[k], L.ZW), tot_ref=L.rms(F[k], L.RW),
                  corr_mult_prim_zone=float(np.corrcoef(M[L.ZW], G[k][L.ZW])[0, 1]), corr_mult_prim_ref=float(np.corrcoef(M[L.RW], G[k][L.RW])[0, 1]))
    dec[k]["mult_over_prim_zone"] = dec[k]["mult_zone"] / dec[k]["prim_zone"]
    dec[k]["mult_over_prim_ref"] = dec[k]["mult_ref"] / dec[k]["prim_ref"]
    dec[k]["tot_over_prim_zone"] = dec[k]["tot_zone"] / dec[k]["prim_zone"]
    dec[k]["tot_over_prim_ref"] = dec[k]["tot_ref"] / dec[k]["prim_ref"]
R["decomposition_trend"] = dec
print("\n=== decomposition (trend): RMS multiples-only / primaries, total / primaries ===")
for k, dd in dec.items():
    print(k, {a: round(b, 3) for a, b in dd.items() if "over" in a or "corr" in a})

# effect of Vs on the primaries-only far stack inside 1380-1500 ms
w = L.win(1.380, 1.500)
vs_eff = {}
for case in ("m07", "m15", "m20"):
    A, B = ST[("G", "trend", "w25", "a48")], ST[("G", case, "w25", "a48")]
    vs_eff[case] = {k: dict(rms_change_rel=L.rms(B[k] - A[k], w) / L.rms(A[k], w), rms_level_ratio=L.rms(B[k], w) / L.rms(A[k], w)) for k in ("near", "mid", "far")}
R["vs_effect_G_1380_1500"] = vs_eff
print("\n=== Vs effect on primaries-only stacks, 1380-1500 ms ===")
for case, dd in vs_eff.items():
    print(case, {k: (round(x["rms_change_rel"], 3), round(x["rms_level_ratio"], 3)) for k, x in dd.items()})

# ---------------------------------------------------------------- 2. copies
cop = {}
wsel = L.win(1.10, 1.36)
order = sorted([i for i in pk_true if 1.10 <= L.T0[i] <= 1.36], key=lambda i: -abs(P0["near"][i]))[:3]
for ic in order:
    key = "parent_%dms" % round(1000 * L.T0[ic])
    cop[key] = {}
    for k in ("near", "mid", "far"):
        Mk = F[k] - G[k]
        cop[key][k] = dict(
            mult_only_first=L.copy_fit(G[k], Mk, ic, 0.085, 0.175, -1), mult_only_second=L.copy_fit(G[k], Mk, ic, 0.185, 0.340, +1),
            in_FS_first=L.copy_fit(F[k], F[k], ic, 0.085, 0.175, -1), in_FS_second=L.copy_fit(F[k], F[k], ic, 0.185, 0.340, +1),
            in_G_first=L.copy_fit(G[k], G[k], ic, 0.085, 0.175, -1))
R["copies"] = cop
print("\n=== copies below the strongest reflectors (trend) ===")
for key, dd in cop.items():
    for k, x in dd.items():
        print(key, k, "mult-only:", x["mult_only_first"], x["mult_only_second"], "| inside FS:", x["in_FS_first"], x["in_FS_second"], "| inside G:", x["in_G_first"])

# theory for the first and second peg-leg of those parents: lag after NMO and amplitude per offset, stacked
import e3_zoep as Z
sb = dict(vp=float(np.mean(m["vp"][81:85])), vs=float(np.mean(m["vs"][81:85])), rho=float(np.mean(m["rho"][81:85])))
R["seabed"] = dict(sb, R0=(sb["vp"] * sb["rho"] - 1480 * 1.03) / (sb["vp"] * sb["rho"] + 1480 * 1.03))
th = {}
tabs = {nw: L.raytab(kin.v, off, n_water=nw) for nw in (1, 2)}
for ic in order:
    key = "parent_%dms" % round(1000 * L.T0[ic])
    k0 = int(np.argmin(np.abs(kin.t0 - L.T0[ic])))
    th[key] = {}
    for nw in (1, 2):
        txm, pxm, dxm = tabs[nw]
        lag = np.full(len(off), np.nan)
        amp = np.full(len(off), np.nan)
        for ix in range(len(off)):
            tm = txm[k0, ix]
            col = kin.tx[:, ix]
            ok = np.isfinite(col)
            if not np.isfinite(tm) or ok.sum() < 10 or not np.isfinite(kin.tx[k0, ix]):
                continue
            t0eq = np.interp(tm, col[ok], kin.t0[ok])
            lag[ix] = t0eq - kin.t0[k0]
            pm, pp = pxm[k0, ix], kin.px[k0, ix]
            rsb = Z.fluid_solid_rpp(1480.0, 1.03, sb["vp"], sb["vs"], sb["rho"], pm)
            cm, cp = np.sqrt(1 - (pm * 1480.0) ** 2), np.sqrt(1 - (pp * 1480.0) ** 2)
            spread = np.sqrt((kin.dxdp[k0, ix] * cp * cp) / (dxm[k0, ix] * cm * cm))
            amp[ix] = (nw + 1) * (-rsb) ** nw * spread
        for sk, (a, b) in L.RANGES.items():
            sel = (off >= a - 0.5) & (off <= b + 0.5) & GATH["G"][1][:, ic] & np.isfinite(lag)
            if sel.sum() == 0:
                continue
            # stack the predicted copies (parent waveform from the primaries stack, shifted by each offset's lag)
            par = G[sk]
            pred = np.zeros(len(L.T0))
            for ix in np.where(sel)[0]:
                pred += amp[ix] * np.interp(L.T0 - lag[ix], L.T0, par, left=0, right=0)
            pred /= sel.sum()
            only = np.zeros(len(L.T0))
            only[ic - 20:ic + 41] = par[ic - 20:ic + 41]
            pred1 = np.zeros(len(L.T0))
            for ix in np.where(sel)[0]:
                pred1 += amp[ix] * np.interp(L.T0 - lag[ix], L.T0, only, left=0, right=0)
            pred1 /= sel.sum()
            lo, hi, sg = ((0.085, 0.175, -1) if nw == 1 else (0.185, 0.340, +1))
            th[key].setdefault(sk, {})["order%d" % nw] = dict(
                lag_ms_min_mean_max=[round(1000 * float(np.nanmin(lag[sel])), 1), round(1000 * float(np.nanmean(lag[sel])), 1), round(1000 * float(np.nanmax(lag[sel])), 1)],
                amp_per_offset_min_mean_max=[round(float(np.nanmin(amp[sel])), 3), round(float(np.nanmean(amp[sel])), 3), round(float(np.nanmax(amp[sel])), 3)],
                stacked_isolated_parent=L.copy_fit(only, pred1, ic, lo, hi, sg), n_offsets=int(sel.sum()))
R["copies_theory"] = th
print("\n=== theory for the copies (own rays + own fluid-solid Zoeppritz) ===")
for key, dd in th.items():
    for sk, x in dd.items():
        for o_, y in x.items():
            print(key, sk, o_, y)

# ---------------------------------------------------------------- 3. NMO flatness, fold, time zero
fl = {}
for cfg in ("G", "FS"):
    c, live, raw = GATH[cfg]
    S = ST[(cfg, "trend", "w25", "a48")]
    for ic in order[:2]:
        ref = S["near"][ic - 25:ic + 26]
        res = {}
        for x in (150, 300, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000):
            ix = int(np.argmin(np.abs(off - x)))
            if not live[ix, ic - 40:ic + 41].all():
                continue
            lags = np.arange(-15, 16)
            cc = np.array([np.sum(ref * c[ix, ic - 25 + l:ic + 26 + l]) for l in lags])
            j = int(np.clip(np.argmax(cc), 1, len(lags) - 2))
            den = cc[j - 1] - 2 * cc[j] + cc[j + 1]
            res[str(x)] = round(float(lags[j] + (0.5 * (cc[j - 1] - cc[j + 1]) / den if den != 0 else 0.0)), 2)
        fl["%s_event_%dms" % (cfg, round(1000 * L.T0[ic]))] = res
R["nmo_residual_ms_vs_near_stack"] = fl
print("\n=== NMO residual (ms) per offset against the near stack ===")
for k, x in fl.items():
    print(k, x)

fold = {}
for mn in MUTES:
    S = ST[("G", "trend", "w25", mn)]
    fold[mn] = {k: {str(int(1000 * tt_)): int(S["n_" + k][int(round(tt_ / L.DT))]) for tt_ in (1.10, 1.20, 1.30, 1.38, 1.408, 1.45, 1.492, 1.56, 1.70, 1.90)} for k in ("near", "mid", "far")}
    nf = S["n_far"]
    fold[mn]["far_fold_monotonic_1100_1900"] = bool(np.all(np.diff(nf[1100:1901]) >= 0))
R["fold"] = fold
print("\n=== fold (live traces) ===")
for mn, x in fold.items():
    print(mn, "far:", x["far"], "monotonic:", x["far_fold_monotonic_1100_1900"])
    print(mn, "mid:", x["mid"])

# sea-bed reflection at the nearest offset (150 m): expected centre = surface-to-surface ray time sqrt(x^2 + (2h)^2)/c
tz = {}
for cfg in ("G", "FS"):
    d, tt, o = L.load_gather(cfg, "trend")
    d = L.to_wavelet(d - water[cfg][0], 35.0, 35.0, point=True)[0]
    te = np.sqrt(150.0 ** 2 + 162.0 ** 2) / 1480.0
    wv = (tt > te - 0.03) & (tt < te + 0.03)
    env = np.abs(__import__("scipy.signal").signal.hilbert(d))
    ipk = np.where(wv)[0][np.argmax(np.abs(d[wv]))]
    ie = np.where(wv)[0][np.argmax(env[wv])]
    tz[cfg] = dict(expected_ms=round(1000 * te, 2), abs_peak_ms=round(1000 * tt[ipk], 1), envelope_peak_ms=round(1000 * tt[ie], 1), peak_value=float(d[ipk]))
R["seabed_time_150m"] = tz
print("\nsea-bed 150 m:", tz)
print("clipping check: max |d| per gather / float32 max is irrelevant; NaN/Inf:", {"%s_%s" % (c_, k_): bool(np.isfinite(L.load_gather(c_, k_)[0]).all()) for c_ in ("G", "FS") for k_ in CASES + ("water",)})
L.jdump(R, "e2_results.json")
print("saved")
