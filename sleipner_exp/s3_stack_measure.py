"""Steps 3-5: direct-wave removal, NMO, near/mid/far stacks, far/near per reflector, copies, multiple filter.

Reads gather_<cfg>_<case>.npz (s2_fd.py) and model_1d.npz (s1_model.py) from SLEIPNER_EXP_ROOT.
Writes stacks.npz, results.json and prints the tables.

Run:  SLEIPNER_EXP_ROOT=... python -B s3_stack_measure.py
"""
import json
import os

import numpy as np

import slp

ROOT = slp.ROOT
F_SRC = 35.0
VARIANTS = {                       # name: (wavelet peak Hz of the source, NMO mode, mute angle, stretch mute)
    "main": (25.0, "exact", 48.0, None),
    "stretch50": (25.0, "exact", None, 0.5),
    "hyper": (25.0, "hyper", 48.0, None),
    "w35": (35.0, "exact", 48.0, None),
    "reallag": (25.0, "exact", 48.0, None),      # as main, but offset ranges that match the multiple lags of the real stacks
}
VRANGE = {"reallag": slp.RANGES_LAG}             # near 150-570, mid 575-1000, far 1200-1800 m (all others: 150-760 / 760-1760 / 1760-3150)
SCALE = dict(nea=0.73, mid=0.62, far=0.55)       # copy strength of the real stacks relative to the model's (checker's p2_scale.py)
CASES = ("trend", "m07", "m15", "m20")
CASES_C = CASES + ("m30", "m40", "sb07", "sb15", "sb20")          # the ray-theory synthetics are cheap: two stronger Vs anomalies as well
CASES_FD = CASES + ("trend_sbl", "trend_sbh", "sb07", "sb15")   # FD only: sea-bed density 1.80 / 2.10 g/cm3 (sea-bed reflection 0.37 / 0.43)


def load(cfg, case, suf=""):
    f = os.path.join(ROOT, "gather_%s_%s%s.npz" % (cfg, case, suf))
    if not os.path.exists(f):
        return None
    z = np.load(f)
    return z["d"].astype(np.float64), z["t"], z["offsets"], json.loads(str(z["info"]))


def main():
    m0 = slp.model("trend")
    off = np.arange(150.0, 3150.1, 25.0)
    rays = slp.Rays(m0, off)
    rays_m = [(slp.Rays(m0, off, n_water=nw), nw) for nw in (1, 2, 3)]      # water-layer multiples for the ray-theory cross-check
    ST, INFO, NMOG, RAW = {}, {}, {}, {}
    water = {cfg: load(cfg, "water") for cfg in ("G", "FS")}
    for vname, (fw, mode, mdeg, smax) in VARIANTS.items():
        for cfg in ("G", "FS"):
            if water[cfg] is None:
                continue
            for case in CASES_FD:
                g = load(cfg, case)
                if g is None:
                    continue
                d, t, o, info = g
                INFO["%s_%s" % (cfg, case)] = info
                d = d - water[cfg][0]
                d = slp.reshape_wavelet(d, F_SRC, fw)
                c, live = slp.nmo(d, t, rays, mode, mdeg, smax)
                ST["%s_%s_%s" % (cfg, case, vname)] = slp.stacks(c, live, o, VRANGE.get(vname))
                if case == "trend" and vname in ("main", "hyper"):
                    NMOG["%s_%s" % (cfg, vname)] = (c, live)
                if case == "trend" and vname == "main":
                    RAW[cfg] = (d, t)
        for case in CASES_C:                              # primaries-only Zoeppritz synthetic with the same ghosts
            if vname not in ("main", "reallag"):
                continue
            mc = slp.model(case)
            d, t = slp.conv_gather(mc, rays, fw, ghosts=True)
            c, live = slp.nmo(d, t, rays, mode, mdeg, smax)
            ST["C_%s_%s" % (case, vname)] = slp.stacks(c, live, off, VRANGE.get(vname))
            d, t = slp.conv_gather(mc, rays, fw, ghosts=True, mult=rays_m)       # the same plus water-layer multiples (ray theory)
            c, live = slp.nmo(d, t, rays, mode, mdeg, smax)
            ST["CM_%s_%s" % (case, vname)] = slp.stacks(c, live, off, VRANGE.get(vname))
    g1, w1 = load("FS", "trend", "_dx1"), load("FS", "water", "_dx1")
    if g1 is not None and w1 is not None:
        INFO["FS_trend_dx1"] = g1[3]
        for vname in ("main", "w35"):
            fw, mode, mdeg, smax = VARIANTS[vname]
            d = slp.reshape_wavelet(g1[0] - w1[0], F_SRC, fw)
            c, live = slp.nmo(d, g1[1], rays, mode, mdeg, smax)
            ST["FSdx1_trend_%s" % vname] = slp.stacks(c, live, g1[2])
    gs = load("FS", "trend", "_sym")                    # the same case with the exactly mirror-symmetric shear modulus (SLEIPNER_SYMSLIP=1)
    if gs is not None and water["FS"] is not None:
        INFO["FS_trend_sym"] = gs[3]
        for vname in ("main", "reallag"):
            fw, mode, mdeg, smax = VARIANTS[vname]
            d = slp.reshape_wavelet(gs[0] - water["FS"][0], F_SRC, fw)
            c, live = slp.nmo(d, gs[1], rays, mode, mdeg, smax)
            ST["FSsym_trend_%s" % vname] = slp.stacks(c, live, gs[2], VRANGE.get(vname))
    test_mode = "G_trend_main" not in ST
    if test_mode:                                         # FD gathers not there yet: let the Zoeppritz synthetic stand in (pipeline test only)
        for case in CASES:
            ST["G_%s_main" % case] = ST["C_%s_main" % case]
            ST["FS_%s_main" % case] = ST["CM_%s_main" % case]
        print("TEST MODE: no FD gathers found, the primaries-only synthetic stands in for G")
    np.savez_compressed(os.path.join(ROOT, "stacks.npz"), t=slp.T_OUT,
                        **{"%s__%s" % (k, kk): v for k, S in ST.items() for kk, v in S.items()})

    R = dict(test_mode=test_mode, model=m0["info"], runs=INFO, variants={k: list(v) for k, v in VARIANTS.items()}, cases={})
    # ---- checks of the kinematics
    chk = {}
    i_soft = [int(np.searchsorted(m0["ze"], z)) - 1 for z in (1358.0, 1390.0)]
    chk["t0_soft_layer_ms"] = [float(1000 * rays.t0[i]) for i in i_soft]
    for x in (150.0, 750.0, 1775.0, 2500.0, 3150.0):
        ix = int(np.argmin(np.abs(off - x)))
        chk["angle_at_soft_top_deg_x%d" % x] = float(rays.ang[i_soft[0], ix])
    # hyperbolic NMO residual of primaries (exact time against Vrms hyperbola) at the soft layer
    k = i_soft[0]
    th = np.sqrt(rays.t0[k] ** 2 + (off / rays.vrms[k]) ** 2)
    res = 1000 * (rays.t[k] - th)
    chk["hyperbolic_residual_ms_at_soft_top"] = {str(int(x)): float(res[int(np.argmin(np.abs(off - x)))]) for x in (150, 760, 1760, 2500, 3150)}
    # flatness of the strongest primary after NMO in the FD data (lag of each offset trace against the near stack)
    for nm, (c, live) in NMOG.items():
        cfg, vname = nm.split("_")
        S0 = ST["%s_trend_%s" % (cfg, vname)]
        w = (slp.T_OUT > 1.2) & (slp.T_OUT < 1.5)
        ic = int(np.where(w)[0][np.argmax(np.abs(S0["nea"][w]))])
        fl = {}
        for x in (300, 750, 1000, 1500, 2000, 2500, 3000):
            ix = int(np.argmin(np.abs(off - x)))
            if live[ix, ic - 30:ic + 31].all():
                fl[str(x)] = round(1000 * slp.DT * slp.lag_to(c[ix].astype(float), S0["nea"], ic, 16, 12), 2)
        chk["nmo_lag_ms_vs_near_stack_%s_event_%dms" % (nm, int(1000 * ic * slp.DT))] = fl
    # time axis: sea-bed reflection on the 150 m trace (ray time from 6 m to the sea bed at 81 m and up to 8 m; the
    # ghosts centre the recorded pulse 9.5 ms later), and its first water-layer multiple
    for cfg, (d, t) in RAW.items():
        tr = d[0]
        t_sb = np.sqrt(150.0 ** 2 + 162.0 ** 2) / slp.VW          # surface-to-surface ray: where the ghosted pulse is centred
        w = (t > t_sb - 0.03) & (t < t_sb + 0.03)
        chk["seabed_150m_%s" % cfg] = dict(expected_centre_ms=round(1000 * t_sb, 1), zero_offset_twt_ms=109.5,
                                           picked_ms=round(float(1000 * t[w][np.argmax(np.abs(tr[w]))]), 1))
    for k in ("nea", "mid", "far"):
        w = (slp.T_OUT >= 1.1) & (slp.T_OUT <= 1.9)
        a, b = ST.get("C_trend_main"), ST.get("G_trend_main")
        if a is not None and b is not None and not test_mode:
            chk["corr_zoeppritz_synthetic_vs_fd_G_%s" % k] = round(float(np.corrcoef(a[k][w], b[k][w])[0, 1]), 3)
        a, b = ST.get("CM_trend_main"), ST.get("FS_trend_main")
        if a is not None and b is not None:
            chk["corr_ray_multiples_synthetic_vs_fd_FS_%s" % k] = round(float(np.corrcoef(a[k][w], b[k][w])[0, 1]), 3)
    # tie to the real full-stack trace at the well (stored in the hackathon app's wells.json; optional)
    wj = os.environ.get("SLEIPNER_WELLS_JSON", "/mnt/ssd/Data/force-hackathon-2026/sleipner/cubes/sleipner94/wells.json")
    if os.path.exists(wj):
        seis = np.array(json.load(open(wj))["wells"][0]["seis"]["trace"], float)
        tr_t = np.arange(len(seis)) * 0.002
        real = np.interp(slp.T_OUT, tr_t, seis, right=0.0)
        tie = {}
        for nm in ("G_trend_main", "FS_trend_main", "C_trend_main", "CM_trend_main"):
            if nm not in ST:
                continue
            syn = -ST[nm]["ful"]                          # the ghosts turn the pulse upside down; the real data are zero-phased
            row = {}
            for lo, hi in ((1.10, 1.40), (1.38, 1.56), (1.10, 1.90)):
                w = (slp.T_OUT >= lo) & (slp.T_OUT <= hi)
                best = max(((float(np.corrcoef(real[w], np.roll(syn, L)[w])[0, 1]), L) for L in range(-30, 31)), key=lambda v: v[0])
                row["%d-%d ms" % (1000 * lo, 1000 * hi)] = dict(corr_zero_lag=round(float(np.corrcoef(real[w], syn[w])[0, 1]), 2), best_corr=round(best[0], 2), best_lag_ms=best[1])
            tie[nm] = row
        chk["tie_to_real_full_stack_at_well"] = tie
    # one interface alone (base of the single soft layer 1358-1390 m against the 30 m below, blocked media):
    # stack-averaged Zoeppritz reflection coefficient, far over near, for Vs 0 / 7 / 15 / 20 % lower in the layer
    def blk(a, b):
        cc = (m0["z"] >= a) & (m0["z"] < b)
        return float(1 / np.mean(1 / m0["vp"][cc])), float(m0["vs"][cc].mean()), float(m0["rho"][cc].mean())
    up, lo = blk(1358.0, 1390.0), blk(1390.0, 1420.0)
    kb = i_soft[1]
    si = {}
    for X in (0.0, 0.07, 0.15, 0.20):
        r = {}
        for k, (a, b) in slp.RANGES_MAIN.items():
            mm = (off > a) & (off <= b) & (rays.ang[kb] <= slp.MUTE_DEG)
            pp = rays.p[kb, mm]
            A = lambda v: np.full(len(pp), float(v))
            r[k] = float(np.mean(np.real(slp.zoep_pp(A(up[0]), A(up[1] * (1 - X)), A(up[2]), A(lo[0]), A(lo[1]), A(lo[2]), pp))))
        si["vs_minus_%02d_pct" % round(100 * X)] = dict(R_near=r["nea"], R_mid=r["mid"], R_far=r["far"], far_over_near=r["far"] / r["nea"])
    chk["single_interface_base_of_soft_layer"] = dict(upper_vp_vs_rho=up, lower_vp_vs_rho=lo, by_case=si)
    # effect of 15 % lower Vs in the single layer (difference to the trend case), FD without sea surface against the
    # primaries-only Zoeppritz synthetic, 1380-1500 ms. Each difference is given relative to the RMS of its own stack.
    if all(kx in ST for kx in ("G_m15_main", "G_trend_main", "C_m15_main", "C_trend_main")) and not test_mode:
        w = (slp.T_OUT >= 1.38) & (slp.T_OUT <= 1.50)
        rw = lambda x: float(np.sqrt(np.mean(x[w] ** 2)))
        lv = {}
        for k in ("nea", "mid", "far"):
            dF, dC = ST["G_m15_main"][k] - ST["G_trend_main"][k], ST["C_m15_main"][k] - ST["C_trend_main"][k]
            cF, cC = 100 * rw(dF) / rw(ST["G_trend_main"][k]), 100 * rw(dC) / rw(ST["C_trend_main"][k])
            lv[k] = dict(corr=round(float(np.corrcoef(dF[w], dC[w])[0, 1]), 3), change_pct_of_stack_rms_fd=round(cF, 1),
                         change_pct_of_stack_rms_zoeppritz=round(cC, 1), strength_fd_over_zoeppritz=round(cF / cC, 3))
        chk["lowvs15_single_layer_difference_fd_vs_zoeppritz_1380_1500ms"] = lv
    R["checks"] = chk

    for vname in VARIANTS:
        key = lambda cfg, case: "%s_%s_%s" % (cfg, case, vname)
        if key("G", "trend") not in ST:
            continue
        out = {}
        G = ST[key("G", "trend")]
        ev = {"G": slp.pick_events(G["ful"])}
        FS = ST.get(key("FS", "trend"))
        if FS is not None:
            ev["FS"] = slp.pick_events(FS["ful"])
        out["events_ms"] = {k: [int(round(1000 * i * slp.DT)) for i in v] for k, v in ev.items()}
        data = {}
        for cfg in ("G", "FS", "C", "CM", "FSdx1", "FSsym"):
            for case in CASES_C + ("trend_sbl", "trend_sbh"):   # CASES_C holds the soft-bed cases too
                kk = key(cfg, case) if cfg in ("G", "FS") else "%s_%s_%s" % (cfg, case, vname)
                if kk in ST:
                    data["%s_%s" % (cfg, case)] = ST[kk]
        # stack-domain multiple filter (as the parent session's fm_sl_demult.py), coefficients measured on the data itself
        filt = {}
        if FS is not None:
            # parents: the two strongest reflectors of the full stack 30-160 ms above the zone (the analogue of the
            # events at 1276 and 1316 ms that the parent session used)
            evs = [i for i in ev["FS"] if 1.24 <= i * slp.DT <= 1.38]
            par = sorted(evs, key=lambda i: -abs(FS["ful"][i]))[:2]
            win = ((0.080, 0.150), (0.190, 0.300))
            for k in ("nea", "mid", "far"):
                c1 = [slp.copy_coef(FS[k], i, *win[0], sign=-1) for i in par]
                c2 = [slp.copy_coef(FS[k], i, *win[1], sign=+1) for i in par]
                filt[k] = dict(a1=float(np.mean([c[0] for c in c1])), lag1_ms=float(1000 * np.mean([c[1] for c in c1])),
                               a2=float(np.mean([c[0] for c in c2])), lag2_ms=float(1000 * np.mean([c[1] for c in c2])),
                               per_parent=[dict(t_ms=1000 * i * slp.DT, c1=round(a[0], 3), lag1_ms=round(1000 * a[1], 1), c2=round(b[0], 3), lag2_ms=round(1000 * b[1], 1)) for i, a, b in zip(par, c1, c2)])
            out["filter_from_data"] = filt
            for case in CASES_FD + ("dx1",):
                kk = "FS_%s" % case if case != "dx1" else "FSdx1_trend"
                if kk in data:
                    S = {k: slp.three_point_filter(data[kk][k], filt[k]["a1"], filt[k]["lag1_ms"] / 1000, filt[k]["a2"], filt[k]["lag2_ms"] / 1000) for k in ("nea", "mid", "far")}
                    S["ful"] = sum(slp.W_FULL[k] * S[k] for k in slp.W_FULL)
                    data[kk + "_filt"] = S
            # the multiples alone = FS - G (same ghosts in both): strength per stack and window, and the true copy
            Mt = {k: FS[k] - G[k] for k in ("nea", "mid", "far")}
            tt = slp.T_OUT
            zone = (tt >= slp.ZONE[0]) & (tt <= slp.ZONE[1])
            ref = ((tt >= slp.T_EV[0]) & (tt < slp.REF_EXCL[0])) | ((tt > slp.REF_EXCL[1]) & (tt <= slp.T_EV[1]))
            rms = lambda x, w: float(np.sqrt(np.mean(x[w] ** 2)))
            out["multiples_over_primaries_rms"] = {k: dict(zone=rms(Mt[k], zone) / rms(G[k], zone), reference=rms(Mt[k], ref) / rms(G[k], ref)) for k in Mt}
            out["rms_far_over_near"] = dict(G=dict(zone=rms(G["far"], zone) / rms(G["nea"], zone), reference=rms(G["far"], ref) / rms(G["nea"], ref)),
                                            FS=dict(zone=rms(FS["far"], zone) / rms(FS["nea"], zone), reference=rms(FS["far"], ref) / rms(FS["nea"], ref)))
            out["rms_mid_over_near"] = dict(G=dict(zone=rms(G["mid"], zone) / rms(G["nea"], zone), reference=rms(G["mid"], ref) / rms(G["nea"], ref)),
                                            FS=dict(zone=rms(FS["mid"], zone) / rms(FS["nea"], zone), reference=rms(FS["mid"], ref) / rms(FS["nea"], ref)))
            # true copies: the multiples-only trace (FS - G) against the primary wavelet of the strongest reflectors
            gp = sorted([i for i in ev["G"] if 1.15 <= i * slp.DT <= 1.60], key=lambda i: -abs(G["ful"][i]))[:4]
            tc = []
            for i in sorted(gp):
                row = dict(parent_ms=1000 * i * slp.DT)
                for k in ("nea", "mid", "far"):
                    f1 = slp.copy_coef(G[k], i, *win[0], sign=-1, y=Mt[k])
                    f2 = slp.copy_coef(G[k], i, *win[1], sign=+1, y=Mt[k])
                    row[k] = dict(first=[round(f1[0], 3), round(1000 * f1[1], 1)], second=[round(f2[0], 3), round(1000 * f2[1], 1)])
                tc.append(row)
            out["true_copies_of_strongest_primaries"] = tc
            mo = {}
            for case in CASES_FD:
                a, b = data.get("FS_%s" % case), data.get("G_%s" % case)
                if a is not None and b is not None:
                    mo[case] = {k: dict(zone=rms(a[k] - b[k], zone) / rms(b[k], zone), reference=rms(a[k] - b[k], ref) / rms(b[k], ref)) for k in ("nea", "mid", "far")}
            out["multiples_over_primaries_rms_by_case"] = mo
            if "FSdx1_trend" in data:
                w19 = (tt >= 1.1) & (tt <= 1.9)
                A, B = data["FS_trend"], data["FSdx1_trend"]
                out["grid_convergence_2m_vs_1m"] = {k: dict(corr=float(np.corrcoef(A[k][w19], B[k][w19])[0, 1]),
                                                            nrms_pct=float(200 * rms(A[k] - B[k], w19) / (rms(A[k], w19) + rms(B[k], w19))),
                                                            rms_ratio_2m_over_1m=float(rms(A[k], w19) / rms(B[k], w19))) for k in ("nea", "mid", "far")}
        if "FSsym_trend" in data and "FS_trend" in data:   # exact mirror against the standard set-up (same grid, same source)
            w19 = (slp.T_OUT >= 1.1) & (slp.T_OUT <= 1.9)
            rr = lambda x: float(np.sqrt(np.mean(x[w19] ** 2)))
            A, B = data["FS_trend"], data["FSsym_trend"]
            out["exact_mirror_vs_standard"] = {k: dict(corr=float(np.corrcoef(A[k][w19], B[k][w19])[0, 1]), nrms_pct=float(200 * rr(A[k] - B[k]) / (rr(A[k]) + rr(B[k]))),
                                                       rms_ratio_standard_over_exact=float(rr(A[k]) / rr(B[k]))) for k in ("nea", "mid", "far")}
        # multiples scaled per stack to the copy strength of the real data: primaries + SCALE * (with surface - without)
        for case in CASES_FD:
            a, b = data.get("FS_%s" % case), data.get("G_%s" % case)
            if a is not None and b is not None:
                Sx = {k: b[k] + SCALE[k] * (a[k] - b[k]) for k in ("nea", "mid", "far")}
                Sx["ful"] = sum(slp.W_FULL[k] * Sx[k] for k in slp.W_FULL)
                data["FSs_%s" % case] = Sx
        # ---- per-event ratios and zone / reference
        tt = slp.T_OUT
        zone_w = (tt >= slp.ZONE[0]) & (tt <= slp.ZONE[1])
        ref_w = ((tt >= slp.T_EV[0]) & (tt < slp.REF_EXCL[0])) | ((tt > slp.REF_EXCL[1]) & (tt <= slp.T_EV[1]))
        rms_ = lambda x, w: float(np.sqrt(np.mean(np.asarray(x)[w] ** 2)))
        for name, S in data.items():
            res = {}
            fz, fr = rms_(S["far"], zone_w) / rms_(S["nea"], zone_w), rms_(S["far"], ref_w) / rms_(S["nea"], ref_w)
            mz, mr = rms_(S["mid"], zone_w) / rms_(S["nea"], zone_w), rms_(S["mid"], ref_w) / rms_(S["nea"], ref_w)
            res["rms"] = dict(far_near_zone=fz, far_near_ref=fr, far_near_zone_over_ref=fz / fr, mid_near_zone=mz, mid_near_ref=mr, mid_near_zone_over_ref=mz / mr,
                              near_zone_over_ref=rms_(S["nea"], zone_w) / rms_(S["nea"], ref_w), mid_zone_over_ref=rms_(S["mid"], zone_w) / rms_(S["mid"], ref_w), far_zone_over_ref=rms_(S["far"], zone_w) / rms_(S["far"], ref_w))
            for evname, idx in ev.items():
                E = [slp.measure_event(S, i) for i in idx]
                res["events_" + evname] = dict(summary=slp.summarise(E),
                                               table=[{k: (round(v, 3) if isinstance(v, float) else v) for k, v in e.items()} for e in E])
            out[name] = res
        R["cases"][vname] = out
    json.dump(R, open(os.path.join(ROOT, "results.json"), "w"), indent=1, default=float)

    # ---- printed summary
    print(json.dumps(R["checks"], indent=1))
    for vname, out in R["cases"].items():
        print("\n================ variant %s: %s" % (vname, VARIANTS[vname]))
        print("events picked (ms):", out["events_ms"])
        for kx in ("filter_from_data", "multiples_over_primaries_rms", "rms_far_over_near", "rms_mid_over_near", "true_copies_of_strongest_primaries"):
            if kx in out:
                print(kx, json.dumps(out[kx], default=float))
        print("%-16s %-3s | far/near: zone  ref  zone/ref (n zone, n ref, bad) | mid/near: zone  ref  zone/ref | peak far/near zone/ref | peak mid/near zone/ref" % ("data", "ev"))
        for name, res in out.items():
            if not isinstance(res, dict) or "events_G" not in res:
                continue
            for evname in ("G", "FS"):
                if "events_" + evname not in res:
                    continue
                s = res["events_" + evname]["summary"]
                f = lambda d, k: ("%5.2f" % d[k]) if d[k] is not None else "  nan"
                print("%-16s %-3s |          %s %s  %s   (%d, %d, %d) |           %s %s  %s |          %s          |        %s" % (
                    name, evname, f(s["fn"], "zone_median"), f(s["fn"], "ref_median"), f(s["fn"], "zone_over_ref"), s["n_zone"], s["n_ref"], s["fn"]["n_bad"],
                    f(s["mn"], "zone_median"), f(s["mn"], "ref_median"), f(s["mn"], "zone_over_ref"), f(s["pk_fn"], "zone_over_ref"), f(s["pk_mn"], "zone_over_ref")))
        print("RMS in windows (zone 1408-1492 ms against reference 1100-1380 and 1560-1900 ms):")
        for name, res in out.items():
            if isinstance(res, dict) and "rms" in res:
                r = res["rms"]
                print("  %-16s far/near zone %.2f ref %.2f zone/ref %.2f | mid/near zone %.2f ref %.2f zone/ref %.2f | near zone/ref %.2f far zone/ref %.2f" % (
                    name, r["far_near_zone"], r["far_near_ref"], r["far_near_zone_over_ref"], r["mid_near_zone"], r["mid_near_ref"], r["mid_near_zone_over_ref"], r["near_zone_over_ref"], r["far_zone_over_ref"]))


if __name__ == "__main__":
    main()
