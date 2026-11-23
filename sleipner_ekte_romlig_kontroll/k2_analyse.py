#!/usr/bin/env python3
"""Independent check of 'sleipner_ekte_romlig', step 2. Own code on own cell stacks (cells.npz from k1_cells.py).
Blocks of 2 x 3 cells (200 m x 225 m) and, as a grid check, 3 x 4 cells (300 m x 300 m).
Time of sample k in the cell arrays: 1000 + 2 k ms (cube sample index 500 + k, cube starts at 0 ms)."""
import json
import sys
import numpy as np

OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig_kontroll"
DT, T0 = 2.0, 1000.0
CI, CX = 8, 6
WELL_CI, WELL_CX = 57 // CI, 143 // CX
rng = np.random.default_rng(20261119)
Cc = np.load(OUT + "/cells.npz")


def si(t):
    return int(round((t - T0) / DT))


def sl(w, sh=0):
    return slice(si(w[0]) + sh, si(w[1]) + 1 + sh)


def r4(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 4)


def pear(x, y):
    x = x - x.mean()
    y = y - y.mean()
    return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum()))


def rank(x):
    r = np.empty(len(x))
    r[np.argsort(x)] = np.arange(len(x))
    return r


def spear(x, y):
    return pear(rank(x), rank(y))


def build(proc, bi, bx):
    """Block mean traces and two checkerboard halves, flattened on the package 1236-1346 ms with the well block as datum."""
    n = Cc[f"{proc}_nea_n"].astype(float)
    NBI, NBX = n.shape[0] // bi, n.shape[1] // bx
    chk = (np.add.outer(np.arange(bi), np.arange(bx)) % 2).astype(bool)
    G = {}
    for k in ("nea", "mid", "far", "ful"):
        c = Cc[f"{proc}_{k}"].astype(np.float64)[:NBI * bi, :NBX * bx]
        w = n[:NBI * bi, :NBX * bx]
        cw = (c * w[..., None]).reshape(NBI, bi, NBX, bx, -1)
        ww = w.reshape(NBI, bi, NBX, bx)
        M = cw.sum((1, 3)) / np.maximum(ww.sum((1, 3)), 1)[..., None]
        A = cw[:, chk[:, None, :].repeat(1, 1)[:, 0, :]] if False else None
        a_s = np.einsum("ibjxt,bx->ijt", cw, (~chk).astype(float)) / np.maximum(np.einsum("ibjx,bx->ij", ww, (~chk).astype(float)), 1)[..., None]
        b_s = np.einsum("ibjxt,bx->ijt", cw, chk.astype(float)) / np.maximum(np.einsum("ibjx,bx->ij", ww, chk.astype(float)), 1)[..., None]
        G[k] = (M.reshape(NBI * NBX, -1), a_s.reshape(NBI * NBX, -1), b_s.reshape(NBI * NBX, -1))
    cnt = n[:NBI * bi, :NBX * bx].reshape(NBI, bi, NBX, bx)
    full = bi * bx * CI * CX
    valid = ((cnt.sum((1, 3)) >= 0.9 * full) & (cnt.min((1, 3)) > 0)).ravel()
    # flatten: integer shift from cross-correlation of the full-stack block trace with the area mean, 1236-1346 ms
    F = G["ful"][0]
    w0, w1 = si(1236), si(1346) + 1
    lags = np.arange(-25, 26)
    sh = np.zeros(len(F))
    ref = F[valid].mean(0)
    for it in range(3):
        cc = np.array([[np.dot(F[b, w0 + L:w1 + L], ref[w0:w1]) for L in lags] for b in range(len(F))])
        j = cc.argmax(1)
        jj = np.clip(j, 1, len(lags) - 2)
        y0, y1, y2 = cc[np.arange(len(F)), jj - 1], cc[np.arange(len(F)), jj], cc[np.arange(len(F)), jj + 1]
        d = np.where((j == jj) & (y0 - 2 * y1 + y2 != 0), 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 + 1e-30), 0.0)
        sh = lags[j] + d
        ish = np.rint(sh).astype(int)
        ref = np.mean([np.roll(F[b], -ish[b]) for b in np.where(valid)[0]], 0)
    wb = (WELL_CI // bi) * NBX + (WELL_CX // bx)
    sh = sh - sh[wb]
    ish = np.rint(sh).astype(int)
    valid &= np.abs(ish) <= 24
    Fl = {k: tuple(np.array([np.roll(x[b], -ish[b]) for b in range(len(x))]) for x in G[k]) for k in G}
    return Fl, valid, sh, NBI, NBX, wb


def residual_shift(F, valid, w, maxlag=10):
    """Sub-sample shift (ms) of each flattened block trace against the area mean in window w."""
    ref = F[valid].mean(0)
    w0, w1 = si(w[0]), si(w[1]) + 1
    lags = np.arange(-maxlag, maxlag + 1)
    cc = np.array([[np.dot(F[b, w0 + L:w1 + L], ref[w0:w1]) for L in lags] for b in range(len(F))])
    j = np.clip(cc.argmax(1), 1, len(lags) - 2)
    i = np.arange(len(F))
    y0, y1, y2 = cc[i, j - 1], cc[i, j], cc[i, j + 1]
    d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2 - 1e-30)
    edge = (cc.argmax(1) == 0) | (cc.argmax(1) == len(lags) - 1)
    return (lags[j] + np.clip(d, -1, 1)) * DT, edge


def energy(T, wins, mode):
    """Energy per block summed over a list of windows. mode 'rms': block mean squared; 'xp': product of the halves."""
    M, A, B = T
    e = 0.0
    for w in wins:
        s = sl(w)
        e = e + ((M[:, s] ** 2).sum(1) if mode == "rms" else (A[:, s] * B[:, s]).sum(1) if mode == "xp"
                 else (A[:, s] ** 2).sum(1) if mode == "A" else (B[:, s] ** 2).sum(1))
    return e


def copy_fit(X, src, lag_ms=(96, 128)):
    """Least-squares coefficient of X(t+L) on X(t), t in src. Returns per-block num, den, tt for each lag."""
    s = sl(src)
    S = X[:, s]
    Ls = np.arange(int(lag_ms[0] / DT), int(lag_ms[1] / DT) + 1)
    num = np.array([(X[:, s.start + L:s.stop + L] * S).sum(1) for L in Ls]).T
    tt = np.array([(X[:, s.start + L:s.stop + L] ** 2).sum(1) for L in Ls]).T
    den = (S * S).sum(1)
    return Ls * DT, num, den, tt


def superblocks(NBI, NBX, ok, si_=3, sx_=4):
    g = (np.arange(NBI)[:, None] // si_) * 100 + (np.arange(NBX)[None, :] // sx_)
    g = g.ravel()
    return [np.where((g == u) & ok)[0] for u in np.unique(g) if ((g == u) & ok).any()]


def boot(fn, groups, nrep=1000):
    v = []
    for _ in range(nrep):
        pick = rng.integers(0, len(groups), len(groups))
        idx = np.concatenate([groups[p] for p in pick])
        v.append(fn(idx))
    v = np.array(v, float)
    return [r4(np.nanpercentile(v, 2.5)), r4(np.nanpercentile(v, 97.5))]


def lattice(x, y, ok, NBI, NBX, step=3):
    ok2 = ok.reshape(NBI, NBX)
    x2, y2 = x.reshape(NBI, NBX), y.reshape(NBI, NBX)
    rs, ns = [], []
    for a in range(step):
        for b in range(step):
            m = ok2[a::step, b::step]
            if m.sum() >= 10:
                rs.append(pear(x2[a::step, b::step][m], y2[a::step, b::step][m]))
                ns.append(int(m.sum()))
    return dict(mean_r=r4(np.mean(rs)), min_r=r4(np.min(rs)), max_r=r4(np.max(rs)), n_mean=r4(np.mean(ns)))


def corr_block(x, y, ok, geo, ci=True):
    NBI, NBX, groups = geo
    m = ok & np.isfinite(x) & np.isfinite(y)
    o = dict(n=int(m.sum()), pearson=r4(pear(x[m], y[m])), spearman=r4(spear(x[m], y[m])))
    if ci:
        gr = [g[m[g]] for g in groups]
        gr = [g for g in gr if len(g)]
        o["pearson_ci95"] = boot(lambda i: pear(x[i], y[i]), gr)
        o["lattice3"] = lattice(x, y, m, NBI, NBX, 3)
    return o


def analyse(proc, bi, bx, full=True):
    Fl, valid, sh, NBI, NBX, wb = build(proc, bi, bx)
    nb = len(valid)
    res = dict(block_cells=[bi, bx], block_m=[bi * CI * 12.5, bx * CX * 12.5], n_blocks=int(nb), n_valid=int(valid.sum()),
               area_km2=r4(valid.sum() * bi * CI * 12.5 * bx * CX * 12.5 / 1e6))
    res["flatten_shift_ms"] = dict(sd=r4(np.std(sh[valid]) * DT), p5=r4(np.percentile(sh[valid], 5) * DT), p95=r4(np.percentile(sh[valid], 95) * DT))
    Z, RA, RB = (1400, 1500), (1150, 1360), (1580, 1880)
    ok = valid.copy()
    groups = superblocks(NBI, NBX, ok)
    geo = (NBI, NBX, groups)
    E = {m: {k: {nm: energy(Fl[k], w, m) for nm, w in (("Z", [Z]), ("RA", [RA]), ("RB", [RB]), ("R", [RA, RB]))} for k in ("nea", "mid", "far")}
         for m in ("rms", "xp", "A", "B")}

    def gD(idx, a="far", b="nea", ref="R", mode="rms", zone="Z", EE=None):
        e = E[mode] if EE is None else EE
        s = lambda k, nm: e[k][nm][idx].sum()
        v = (s(a, zone) / s(b, zone)) / (s(a, ref) / s(b, ref))
        return float(np.sqrt(v)) if v > 0 else np.nan

    def gL(idx, k, ref="RB", mode="rms", zone="Z", EE=None):
        e = E[mode] if EE is None else EE
        v = e[k][zone][idx].sum() / e[k][ref][idx].sum()
        return float(np.sqrt(v)) if v > 0 else np.nan

    iv = np.where(ok)[0]
    gl = {}
    for mode in ("rms", "xp"):
        for a, b in (("far", "nea"), ("mid", "nea"), ("far", "mid")):
            for ref in ("R", "RA", "RB"):
                key = f"{a}_{b}_ref{ref}_{mode}"
                gl[key] = dict(value=r4(gD(iv, a, b, ref, mode)))
                if full and ref != "RA":
                    gl[key]["ci95"] = boot(lambda i: gD(i, a, b, ref, mode), groups)
    lnDb = 0.5 * (np.log(E["rms"]["far"]["Z"] / E["rms"]["nea"]["Z"]) - np.log(E["rms"]["far"]["R"] / E["rms"]["nea"]["R"]))
    lnDbB = 0.5 * (np.log(E["rms"]["far"]["Z"] / E["rms"]["nea"]["Z"]) - np.log(E["rms"]["far"]["RB"] / E["rms"]["nea"]["RB"]))
    lnMb = 0.5 * (np.log(E["rms"]["mid"]["Z"] / E["rms"]["nea"]["Z"]) - np.log(E["rms"]["mid"]["R"] / E["rms"]["nea"]["R"]))
    gl["far_near_block_median_rms"] = r4(np.exp(np.median(lnDb[ok])))
    gl["mid_near_block_median_rms"] = r4(np.exp(np.median(lnMb[ok])))
    gl["far_near_well_block_rms"] = r4(np.exp(lnDb[wb]))
    gl["frac_blocks_far_near_below_0p8"] = r4(np.mean(np.exp(lnDb[ok]) < 0.8))
    res["global"] = gl

    # ---- copy of the package 1256-1356 ms
    SRC = (1256, 1356)
    cp = {}
    fits = {}
    for k in ("nea", "mid", "far"):
        Ls, num, den, tt = copy_fit(Fl[k][0], SRC)
        fits[k] = (Ls, num, den, tt)
        bg = num[ok].sum(0) / den[ok].sum()
        j = int(np.argmin(bg))
        jj = min(max(j, 1), len(Ls) - 2)
        d = 0.5 * (bg[jj - 1] - bg[jj + 1]) / (bg[jj - 1] - 2 * bg[jj] + bg[jj + 1])
        rg = num[ok].sum(0) / np.sqrt(den[ok].sum() * tt[ok].sum(0))
        cp[k] = dict(coef=r4(bg[j]), lag_ms=r4(Ls[jj] + d * DT), corr_at_lag=r4(rg[j]))
    res["copy_global"] = cp
    Ls, num, den, tt = fits["nea"]
    bb = num / den[:, None]
    jb = bb.argmin(1)
    i_ = np.arange(nb)
    b_blk, r_blk = bb[i_, jb], num[i_, jb] / np.sqrt(den * tt[i_, jb])
    amp_blk = np.sqrt(tt[i_, jb] / den)
    okc = ok & (b_blk < 0)
    res["copy_blocks"] = dict(n_negative=int(okc.sum()), coef_median=r4(np.median(b_blk[okc])), coef_p5=r4(np.percentile(-b_blk[okc], 5)),
                              coef_p95=r4(np.percentile(-b_blk[okc], 95)), corr_median=r4(np.median(r_blk[okc])),
                              corr_p5=r4(np.percentile(-r_blk[okc], 5)), corr_p95=r4(np.percentile(-r_blk[okc], 95)),
                              lag_median_ms=r4(np.median(Ls[jb][okc])), sd_ln_coef=r4(np.std(np.log(-b_blk[okc]))),
                              sd_ln_corr=r4(np.std(np.log(-r_blk[okc]))), sd_ln_amp_ratio=r4(np.std(np.log(amp_blk[okc]))))
    lnb = np.where(okc, np.log(np.where(okc, -b_blk, 1.0)), np.nan)
    lnr = np.where(okc, np.log(np.where(okc, -r_blk, 1.0)), np.nan)
    lna = np.log(amp_blk)

    # ---- Test 2 reproduced and taken apart
    lvl = {k: 0.5 * np.log(E["rms"][k]["Z"] / E["rms"][k]["RB"]) for k in ("nea", "mid", "far")}
    T2 = {}
    T2["lnD_vs_ln_coef"] = corr_block(lnDb, lnb, okc, geo, full)
    T2["lnD_refbelow_vs_ln_coef"] = corr_block(lnDbB, lnb, okc, geo, full)
    T2["lnD_vs_ln_corr_shape_only"] = corr_block(lnDb, lnr, okc, geo, full)
    T2["lnD_vs_ln_amp_ratio_no_shape"] = corr_block(lnDb, lna, okc, geo, full)
    T2["lnD_vs_minus_corr_linear"] = corr_block(lnDb, -r_blk, ok, geo, full)
    T2["lnD_refbelow_vs_near_level"] = corr_block(lnDbB, lvl["nea"], ok, geo, full)
    T2["lnD_refbelow_vs_far_level"] = corr_block(lnDbB, lvl["far"], ok, geo, full)
    T2["near_level_vs_ln_coef"] = corr_block(lvl["nea"], lnb, okc, geo, full)
    T2["mid_level_vs_ln_coef"] = corr_block(lvl["mid"], lnb, okc, geo, full)
    T2["far_level_vs_ln_coef"] = corr_block(lvl["far"], lnb, okc, geo, full)
    T2["near_level_vs_ln_corr_shape_only"] = corr_block(lvl["nea"], lnr, okc, geo, full)
    T2["far_level_vs_ln_corr_shape_only"] = corr_block(lvl["far"], lnr, okc, geo, full)
    T2["mid_level_vs_ln_corr_shape_only"] = corr_block(lvl["mid"], lnr, okc, geo, full)
    T2["mid_near_vs_ln_coef"] = corr_block(lnMb, lnb, okc, geo, full)
    # copy measured on MID (other data than near) as predictor
    Lm, numm, denm, ttm = fits["mid"]
    bm = (numm / denm[:, None]).min(1)
    okm = ok & (bm < 0)
    T2["lnD_vs_ln_coef_measured_on_mid"] = corr_block(lnDb, np.where(okm, np.log(np.where(okm, -bm, 1.0)), np.nan), okm, geo, full)
    # wrong lag (150 +- 16 ms): largest |coef|
    Lw, numw, denw, ttw = copy_fit(Fl["nea"][0], SRC, (134, 166))
    bw = np.abs(numw / denw[:, None]).max(1)
    T2["lnD_vs_ln_abs_coef_wrong_lag_150"] = corr_block(lnDb, np.log(bw), ok, geo, full)

    # groups by copy coefficient: terciles and weakest decile, same data and split halves (select on one half, measure on the other)
    def groups_by(pred, m, EE=None, mode="rms", zone="Z", ref="R", refl="RB"):
        idx = np.where(m)[0]
        o = np.argsort(pred[idx])
        n = len(idx)
        sel = dict(weak_tercile=idx[o[:n // 3]], strong_tercile=idx[o[-(n // 3):]], weak_decile=idx[o[:n // 10]], strong_decile=idx[o[-(n // 10):]])
        out = {}
        for nm, ii in sel.items():
            out[nm] = dict(far_near=r4(gD(ii, "far", "nea", ref, mode, zone, EE)), mid_near=r4(gD(ii, "mid", "nea", ref, mode, zone, EE)),
                           near=r4(gL(ii, "nea", refl, mode, zone, EE)), mid=r4(gL(ii, "mid", refl, mode, zone, EE)), far=r4(gL(ii, "far", refl, mode, zone, EE)))
        out["all"] = dict(far_near=r4(gD(idx, "far", "nea", ref, mode, zone, EE)), mid_near=r4(gD(idx, "mid", "nea", ref, mode, zone, EE)))
        return out
    T2["groups_same_data"] = groups_by(lnb, okc)
    # split halves: copy coefficient from half A (at the block's own best lag), far/near from half B, and the other way round
    half = {}
    for hsel, hmeas, tag in ((1, "B", "select_A_measure_B"), (2, "A", "select_B_measure_A")):
        _, nu, de, _t = copy_fit(Fl["nea"][hsel], SRC)
        bh = (nu / de[:, None]).min(1)
        mh = ok & (bh < 0)
        half[tag] = groups_by(np.where(mh, np.log(np.where(mh, -bh, 1.0)), np.nan), mh, mode=hmeas)
        ampA = np.sqrt(energy(Fl["nea"], [(1366, 1466)], "A" if hsel == 1 else "B") / energy(Fl["nea"], [SRC], "A" if hsel == 1 else "B"))
        half[tag + "_by_amp_ratio_no_shape"] = groups_by(np.log(ampA), ok, mode=hmeas)
    T2["groups_split_halves"] = half
    T2["groups_same_data_by_amp_ratio_no_shape"] = groups_by(lna, ok)
    res["test2"] = T2

    # ---- placebo: the same construction on window pairs away from the zone
    PL = {}
    for T in ((1580, 1680), (1680, 1780), (1780, 1880)):
        src = (T[0] - 144, T[1] - 144)
        refw = [RA] + [w for w in ((1580, T[0]), (T[1], 1880)) if w[1] - w[0] >= 40]
        refb = [w for w in ((1580, T[0]), (T[1], 1880)) if w[1] - w[0] >= 40]
        EE = {k: {"Z": energy(Fl[k], [T], "rms"), "R": energy(Fl[k], refw, "rms"), "RB": energy(Fl[k], refb, "rms")} for k in ("nea", "mid", "far")}
        lD = 0.5 * (np.log(EE["far"]["Z"] / EE["nea"]["Z"]) - np.log(EE["far"]["R"] / EE["nea"]["R"]))
        lDB = 0.5 * (np.log(EE["far"]["Z"] / EE["nea"]["Z"]) - np.log(EE["far"]["RB"] / EE["nea"]["RB"]))
        lv = {k: 0.5 * np.log(EE[k]["Z"] / EE[k]["RB"]) for k in ("nea", "mid", "far")}
        Lp, nu, de, t_ = copy_fit(Fl["nea"][0], src)
        b2 = nu / de[:, None]
        j2 = b2.argmin(1)
        bneg = b2[i_, j2]
        jabs = np.abs(b2).argmax(1)
        babs = np.abs(b2)[i_, jabs]
        mneg = ok & (bneg < 0)
        la_neg = np.log(np.sqrt(t_[i_, j2] / de))
        p = dict(source_ms=list(src), ref_ms=[list(w) for w in refw], far_near_all=r4(gD(iv, EE=EE)),
                 n_negative=int(mneg.sum()),
                 coef_global_min=r4((nu[ok].sum(0) / de[ok].sum()).min()), coef_global_max=r4((nu[ok].sum(0) / de[ok].sum()).max()),
                 lnD_vs_ln_coef_most_negative=corr_block(lD, np.where(mneg, np.log(np.where(mneg, -bneg, 1.0)), np.nan), mneg, geo, full),
                 lnD_vs_ln_abs_coef=corr_block(lD, np.log(babs), ok, geo, full),
                 lnD_vs_ln_amp_ratio_no_shape=corr_block(lD, la_neg, ok, geo, full),
                 lnD_refbelow_vs_near_level=corr_block(lDB, lv["nea"], ok, geo, False),
                 lnD_refbelow_vs_far_level=corr_block(lDB, lv["far"], ok, geo, False),
                 near_level_vs_ln_abs_coef=corr_block(lv["nea"], np.log(babs), ok, geo, False),
                 mid_level_vs_ln_abs_coef=corr_block(lv["mid"], np.log(babs), ok, geo, False),
                 far_level_vs_ln_abs_coef=corr_block(lv["far"], np.log(babs), ok, geo, False),
                 groups_by_abs_coef=groups_by(np.log(babs), ok, EE=EE),
                 groups_by_amp_ratio_no_shape=groups_by(la_neg, ok, EE=EE))
        PL[f"{T[0]}-{T[1]}"] = p
    res["placebo"] = PL

    # ---- predictor that shares no samples with the zone or with the reference of the ratio: package level on near
    def indep(T, S, X, Y):
        e = lambda k, w: energy(Fl[k], w, "rms")
        pk = 0.5 * np.log(e("nea", [S]) / e("nea", X))
        pkf = 0.5 * np.log(e("far", [S]) / e("far", X))
        lD = 0.5 * (np.log(e("far", [T]) / e("nea", [T])) - np.log(e("far", Y) / e("nea", Y)))
        ln_ = 0.5 * np.log(e("nea", [T]) / e("nea", Y))
        lf_ = 0.5 * np.log(e("far", [T]) / e("far", Y))
        return dict(target=list(T), source=list(S), norm_predictor=[list(w) for w in X], ref_ratio=[list(w) for w in Y],
                    lnD_vs_source_level_near=corr_block(lD, pk, ok, geo, full),
                    near_target_level_vs_near_source_level=corr_block(ln_, pk, ok, geo, full),
                    far_target_level_vs_far_source_level=corr_block(lf_, pkf, ok, geo, full),
                    far_target_level_vs_near_source_level=corr_block(lf_, pk, ok, geo, False))
    res["independent_predictor"] = dict(
        zone=indep(Z, SRC, [(1700, 1880)], [(1150, 1240), (1580, 1690)]),
        placebo_1580_1680=indep((1580, 1680), (1436, 1536), [(1150, 1250)], [(1256, 1356), (1780, 1880)]),
        placebo_1680_1780=indep((1680, 1780), (1536, 1636), [(1150, 1250)], [(1256, 1356), (1800, 1880)]),
        placebo_1780_1880=indep((1780, 1880), (1636, 1736), [(1150, 1250)], [(1256, 1356), (1580, 1630)]))

    # ---- Test 3: where does the block-to-block variation of far/near come from, zone against other windows
    T3 = {}
    for nm, T in (("1150-1250", (1150, 1250)), ("1256-1356", (1256, 1356)), ("1400-1500 zone", Z), ("1580-1680", (1580, 1680)), ("1680-1780", (1680, 1780)), ("1780-1880", (1780, 1880))):
        refw = [w for w in ((1150, 1360), (1580, 1880)) if not (w[0] <= T[0] and T[1] <= w[1])]
        if len(refw) == 2 and nm != "1400-1500 zone":
            pass
        if nm in ("1150-1250", "1256-1356"):
            refw = [RB]
        elif nm != "1400-1500 zone":
            refw = [RA]
        e = lambda k, w: energy(Fl[k], w, "rms")
        ln_ = 0.5 * np.log(e("nea", [T]) / e("nea", refw))[ok]
        lf_ = 0.5 * np.log(e("far", [T]) / e("far", refw))[ok]
        lm_ = 0.5 * np.log(e("mid", [T]) / e("mid", refw))[ok]
        lD = lf_ - ln_
        vD = np.var(lD)
        T3[nm] = dict(ref=[list(w) for w in refw], sd_ln_far_near=r4(np.sqrt(vD)), r_near_far=r4(pear(ln_, lf_)), r_near_mid=r4(pear(ln_, lm_)), r_mid_far=r4(pear(lm_, lf_)),
                      share_from_far=r4(np.cov(lD, lf_)[0, 1] / np.var(lD, ddof=1)), share_from_near=r4(-np.cov(lD, ln_)[0, 1] / np.var(lD, ddof=1)),
                      sd_ln_near=r4(np.std(ln_)), sd_ln_mid=r4(np.std(lm_)), sd_ln_far=r4(np.std(lf_)))
        o = np.argsort(lD)
        q = len(o) // 4
        ii = np.where(ok)[0]
        T3[nm]["deepest_over_shallowest_quarter_same_data"] = {k: r4(np.sqrt((e(k, [T])[ii[o[:q]]].sum() / e(k, refw)[ii[o[:q]]].sum()) / (e(k, [T])[ii[o[-q:]]].sum() / e(k, refw)[ii[o[-q:]]].sum()))) for k in ("nea", "mid", "far")}
    res["test3"] = T3

    # ---- Test 1: are the package and the event near 1386 ms parallel
    d1386, e1 = residual_shift(Fl["ful"][0], ok, (1370, 1402))
    d1276, e2 = residual_shift(Fl["ful"][0], ok, (1262, 1290))
    d1316, e3 = residual_shift(Fl["ful"][0], ok, (1302, 1330))
    m1 = ok & ~e1 & ~e3
    iso = (d1386 - d1316)[m1]
    res["test1"] = dict(n=int(m1.sum()), interval_1316_to_1386_sd_ms=r4(np.std(iso)), p5_p95_span_ms=r4(np.percentile(iso, 95) - np.percentile(iso, 5)),
                        robust_sd_ms=r4(1.4826 * np.median(np.abs(iso - np.median(iso)))), interval_1276_to_1316_sd_ms=r4(np.std((d1316 - d1276)[ok & ~e2 & ~e3])))

    if full:
        # ---- time profile of far/near and mid/near (100 ms window, 10 ms step), normalised by 1150-1360 + 1580-1880
        prof = []
        for tc in range(1200, 1840, 10):
            w = [(tc - 50, tc + 50)]
            row = dict(t=tc)
            for a, b in (("far", "nea"), ("mid", "nea")):
                v = (energy(Fl[a], w, "rms")[ok].sum() / energy(Fl[b], w, "rms")[ok].sum()) / (E["rms"][a]["R"][ok].sum() / E["rms"][b]["R"][ok].sum())
                row[f"{a}_{b}"] = r4(np.sqrt(v))
            for k in ("nea", "mid", "far"):
                row[f"level_{k}"] = r4(np.sqrt(energy(Fl[k], w, "rms")[ok].sum() / E["rms"][k]["R"][ok].sum() * (sum(si(x[1]) - si(x[0]) + 1 for x in (RA, RB)) / 51.0)))
            prof.append(row)
        res["time_profile"] = prof
        outside = [p["far_nea"] for p in prof if p["t"] <= 1330 or p["t"] >= 1610]
        inside = [p for p in prof if 1450 <= p["t"] <= 1450]
        res["time_profile_summary"] = dict(min_far_near=min(p["far_nea"] for p in prof), t_of_min=[p["t"] for p in prof if p["far_nea"] == min(q["far_nea"] for q in prof)][0],
                                           min_outside_1380_1560_windows=min(outside), max_outside=max(outside), median_outside=r4(np.median(outside)),
                                           at_1450=inside[0]["far_nea"], min_mid_near=min(p["mid_nea"] for p in prof),
                                           t_of_min_mid_near=[p["t"] for p in prof if p["mid_nea"] == min(q["mid_nea"] for q in prof)][0])
        # ---- time shift of far and mid against near (sum of block cross-correlations), 100 ms windows
        lagp = []
        Lr = np.arange(-12, 13)
        for tc in range(1220, 1840, 30):
            s = sl((tc - 50, tc + 50))
            row = dict(t=tc)
            for k in ("far", "mid"):
                N, Fk = Fl["nea"][0][ok], Fl[k][0][ok]
                xc = np.array([(N[:, s] * Fk[:, s.start + L:s.stop + L]).sum() for L in Lr])
                nrm = np.sqrt((N[:, s] ** 2).sum() * (Fk[:, s] ** 2).sum())
                j = int(np.clip(xc.argmax(), 1, len(Lr) - 2))
                d = 0.5 * (xc[j - 1] - xc[j + 1]) / (xc[j - 1] - 2 * xc[j] + xc[j + 1])
                row[f"{k}_after_near_ms"] = r4((Lr[j] + d) * DT)
                row[f"{k}_cc_peak"] = r4(xc[j] / nrm)
                row[f"{k}_cc_zero"] = r4(xc[12] / nrm)
            lagp.append(row)
        res["lag_profile"] = lagp
        # polarity and copy of the area mean near trace
        mn = Fl["nea"][0][ok].mean(0)
        res["area_mean_near_trace"] = {str(t): r4(mn[si(t)] / np.abs(mn).max()) for t in (1276, 1316, 1386, 1426, 1496, 1536)}
        j0, j1 = si(1250), si(1340)
        pk = j0 + int(np.argmax(np.abs(mn[j0:j1])))
        res["area_mean_near_trace"]["strongest_1250_1340_ms"] = [T0 + DT * pk, r4(mn[pk] / np.abs(mn).max())]
        np.savez(OUT + f"/maps_{proc}.npz", lnD=lnDb.reshape(NBI, NBX), lnb=lnb.reshape(NBI, NBX), lnr=lnr.reshape(NBI, NBX), lna=lna.reshape(NBI, NBX),
                 near=lvl["nea"].reshape(NBI, NBX), far=lvl["far"].reshape(NBI, NBX), ok=ok.reshape(NBI, NBX))
    return res


RES = {"setup": dict(zone_ms=[1400, 1500], ref_above_ms=[1150, 1360], ref_below_ms=[1580, 1880], copy_source_ms=[1256, 1356], copy_lag_search_ms=[96, 128],
                     estimator="rms of block mean trace ('rms') and product of two checkerboard halves ('xp')", flatten="cross-correlation of full stack 1236-1346 ms, datum = block of well 15/9-13",
                     bootstrap="1000 resamples of fields of 3 x 4 blocks (600 m x 900 m)")}
for proc in ("p07", "p10"):
    RES[proc] = analyse(proc, 2, 3, True)
    print(proc, "done", flush=True)
    RES[proc + "_blocks300m"] = analyse(proc, 3, 4, False)
top = {k: Cc[f"p07_{k}_top"] for k in ("nea", "ful")}
RES["time_origin_check"] = {k: dict(first_sample_above_5pct_of_max_ms=float(np.argmax(v > 0.05 * v.max()) * DT), t_of_max_first_400ms=float(np.argmax(v) * DT)) for k, v in top.items()}
json.dump(RES, open(OUT + "/results.json", "w"), indent=1)
print("written")
