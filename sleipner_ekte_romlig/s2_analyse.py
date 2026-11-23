#!/usr/bin/env python3
"""Spatial behaviour of the far/near deficit at 1400-1500 ms on Sleipner 94p07 and 94p10.

Input: blockstack.npz from s1_blockstack.py (mean trace per cell of 8 x 9 traces). Unit of analysis: blocks of
2 x 2 cells (16 x 18 traces, 200 m x 225 m). Each block is split in two halves (diagonal cells A, anti-diagonal
cells B); signal energy in a window is the cross product sum(A*B), which random noise does not inflate.
Amplitude ratio far/near = sqrt(E_far / E_near): no alignment between stacks is needed.
Writes results.json, maps.npz and two figures.
"""
import json
import numpy as np

OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig"
DT = 2.0
NS = 1001
NBI, NBX = 16, 26
KS = ("nea", "mid", "far")
ZD = np.load(OUT + "/blockstack.npz")
rng = np.random.default_rng(11)

# nominal windows (ms)
ZONE = (1400.0, 1500.0)
REF_A = (1150.0, 1360.0)
REF_B = (1580.0, 1880.0)
SRC = (-60.0, 40.0)             # copy source window relative to the 1316 ms reflector (1256-1356 ms)
LAGS = np.arange(46, 68)        # copy lag search 92-134 ms


def medfilt3(a):
    p = np.pad(a, 1, mode="edge")
    st = np.stack([p[i:i + a.shape[0], j:j + a.shape[1]] for i in range(3) for j in range(3)], 0)
    return np.nanmedian(st, 0)


def take(X, idx):
    return np.take_along_axis(X, np.clip(idx, 0, NS - 1), axis=2)


def coarse(F, N):
    Fw = F * N[..., None]

    def comb(mask):
        num = np.zeros((NBI, NBX, NS))
        den = np.zeros((NBI, NBX))
        for a in range(2):
            for b in range(2):
                if mask[a][b]:
                    num += Fw[a::2, b::2]
                    den += N[a::2, b::2]
        return num / np.maximum(den, 1)[..., None]

    A = comb([[1, 0], [0, 1]])
    B = comb([[0, 1], [1, 0]])
    T = comb([[1, 1], [1, 1]])
    valid = (N[0::2, 0::2] >= 36) & (N[1::2, 0::2] >= 36) & (N[0::2, 1::2] >= 36) & (N[1::2, 1::2] >= 36)
    return A, B, T, valid


def track(T, valid, T0, half=15, maxlag=20, iters=2):
    ic = int(round(T0 / DT))
    lags = np.arange(-maxlag, maxlag + 1)
    seg = lambda L: T[:, :, ic - half + L: ic + half + 1 + L]
    pilot = seg(0)[valid].mean(0)
    lag = ccb = None
    for _ in range(iters):
        pn = np.sqrt((pilot ** 2).sum())
        cc = np.stack([(seg(L) * pilot).sum(-1) / (np.sqrt((seg(L) ** 2).sum(-1)) * pn + 1e-30) for L in lags], -1)
        g = medfilt3(np.where(valid, lags[cc.argmax(-1)].astype(float), np.nan))
        g = np.where(np.isfinite(g), g, 0.0)
        m = np.where(np.abs(lags[None, None, :] - g[..., None]) <= 3, cc, -9.0)
        j = np.clip(m.argmax(-1), 1, len(lags) - 2)
        y0, y1, y2 = (np.take_along_axis(cc, (j + d)[..., None], 2)[..., 0] for d in (-1, 0, 1))
        den = y0 - 2 * y1 + y2
        d = np.clip(np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0), -1, 1)
        lag, ccb = lags[j] + d, y1
        idx = (ic - half + lags[j])[..., None] + np.arange(2 * half + 1)
        sel = valid & (ccb > 0.5)
        pilot = take(T, idx)[sel].mean(0)
    return T0 + lag * DT, ccb


def win_energy(X, Y, datum, w, d0):
    """sum of X*Y over samples in [w0, w1] ms, the window shifted by (datum - d0)."""
    cs = np.concatenate([np.zeros(X.shape[:2] + (1,)), np.cumsum(X * Y, -1)], -1)
    i0 = np.rint((w[0] + datum - d0) / DT).astype(int)
    i1 = np.rint((w[1] + datum - d0) / DT).astype(int) + 1
    i0, i1 = np.clip(i0, 0, NS), np.clip(i1, 0, NS)
    return np.take_along_axis(cs, i1[..., None], 2)[..., 0] - np.take_along_axis(cs, i0[..., None], 2)[..., 0]


def copy_strength(A, B, T, ts, fixed_lag=None):
    i0 = np.rint((ts + SRC[0]) / DT).astype(int)
    n = int((SRC[1] - SRC[0]) / DT) + 1
    idx = i0[..., None] + np.arange(n)
    sA, sB, sT = take(A, idx), take(B, idx), take(T, idx)
    den = (sA * sB).sum(-1)
    denT = (sT ** 2).sum(-1)
    num = np.empty((NBI, NBX, len(LAGS)))
    rho = np.empty((NBI, NBX, len(LAGS)))
    for j, L in enumerate(LAGS):
        tA, tB, tT = take(A, idx + L), take(B, idx + L), take(T, idx + L)
        num[..., j] = 0.5 * ((tA * sB).sum(-1) + (tB * sA).sum(-1))
        rho[..., j] = (tT * sT).sum(-1) / np.sqrt((tT ** 2).sum(-1) * denT + 1e-30)
    c = num / den[..., None]
    if fixed_lag is None:
        g = medfilt3(LAGS[rho.argmin(-1)].astype(float))
        m = np.where(np.abs(LAGS[None, None, :] - g[..., None]) <= 2, rho, 9.0)
        j = np.clip(m.argmin(-1), 1, len(LAGS) - 2)
    else:
        j = np.full((NBI, NBX), int(np.argmin(np.abs(LAGS - fixed_lag))))
    y0, y1, y2 = (np.take_along_axis(rho, (j + d)[..., None], 2)[..., 0] for d in (-1, 0, 1))
    dd = y0 - 2 * y1 + y2
    d = np.clip(np.where(np.abs(dd) > 1e-12, 0.5 * (y0 - y2) / dd, 0.0), -1, 1)
    lag_ms = (LAGS[j] + d) * DT
    cj = np.take_along_axis(c, j[..., None], 2)[..., 0]
    rj = y1
    Lj = LAGS[j][..., None]
    cA = (take(A, idx + Lj) * sA).sum(-1) / (sA ** 2).sum(-1)
    cB = (take(B, idx + Lj) * sB).sum(-1) / (sB ** 2).sum(-1)
    return dict(c=cj, rho=rj, lag=lag_ms, cA=cA, cB=cB, num=num, den=den)


def pearson(x, y):
    x, y = x - x.mean(), y - y.mean()
    return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum() + 1e-300))


def rank(x):
    r = np.empty(len(x))
    r[np.argsort(x)] = np.arange(len(x))
    return r


def spearman(x, y):
    return pearson(rank(x), rank(y))


def slope(x, y):
    x = x - x.mean()
    return float((x * (y - y.mean())).sum() / (x * x).sum())


def sb_groups(ok, sb=4):
    ids = (np.arange(NBI)[:, None] // sb) * 100 + (np.arange(NBX)[None, :] // sb)
    return [np.flatnonzero((ids == g) & ok) for g in np.unique(ids[ok])]


def sb_boot(fn, ok, nrep=2000, sb=4):
    """Bootstrap over super-blocks of sb x sb blocks (800 m x 900 m for sb=4). fn takes flat block indices."""
    G = sb_groups(ok, sb)
    out = []
    for _ in range(nrep):
        pick = rng.integers(0, len(G), len(G))
        out.append(fn(np.concatenate([G[p] for p in pick])))
    out = np.array(out, float)
    return [float(v) for v in np.nanpercentile(out, [2.5, 97.5], axis=0).T.ravel()] if out.ndim > 1 else [float(v) for v in np.nanpercentile(out, [2.5, 97.5])]


def lattice(x, y, ok, step=3):
    rs, ns = [], []
    for a in range(step):
        for b in range(step):
            m = np.zeros((NBI, NBX), bool)
            m[a::step, b::step] = True
            m &= ok
            if m.sum() >= 8:
                rs.append(pearson(x[m], y[m]))
                ns.append(int(m.sum()))
    return dict(mean_r=float(np.mean(rs)), min_r=float(np.min(rs)), max_r=float(np.max(rs)), n_per_lattice=float(np.mean(ns)),
                approx_se=float(1 / np.sqrt(np.mean(ns) - 3)))


def r3(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 4)


def analyse(proc):
    res = {}
    S = {}
    for k in KS + ("ful",):
        A, B, T, valid = coarse(ZD[f"{proc}_{k}"].astype(float), ZD[f"{proc}_{k}_n"].astype(float))
        S[k] = (A, B, T)
    res["n_blocks_valid"] = int(valid.sum())
    res["area_km2_valid"] = round(float(valid.sum() * 0.2 * 0.225), 2)

    # ---- horizons on the full stack
    hz, cc = {}, {}
    for T0 in (1276, 1316, 1386, 1456):
        tt, c = track(S["ful"][2], valid, T0)
        hz[T0] = medfilt3(np.where(valid, tt, np.nan))
        cc[T0] = c
        hz[T0] = np.where(np.isfinite(hz[T0]), hz[T0], np.nanmedian(tt[valid]))
        v = tt[valid]
        res[f"horizon_{T0}"] = dict(median_ms=r3(np.median(v)), p5=r3(np.percentile(v, 5)), p95=r3(np.percentile(v, 95)),
                                    min=r3(hz[T0][valid].min()), max=r3(hz[T0][valid].max()), cc_median=r3(np.median(c[valid])),
                                    frac_cc_above_0p6=r3(np.mean(c[valid] > 0.6)))
    mn = S["nea"][2][valid].mean(0)
    isb = int(np.argmax(np.abs(mn[40:90]))) + 40
    tsb, csb = track(S["nea"][2], valid, isb * DT, half=6, maxlag=8)
    hz["sb"] = np.where(valid, tsb, np.nan)
    res["seabed_near"] = dict(median_ms=r3(np.median(tsb[valid])), p5=r3(np.percentile(tsb[valid], 5)), p95=r3(np.percentile(tsb[valid], 95)),
                              cc_median=r3(np.median(csb[valid])))
    d1316, d1386 = float(np.median(hz[1316][valid])), float(np.median(hz[1386][valid]))
    hist = np.histogram((hz[1456] - hz[1316])[valid], bins=np.arange(100, 172, 6))[0]
    res["pick_1456_minus_1316_histogram_6ms_bins_from_100"] = [int(v) for v in hist]

    # ---- Test 1 input: isochrons
    iso = {}
    for a, b in ((1316, 1456), (1316, 1386), (1276, 1456), (1276, 1386), (1276, 1316), (1386, 1456)):
        I = hz[b] - hz[a]
        good = valid & (cc[a] > 0.6) & (cc[b] > 0.6)
        v = I[good]
        iso[f"{a}_{b}"] = dict(n=int(good.sum()), median=r3(np.median(v)), p5=r3(np.percentile(v, 5)), p95=r3(np.percentile(v, 95)),
                               p95_minus_p5=r3(np.percentile(v, 95) - np.percentile(v, 5)), sd=r3(v.std()), min=r3(v.min()), max=r3(v.max()))
    res["isochron_ms"] = iso
    res["structure_relief_1316_ms_p95_minus_p5"] = r3(np.percentile(hz[1316][valid], 95) - np.percentile(hz[1316][valid], 5))

    # ---- copy on near (and mid, far for information)
    cp = {k: copy_strength(*S[k], hz[1316]) for k in KS}
    gl = {}
    for k in KS:
        cg = cp[k]["num"][valid].sum(0) / cp[k]["den"][valid].sum()
        j = int(np.argmin(cg))
        gl[k] = dict(coef=r3(cg[j]), lag_ms=float(LAGS[j] * DT))
    res["copy_global"] = gl
    cn = cp["nea"]
    okc = valid & (cn["den"] > 0) & (cn["c"] < 0)
    res["copy_near_blocks"] = dict(n=int(okc.sum()), coef_median=r3(np.median(cn["c"][okc])), coef_p5=r3(np.percentile(cn["c"][okc], 5)),
                                   coef_p95=r3(np.percentile(cn["c"][okc], 95)), coef_sd=r3(cn["c"][okc].std()),
                                   rho_median=r3(np.median(cn["rho"][okc])), lag_median_ms=r3(np.median(cn["lag"][okc])),
                                   lag_p5=r3(np.percentile(cn["lag"][okc], 5)), lag_p95=r3(np.percentile(cn["lag"][okc], 95)),
                                   halves_reliability_r=r3(pearson(cn["cA"][okc], cn["cB"][okc])))
    rel = res["copy_near_blocks"]["halves_reliability_r"]
    res["copy_near_blocks"]["reliability_full_spearman_brown"] = r3(2 * rel / (1 + rel))
    cfix = copy_strength(*S["nea"], hz[1316], fixed_lag=gl["nea"]["lag_ms"] / DT)

    # ---- energies
    def energies(X, Y):
        E = {}
        for k in KS:
            x, y = S[k][X], S[k][Y]
            E[k] = dict(Z=win_energy(x, y, hz[1386], ZONE, d1386), Zs=win_energy(x, y, hz[1316], ZONE, d1316),
                        Ra=win_energy(x, y, hz[1316], REF_A, d1316), Rb=win_energy(x, y, hz[1386], REF_B, d1386))
        return E

    E = energies(0, 1)          # cross halves: noise-unbiased
    EA, EB = energies(0, 0), energies(1, 1)

    def lnL(E, k, zone="Z", ref="ab"):
        R = (E[k]["Ra"] if "a" in ref else 0) + (E[k]["Rb"] if "b" in ref else 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return 0.5 * (np.log(E[k][zone]) - np.log(R))

    def lnD(E, a="far", zone="Z", ref="ab"):
        return lnL(E, a, zone, ref) - lnL(E, "nea", zone, ref)

    def gD(E, idx, a="far", zone="Z", ref="ab"):
        def s(k, key):
            return E[k][key].ravel()[idx].sum()
        R = lambda k: (s(k, "Ra") if "a" in ref else 0) + (s(k, "Rb") if "b" in ref else 0)
        return float(np.sqrt(s(a, zone) / s("nea", zone)) / np.sqrt(R(a) / R("nea")))

    def gL(E, idx, k, zone="Z", ref="ab"):
        R = (E[k]["Ra"].ravel()[idx].sum() if "a" in ref else 0) + (E[k]["Rb"].ravel()[idx].sum() if "b" in ref else 0)
        return float(np.sqrt(E[k][zone].ravel()[idx].sum() / R))

    lD = lnD(E)
    lDm = lnD(E, "mid")
    ok = valid & np.isfinite(lD) & np.isfinite(lDm) & okc
    for k in KS:
        ok &= np.isfinite(lnL(E, k)) & np.isfinite(lnL(E, k, ref="b")) & np.isfinite(lnL(E, k, zone="Zs"))
        ok &= np.isfinite(lnL(EA, k)) & np.isfinite(lnL(EB, k)) & np.isfinite(lnL(EA, k, ref="b")) & np.isfinite(lnL(EB, k, ref="b"))
    idx_all = np.flatnonzero(ok)
    res["n_blocks_used"] = int(ok.sum())

    G = {}
    for name, a, zone, ref in (("far_near_zone_follow_zone", "far", "Z", "ab"), ("mid_near_zone_follow_zone", "mid", "Z", "ab"),
                               ("far_near_zone_follow_strong", "far", "Zs", "ab"), ("mid_near_zone_follow_strong", "mid", "Zs", "ab"),
                               ("far_near_ref_above_only", "far", "Z", "a"), ("far_near_ref_below_only", "far", "Z", "b"),
                               ("mid_near_ref_above_only", "mid", "Z", "a"), ("mid_near_ref_below_only", "mid", "Z", "b")):
        ci = sb_boot(lambda i: gD(E, i, a, zone, ref), ok, 1000)
        G[name] = dict(zone_over_ref=r3(gD(E, idx_all, a, zone, ref)), ci95=[r3(ci[0]), r3(ci[1])])
    for k in KS:
        ci = sb_boot(lambda i: gL(E, i, k), ok, 1000)
        G[f"level_{k}_zone_over_ref"] = dict(value=r3(gL(E, idx_all, k)), ci95=[r3(ci[0]), r3(ci[1])])
        G[f"level_{k}_zone_over_ref_below"] = r3(gL(E, idx_all, k, ref="b"))
    # absolute ratios (unknown scaling between stacks)
    for a in ("far", "mid"):
        G[f"abs_{a}_near_zone"] = r3(np.sqrt(E[a]["Z"].ravel()[idx_all].sum() / E["nea"]["Z"].ravel()[idx_all].sum()))
        G[f"abs_{a}_near_ref"] = r3(np.sqrt((E[a]["Ra"] + E[a]["Rb"]).ravel()[idx_all].sum() / (E["nea"]["Ra"] + E["nea"]["Rb"]).ravel()[idx_all].sum()))
        G[f"abs_{a}_near_ref_above"] = r3(np.sqrt(E[a]["Ra"].ravel()[idx_all].sum() / E["nea"]["Ra"].ravel()[idx_all].sum()))
        G[f"abs_{a}_near_ref_below"] = r3(np.sqrt(E[a]["Rb"].ravel()[idx_all].sum() / E["nea"]["Rb"].ravel()[idx_all].sum()))
    res["global"] = G
    PW = {}
    for est, (X, Y) in (("split_half", (0, 1)), ("plain_rms_block_traces", (2, 2))):
        e = {}
        for k in KS:
            z = win_energy(S[k][X], S[k][Y], hz[1386], (1408.0, 1492.0), d1386).ravel()[idx_all].sum()
            r_ = (win_energy(S[k][X], S[k][Y], hz[1316], (1100.0, 1380.0), d1316) + win_energy(S[k][X], S[k][Y], hz[1386], (1560.0, 1900.0), d1386)).ravel()[idx_all].sum()
            e[k] = float(np.sqrt(z / r_))
        PW[est] = dict(near_alene=r3(e["nea"]), mid_alene=r3(e["mid"]), far_alene=r3(e["far"]), far_over_near=r3(e["far"] / e["nea"]), mid_over_near=r3(e["mid"] / e["nea"]))
    res["same_windows_as_prediksjoner_json"] = PW
    res["block_lnD_far_near"] = dict(median_D=r3(np.exp(np.median(lD[ok]))), p5_D=r3(np.exp(np.percentile(lD[ok], 5))),
                                     p95_D=r3(np.exp(np.percentile(lD[ok], 95))), sd_ln=r3(lD[ok].std()),
                                     frac_blocks_below_1=r3(np.mean(lD[ok] < 0)))

    # ---- Test 2: deficit against copy strength
    def test2(lDx, cabs, tag):
        x, y = np.log(cabs), lDx
        f = lambda i: [pearson(x.ravel()[i], y.ravel()[i]), spearman(x.ravel()[i], y.ravel()[i]), slope(x.ravel()[i], y.ravel()[i])]
        v = f(idx_all)
        ci = sb_boot(f, ok, 2000)
        ci2 = sb_boot(f, ok, 1000, sb=6)
        return dict(what=tag, n=int(ok.sum()), pearson_r=r3(v[0]), pearson_ci95_superblock4=[r3(ci[0]), r3(ci[1])],
                    pearson_ci95_superblock6=[r3(ci2[0]), r3(ci2[1])],
                    spearman_r=r3(v[1]), spearman_ci95=[r3(ci[2]), r3(ci[3])],
                    slope_lnD_per_ln_copy=r3(v[2]), slope_ci95=[r3(ci[4]), r3(ci[5])],
                    lattice_600m=lattice(x, y, ok, 3), lattice_800m=lattice(x, y, ok, 4))

    cabs = np.where(ok, -cn["c"], 1.0)
    cabs_fix = np.where(ok & (cfix["c"] < 0), -cfix["c"], np.nan)
    T2 = {}
    T2["far_near_vs_copy"] = test2(lD, cabs, "ln far/near zone/ref (zone follows zone reflectors, ref above+below) against ln |copy coefficient on near|")
    T2["far_near_refbelow_vs_copy"] = test2(lnD(E, ref="b"), cabs, "same, reference below only (does not contain the copy source)")
    T2["far_near_followstrong_vs_copy"] = test2(lnD(E, zone="Zs"), cabs, "same, zone window follows the strong reflector")
    T2["mid_near_vs_copy"] = test2(lDm, cabs, "ln mid/near zone/ref against ln |copy|")
    T2["near_level_vs_copy"] = test2(lnL(E, "nea"), cabs, "ln near zone/ref against ln |copy| (power check: multiples must give a positive relation)")
    T2["near_level_refbelow_vs_copy"] = test2(lnL(E, "nea", ref="b"), cabs, "ln near zone/ref(below) against ln |copy|")
    T2["far_level_vs_copy"] = test2(lnL(E, "far"), cabs, "ln far zone/ref against ln |copy|")
    T2["far_level_refbelow_vs_copy"] = test2(lnL(E, "far", ref="b"), cabs, "ln far zone/ref(below) against ln |copy|")
    T2["far_near_vs_rho"] = test2(lD, np.where(ok, np.maximum(-cn["rho"], 1e-3), 1.0), "ln far/near zone/ref against ln |copy correlation on near|")
    if np.isfinite(cabs_fix[ok]).all():
        T2["far_near_vs_copy_fixedlag"] = test2(lD, np.where(ok, cabs_fix, 1.0), "same as first, copy coefficient at one fixed lag")
    def zero_copy(i, ref="ab"):
        x, y = cabs.ravel()[i], lnD(E, ref=ref).ravel()[i]
        b = slope(x, y)
        return [float(np.exp(y.mean() - b * x.mean())), b]
    for ref in ("ab", "b"):
        v = zero_copy(idx_all, ref)
        ci = sb_boot(lambda i: zero_copy(i, ref), ok, 2000)
        T2[f"extrapolation_to_zero_copy_ref_{ref}"] = dict(model="ln(far/near zone/ref) = a + b*|copy coefficient|; value at zero copy = exp(a)",
                                                          D_at_zero_copy=r3(v[0]), D_ci95=[r3(ci[0]), r3(ci[1])], b=r3(v[1]), b_ci95=[r3(ci[2]), r3(ci[3])],
                                                          copy_range_in_data=[r3(np.percentile(cabs[ok], 5)), r3(np.percentile(cabs[ok], 95))])
    def contrasts(i):
        c = cabs.ravel()[i]
        q10, q33, q67 = np.percentile(c, [10, 33.33, 66.67])
        out = []
        for sel in (i[c <= q33], i[c >= q67], i[c <= q10]):
            out += [gD(E, sel), gD(E, sel, ref="b"), gD(E, sel, "mid"), gL(E, sel, "nea", ref="b"), gL(E, sel, "mid", ref="b"), gL(E, sel, "far", ref="b"), gD(E, sel) / gD(E, sel, "mid")]
        return out
    v = contrasts(idx_all)
    ci = sb_boot(contrasts, ok, 2000)
    names = ("D_far_near", "D_far_near_refbelow", "D_mid_near", "level_near_refbelow", "level_mid_refbelow", "level_far_refbelow", "D_far_mid")
    NN = len(names)
    T2["copy_strength_groups_with_ci"] = {g: {n: dict(value=r3(v[NN * a + b]), ci95=[r3(ci[2 * (NN * a + b)]), r3(ci[2 * (NN * a + b) + 1])]) for b, n in enumerate(names)}
                                          for a, g in enumerate(("weak_tercile", "strong_tercile", "weakest_decile"))}
    def split_sel(i):
        out = []
        for cc_, Ee in ((cn["cA"], EB), (cn["cB"], EA)):
            c = -cc_.ravel()[i]
            q10, q33, q67 = np.percentile(c, [10, 33.33, 66.67])
            out.append([gD(Ee, i[c <= q10]), gD(Ee, i[c <= q33]), gD(Ee, i[c >= q67]), gD(Ee, i), gD(Ee, i[c <= q10], ref="b"),
                        gD(Ee, i[c <= q10]) / gD(Ee, i[c <= q10], "mid"), gD(Ee, i) / gD(Ee, i, "mid")])
        return [float(x) for x in np.exp(np.mean(np.log(np.array(out)), 0))]
    v = split_sel(idx_all)
    ci = sb_boot(split_sel, ok, 2000)
    T2["independent_halves_selection"] = dict(
        what="blocks grouped by the copy coefficient of one half of each block; far/near zone/ref measured on the other half (plain energies of that half); mean of the two ways",
        **{n: dict(value=r3(v[j]), ci95=[r3(ci[2 * j]), r3(ci[2 * j + 1])]) for j, n in enumerate(
            ("D_far_near_weakest_decile", "D_far_near_weak_tercile", "D_far_near_strong_tercile", "D_far_near_all_blocks", "D_far_near_refbelow_weakest_decile",
             "D_far_mid_weakest_decile", "D_far_mid_all_blocks"))})
    ci = sb_boot(lambda i: gD(E, i) / gD(E, i, "mid"), ok, 1000)
    T2["far_mid_all_blocks"] = dict(value=r3(gD(E, idx_all) / gD(E, idx_all, "mid")), ci95=[r3(ci[0]), r3(ci[1])])
    T2["far_mid_vs_copy"] = test2(lD - lDm, cabs, "ln far/mid zone/ref against ln |copy on near| (near's zone energy is not in this ratio)")
    T2["mid_level_refbelow_vs_copy"] = test2(lnL(E, "mid", ref="b"), cabs, "ln mid zone/ref(below) against ln |copy on near| (mid is independent data)")
    dq = np.percentile(cabs[ok], 10)
    wk = np.flatnonzero(ok & (cabs <= dq))
    T2["weakest_copy_decile"] = dict(n=int(len(wk)), copy_median=r3(np.median(cabs.ravel()[wk])), D_far_near=r3(gD(E, wk)), D_far_near_refbelow=r3(gD(E, wk, ref="b")),
                                     D_mid_near=r3(gD(E, wk, "mid")))
    # tercile contrast
    q1, q2 = np.percentile(cabs[ok], [33.33, 66.67])
    lo, hi = ok & (cabs <= q1), ok & (cabs >= q2)
    T2["tercile_contrast"] = dict(copy_weak_median=r3(np.median(cabs[lo])), copy_strong_median=r3(np.median(cabs[hi])),
                                  D_far_near_weak_copy=r3(gD(E, np.flatnonzero(lo))), D_far_near_strong_copy=r3(gD(E, np.flatnonzero(hi))),
                                  D_mid_near_weak_copy=r3(gD(E, np.flatnonzero(lo), "mid")), D_mid_near_strong_copy=r3(gD(E, np.flatnonzero(hi), "mid")),
                                  D_far_near_refbelow_weak_copy=r3(gD(E, np.flatnonzero(lo), ref="b")), D_far_near_refbelow_strong_copy=r3(gD(E, np.flatnonzero(hi), ref="b")),
                                  **{f"level_{k}_weak_copy": r3(gL(E, np.flatnonzero(lo), k)) for k in KS},
                                  **{f"level_{k}_strong_copy": r3(gL(E, np.flatnonzero(hi), k)) for k in KS},
                                  **{f"levelbelow_{k}_weak_copy": r3(gL(E, np.flatnonzero(lo), k, ref="b")) for k in KS},
                                  **{f"levelbelow_{k}_strong_copy": r3(gL(E, np.flatnonzero(hi), k, ref="b")) for k in KS})
    def zero_copy_indep(i):
        a_ = []
        for cc_, Ee in ((cn["cA"], EB), (cn["cB"], EA)):
            x, y = -cc_.ravel()[i], lnD(Ee).ravel()[i]
            b = slope(x, y)
            a_.append(y.mean() - b * x.mean())
        return float(np.exp(np.mean(a_)))
    st = rng.bit_generator.state
    ci = sb_boot(zero_copy_indep, ok, 2000)
    rng.bit_generator.state = st
    T2["extrapolation_to_zero_copy_independent_halves"] = dict(
        what="as extrapolation_to_zero_copy_ref_ab, but the copy coefficient comes from one half of each block and far/near from the other half",
        D_at_zero_copy=r3(zero_copy_indep(idx_all)), D_ci95=[r3(ci[0]), r3(ci[1])])
    res["test2"] = T2

    # ---- Test 3: is the deficit made by high near or by low far?
    T3 = {}
    for ref in ("ab", "b"):
        lDr = lnD(E, ref=ref)
        q25, q75 = np.percentile(lDr[ok], [25, 75])
        deep, shal = ok & (lDr <= q25), ok & (lDr >= q75)
        d = dict(D_deep=r3(gD(E, np.flatnonzero(deep), ref=ref)), D_shallow=r3(gD(E, np.flatnonzero(shal), ref=ref)))
        for k in KS:
            a, b = gL(E, np.flatnonzero(deep), k, ref=ref), gL(E, np.flatnonzero(shal), k, ref=ref)
            d[f"level_{k}_deep"], d[f"level_{k}_shallow"], d[f"level_{k}_deep_over_shallow"] = r3(a), r3(b), r3(a / b)
        T3[f"same_data_quartiles_ref_{ref}"] = d
        # split: rank on one half, evaluate on the other
        acc = {k: [] for k in KS}
        accD = []
        for Er, Ee in ((EA, EB), (EB, EA)):
            lr = lnD(Er, ref=ref)
            o2 = ok & np.isfinite(lr)
            q25, q75 = np.percentile(lr[o2], [25, 75])
            deep, shal = o2 & (lr <= q25), o2 & (lr >= q75)
            for k in KS:
                acc[k].append(np.log(gL(Ee, np.flatnonzero(deep), k, ref=ref) / gL(Ee, np.flatnonzero(shal), k, ref=ref)))
            accD.append(np.log(gD(Ee, np.flatnonzero(deep), ref=ref) / gD(Ee, np.flatnonzero(shal), ref=ref)))
        d2 = {f"level_{k}_deep_over_shallow": r3(np.exp(np.mean(acc[k]))) for k in KS}
        d2["D_deep_over_shallow"] = r3(np.exp(np.mean(accD)))
        T3[f"split_half_quartiles_ref_{ref}"] = d2

        # variance shares with independent halves
        def shares(i):
            nf = nn = dn = 0.0
            for Er, Ee in ((EA, EB), (EB, EA)):
                a = lnD(Er, ref=ref).ravel()[i]
                f_, n_ = lnL(Ee, "far", ref=ref).ravel()[i], lnL(Ee, "nea", ref=ref).ravel()[i]
                m = np.isfinite(a) & np.isfinite(f_) & np.isfinite(n_)
                a, f_, n_ = a[m] - a[m].mean(), f_[m], n_[m]
                nf += (a * (f_ - f_.mean())).sum()
                nn += -(a * (n_ - n_.mean())).sum()
                dn += (a * ((f_ - n_) - (f_ - n_).mean())).sum()
            return [nf / dn, nn / dn]
        v = shares(idx_all)
        ci = sb_boot(shares, ok, 2000)
        T3[f"variance_share_ref_{ref}"] = dict(far_low=r3(v[0]), far_ci95=[r3(ci[0]), r3(ci[1])], near_high=r3(v[1]), near_ci95=[r3(ci[2]), r3(ci[3])])
        T3[f"sd_ln_levels_ref_{ref}"] = {k: r3(lnL(E, k, ref=ref)[ok].std()) for k in KS}
        T3[f"corr_ln_near_far_levels_ref_{ref}"] = r3(pearson(lnL(E, "nea", ref=ref)[ok], lnL(E, "far", ref=ref)[ok]))
    res["test3"] = T3

    # ---- Test 1: does the far/near minimum follow the strong reflectors or the zone reflectors?
    I = hz[1386] - hz[1316]
    okI = ok & (cc[1316] > 0.6) & (cc[1386] > 0.6)
    tau = np.arange(-60, 151)                       # samples relative to the 1316 reflector: -120..+300 ms
    idx = np.rint(hz[1316] / DT).astype(int)[..., None] + tau
    P = {k: take(S[k][0], idx) * take(S[k][1], idx) for k in ("nea", "far", "mid")}
    box = np.ones(25) / 25

    def profile(sel, a="far"):
        ef = np.convolve(P[a].reshape(-1, len(tau))[sel].sum(0), box, mode="same")
        en = np.convolve(P["nea"].reshape(-1, len(tau))[sel].sum(0), box, mode="same")
        with np.errstate(invalid="ignore", divide="ignore"):
            return 0.5 * np.log(ef / en)

    def tmin(pr):
        w = (tau * DT >= 60) & (tau * DT <= 230)
        p = np.where(w & np.isfinite(pr), pr, np.inf)
        j = int(np.argmin(p))
        if 0 < j < len(tau) - 1 and np.isfinite(pr[j - 1]) and np.isfinite(pr[j + 1]):
            dd = pr[j - 1] - 2 * pr[j] + pr[j + 1]
            d = 0.5 * (pr[j - 1] - pr[j + 1]) / dd if abs(dd) > 1e-12 else 0.0
        else:
            d = 0.0
        return (tau[j] + float(np.clip(d, -1, 1))) * DT

    qs = np.percentile(I[okI], [0, 25, 50, 75, 100])
    grp = np.clip(np.searchsorted(qs, I, side="right") - 1, 0, 3)
    gi = [np.flatnonzero(okI & (grp == g)) for g in range(4)]

    def shift_of(pr, ref):
        w = np.flatnonzero((tau * DT >= 30) & (tau * DT <= 250))
        cs = np.array([pearson(pr[w + s], ref[w]) for s in range(-15, 16)])
        j = int(np.clip(np.argmax(cs), 1, len(cs) - 2))
        dd = cs[j - 1] - 2 * cs[j] + cs[j + 1]
        d = 0.5 * (cs[j - 1] - cs[j + 1]) / dd if abs(dd) > 1e-12 else 0.0
        return (j - 15 + float(np.clip(d, -1, 1))) * DT

    def t1(sel_idx):
        s = set(sel_idx.tolist())
        cnt = np.bincount(sel_idx, minlength=NBI * NBX)
        xs, ys = [], []
        for g in range(4):
            ii = np.repeat(gi[g], cnt[gi[g]])
            if len(ii) < 5:
                return [np.nan, np.nan]
            xs.append(I.ravel()[ii].mean())
            ys.append(shift_of(profile(ii), prof_ref))
        xs, ys = np.array(xs), np.array(ys)
        return [slope(xs, ys), ys.mean()]

    xs = [float(I.ravel()[g].mean()) for g in gi]
    prof_all = profile(np.flatnonzero(okI))
    prof_ref = prof_all
    ys = [float(shift_of(profile(g), prof_ref)) for g in gi]
    ci = sb_boot(t1, okI, 1000)
    res["test1"] = dict(isochron_1316_1386_group_means_ms=[r3(v) for v in xs], shift_of_far_near_profile_relative_to_1316_ms=[r3(v) for v in ys],
                        slope=r3(slope(np.array(xs), np.array(ys))), slope_ci95=[r3(ci[0]), r3(ci[1])],
                        slope_expected="0 if the minimum follows the strong reflectors (multiples), 1 if it follows the zone reflectors (clay)",
                        isochron_range_between_group_means_ms=r3(xs[-1] - xs[0]), n_blocks=int(okI.sum()),
                        tmin_all_blocks_after_1316_ms=r3(tmin(prof_all)))
    maps = dict(valid=valid, ok=ok, t1276=hz[1276], t1316=hz[1316], t1386=hz[1386], t1456=hz[1456], seabed=hz["sb"], iso=I,
                copy=np.where(ok, cn["c"], np.nan), copy_lag=np.where(ok, cn["lag"], np.nan), copy_rho=np.where(ok, cn["rho"], np.nan),
                D_fn=np.where(ok, np.exp(lD), np.nan), D_mn=np.where(ok, np.exp(lDm), np.nan),
                D_fn_strong=np.where(ok, np.exp(lnD(E, zone="Zs")), np.nan), D_fn_refbelow=np.where(ok, np.exp(lnD(E, ref="b")), np.nan),
                L_nea=np.where(ok, np.exp(lnL(E, "nea")), np.nan), L_mid=np.where(ok, np.exp(lnL(E, "mid")), np.nan),
                L_far=np.where(ok, np.exp(lnL(E, "far")), np.nan), prof_tau_ms=tau * DT, prof_all=prof_all,
                prof_groups=np.array([profile(g) for g in gi]), prof_mid_all=profile(np.flatnonzero(okI), "mid"))
    # ---- extra 1: time lag between near and far / mid (block traces flattened on 1316), per 100 ms window
    tau2 = np.arange(-90, 291)
    idx2 = np.rint(hz[1316] / DT).astype(int)[..., None] + tau2
    Tn, Tm, Tf = (take(S[k][2], idx2)[ok] for k in ("nea", "mid", "far"))
    ell = np.arange(-12, 13)
    lagprof = []
    for c_ms in (-100, -40, 0, 60, 90, 110, 130, 160, 200, 260, 320, 380, 440, 500):
        jc = int(np.argmin(np.abs(tau2 * DT - c_ms)))
        sl = slice(jc - 25, jc + 26)
        row = dict(centre_ms_after_1316=c_ms, approx_time_ms=r3(d1316 + c_ms))
        for nm, X in (("far", Tf), ("mid", Tm)):
            xcv = np.array([(Tn[:, sl] * X[:, jc - 25 + l: jc + 26 + l]).sum() for l in ell])
            nrm = np.sqrt((Tn[:, sl] ** 2).sum() * (X[:, sl] ** 2).sum())
            j = int(np.clip(np.argmax(xcv), 1, len(ell) - 2))
            dd = xcv[j - 1] - 2 * xcv[j] + xcv[j + 1]
            d = 0.5 * (xcv[j - 1] - xcv[j + 1]) / dd if abs(dd) > 1e-12 else 0.0
            row[f"lag_{nm}_after_near_ms"] = r3((ell[j] + float(np.clip(d, -1, 1))) * DT)
            row[f"cc_near_{nm}_best"] = r3(xcv[j] / nrm)
            row[f"cc_near_{nm}_zero_lag"] = r3(xcv[12] / nrm)
        lagprof.append(row)
    res["lag_profile"] = lagprof

    # ---- extra 2: do near and far levels vary together between blocks? zone against control windows of 100 ms
    def wE(k, w, datum, d0):
        return win_energy(S[k][0], S[k][1], datum, w, d0)
    ctrl = {}
    for name, w, dat, d0, inside in (("1150-1250", (1150.0, 1250.0), hz[1316], d1316, "a"), ("1256-1356 strong package", (1256.0, 1356.0), hz[1316], d1316, "a"),
                                     ("1400-1500 zone", ZONE, hz[1386], d1386, ""), ("1580-1680", (1580.0, 1680.0), hz[1386], d1386, "b"),
                                     ("1680-1780", (1680.0, 1780.0), hz[1386], d1386, "b"), ("1780-1880", (1780.0, 1880.0), hz[1386], d1386, "b")):
        lv = {}
        for k in KS:
            ew = wE(k, w, dat, d0)
            ref = E[k]["Ra"] + E[k]["Rb"] - (ew if inside else 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                lv[k] = 0.5 * (np.log(ew) - np.log(ref))
        m = ok & np.isfinite(lv["nea"]) & np.isfinite(lv["far"]) & np.isfinite(lv["mid"])
        f = lambda i: [pearson(lv["nea"].ravel()[i], lv["far"].ravel()[i]), pearson(lv["nea"].ravel()[i], lv["mid"].ravel()[i]), pearson(lv["mid"].ravel()[i], lv["far"].ravel()[i])]
        v = f(np.flatnonzero(m))
        ci = sb_boot(f, m, 1000)
        ctrl[name] = dict(n=int(m.sum()), r_near_far=r3(v[0]), r_near_far_ci95=[r3(ci[0]), r3(ci[1])], r_near_mid=r3(v[1]), r_near_mid_ci95=[r3(ci[2]), r3(ci[3])],
                          r_mid_far=r3(v[2]), r_mid_far_ci95=[r3(ci[4]), r3(ci[5])],
                          sd_ln_level={k: r3(lv[k][m].std()) for k in KS})
    res["level_covariation_between_blocks"] = ctrl
    # copy lag against sea-bed time (peg-leg check)
    sbok = ok & np.isfinite(hz["sb"])
    ci_sb = sb_boot(lambda i: [pearson(hz["sb"].ravel()[i], cn["lag"].ravel()[i]), slope(hz["sb"].ravel()[i], cn["lag"].ravel()[i])], sbok, 2000)
    iso2 = hz[1386] - hz[1276]
    res["copy_lag_vs_seabed"] = dict(r_ci95=[r3(ci_sb[0]), r3(ci_sb[1])], slope_ci95=[r3(ci_sb[2]), r3(ci_sb[3])], r_isochron_1276_1386_vs_seabed=r3(pearson(hz["sb"][sbok], iso2[sbok])), r=r3(pearson(hz["sb"][sbok], cn["lag"][sbok])), slope_ms_per_ms=r3(slope(hz["sb"][sbok], cn["lag"][sbok])),
                                     seabed_sd_ms=r3(hz["sb"][sbok].std()), lag_sd_ms=r3(cn["lag"][sbok].std()))
    return res, maps


RES, MAPS = {}, {}
for proc in ("p07", "p10"):
    RES[proc], m = analyse(proc)
    for k, v in m.items():
        MAPS[f"{proc}_{k}"] = v
# agreement between the processings, block by block
ok = MAPS["p07_ok"] & MAPS["p10_ok"]
RES["p07_vs_p10"] = dict(n=int(ok.sum()), r_ln_D_far_near=r3(pearson(np.log(MAPS["p07_D_fn"][ok]), np.log(MAPS["p10_D_fn"][ok]))),
                         r_copy=r3(pearson(MAPS["p07_copy"][ok], MAPS["p10_copy"][ok])))
RES["setup"] = dict(block_m=[200, 225], zone_window_ms=ZONE, ref_above_ms=REF_A, ref_below_ms=REF_B, copy_source_rel_1316_ms=SRC,
                    copy_lag_search_ms=[float(LAGS[0] * DT), float(LAGS[-1] * DT)], bootstrap="super-blocks of 4 x 4 blocks (800 m x 900 m), resampled with replacement",
                    estimator="amplitude ratio = sqrt of ratio of split-half cross energies sum(A*B); zone/ref = ratio in zone divided by ratio in reference")
json.dump(RES, open(OUT + "/results.json", "w"), indent=1)
np.savez(OUT + "/maps.npz", **MAPS)
print(json.dumps(RES, indent=1))
