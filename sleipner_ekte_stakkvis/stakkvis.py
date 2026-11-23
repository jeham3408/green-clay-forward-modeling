#!/usr/bin/env python3
"""stakkvis.py: amplitude level of near, mid and far, each stack on its own, in the zone 1408-1492 ms relative to
the reference (1100-1900 ms outside 1380-1560 ms), on the open Sleipner 1994 stacks (94p07 and 94p10) and on the
demultipled cubes. Read-only on the cubes (memory-mapped, 8 inlines at a time). Writes results.json and profiles.npz.

Run:  run_capped.sh 700 <geoscience python> -B stakkvis.py [version ...]
"""
import json
import os
import sys
import time
import numpy as np
from scipy.signal import hilbert

SL = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis"
VERS = {
    "p07": (SL + "/avo/{}.npy", SL + "/avo/ful.npy"),
    "p10": (SL + "/avo10/{}.npy", SL + "/avo10/ful.npy"),
    "p07_dmE": (SL + "/avo/dmE_{}.npy", SL + "/avo/ful.npy"),
    "p10_dmE": (SL + "/avo10/dmE_{}.npy", SL + "/avo10/ful.npy"),
    "p07_dm": (SL + "/avo/dm_{}.npy", SL + "/avo/ful.npy"),
    "p07_dmD": (SL + "/avo/dmD_{}.npy", SL + "/avo/ful.npy"),
    "p07_dmA": (SL + "/avo/dmA_{}.npy", SL + "/avo/ful.npy"),
    "p07_dmB": (SL + "/avo/dmB_{}.npy", SL + "/avo/ful.npy"),
}
KS = ("nea", "mid", "far")
NIL, NXL, NS = 249, 468, 1001
S0, S1 = 400, 1000                       # samples 400..999 = 800..1998 ms (sample i is at 2*i ms)
NT = S1 - S0
tms = 2.0 * np.arange(S0, S1)
ZM = (tms >= 1408) & (tms <= 1492)
RM = ((tms >= 1100) & (tms < 1380)) | ((tms > 1560) & (tms <= 1900))
NZ, NR = int(ZM.sum()), int(RM.sum())
BI, BX = 8, 9                            # small blocks for block-mean traces: 100 m x 112.5 m
NBI, NBX = (NIL + BI - 1) // BI, NXL // BX
WELL_IL, WELL_XL = 57, 143
# trace masks
il = np.arange(NIL)[:, None] + np.zeros((1, NXL), int)
xl = np.arange(NXL)[None, :] + np.zeros((NIL, 1), int)
AREA = {"kube": np.ones((NIL, NXL), bool),
        "bronn": (il <= WELL_IL + 80) & (xl >= WELL_XL - 80) & (xl <= WELL_XL + 80)}
SUB = (il // 83) * 6 + xl // 78          # 18 sub-areas of 83 x 78 traces = 1.04 km x 0.98 km
NSUB = 18
bil = (np.arange(NBI) * BI + (BI - 1) / 2)[:, None] + np.zeros((1, NBX))
bxl = (np.arange(NBX) * BX + (BX - 1) / 2)[None, :] + np.zeros((NBI, 1))
BSUB = (np.minimum(bil, NIL - 1) // 83).astype(int) * 6 + (bxl // 78).astype(int)
BAREA = {"kube": np.ones((NBI, NBX), bool),
         "bronn": (bil <= WELL_IL + 80) & (bxl >= WELL_XL - 80) & (bxl <= WELL_XL + 80)}


def st(v):
    v = np.asarray([x for x in v if x is not None and np.isfinite(x)], float)
    if len(v) < 2:
        return None
    med = float(np.median(v))
    sd = float(v.std(ddof=1))
    return dict(n=int(len(v)), median=med, mean=float(v.mean()), sd=sd, se=sd / np.sqrt(len(v)),
                rsd=float(1.4826 * np.median(np.abs(v - med))), min=float(v.min()), max=float(v.max()))


def rnd(o):
    if isinstance(o, dict):
        return {str(k): rnd(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [rnd(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 4)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return rnd(o.tolist())
    return o


def sdiv(a, b):
    return float(a / b) if (b != 0 and np.isfinite(a) and np.isfinite(b)) else float("nan")


def ssqrt(x):
    return float(np.sqrt(x)) if x > 0 else float("nan")


def levels(T, k, m):
    """Zone/reference amplitude level of stack k over the traces in mask m: pooled and trace by trace."""
    t = T[k]
    pooled = dict(rms=ssqrt(sdiv(t["ssz"][m].sum() / NZ, t["ssr"][m].sum() / NR)),
                  abs=sdiv(t["abz"][m].sum(), t["abr"][m].sum()),
                  env=sdiv(t["enz"][m].sum(), t["enr"][m].sum()),
                  coh=ssqrt(sdiv(t["coz"][m].sum() / NZ, t["cor"][m].sum() / NR)))
    return pooled


def per_trace(T, k, m):
    t = T[k]
    return dict(rms=np.sqrt((t["ssz"][m] / NZ) / (t["ssr"][m] / NR)), abs=t["abz"][m] / t["abr"][m],
                env=t["enz"][m] / t["enr"][m])


def pct(x):
    x = x[np.isfinite(x)]
    q = np.percentile(x, [16, 50, 84])
    return dict(median=float(q[1]), p16=float(q[0]), p84=float(q[2]), n=int(len(x)))


OFF = np.arange(-8, 9)                   # +-16 ms window around a pick
LAGS = np.arange(-4, 5)                  # alignment search +-8 ms
PLO, PHI = (1100 - 800) // 2, (1900 - 800) // 2
COLS = ("r", "c", "t", "nn", "mm", "ff", "nm", "nf", "mf", "en", "em", "ef", "Lm", "Lf", "amp")


def reflectors(B, E, src, bok):
    """Picks = local extrema of the block-mean trace of stack `src` (all above 0.3 x its rms). For each pick: windows of
    near, mid and far, mid and far aligned to near by the lag of maximum cross-correlation; products and envelope peaks."""
    recs = []
    i = np.arange(PLO, PHI + 1)
    for r in range(NBI):
        for c in range(NBX):
            if not bok[r, c]:
                continue
            p = B[src][r, c].astype(float)
            n, m, f = (B[k][r, c].astype(float) for k in KS)
            thr = 0.3 * np.sqrt((p[PLO:PHI + 1] ** 2).mean())
            a = np.abs(p)
            ext = (a[i] >= a[i - 1]) & (a[i] > a[i + 1]) & (a[i] > thr) & (p[i] * p[i - 1] > 0) & (p[i] * p[i + 1] > 0)
            idx = i[ext]
            if len(idx) == 0:
                continue
            w = idx[:, None] + OFF
            wn = n[w]

            def align(y):
                cc = np.stack([(wn * y[w + L]).sum(-1) for L in LAGS], 1)
                j = cc.argmax(1)
                return LAGS[j], y[w + LAGS[j][:, None]]
            Lm, wm = align(m)
            Lf, wf = align(f)
            e4 = idx[:, None] + np.arange(-4, 5)
            en, em, ef = (E[k][r, c][e4].max(1) for k in KS)
            recs.append(np.column_stack([np.full(len(idx), r), np.full(len(idx), c), tms[idx], (wn * wn).sum(-1),
                                         (wm * wm).sum(-1), (wf * wf).sum(-1), (wn * wm).sum(-1), (wn * wf).sum(-1),
                                         (wm * wf).sum(-1), en, em, ef, Lm, Lf, p[idx]]))
    return np.concatenate(recs)


def est(rec):
    """Amplitude ratios from a set of picks (rows of rec)."""
    if len(rec) < 5:
        return None
    g = {k: rec[:, j] for j, k in enumerate(COLS)}
    S = {k: g[k].sum() for k in ("nn", "mm", "ff", "nm", "nf", "mf", "en", "em", "ef")}
    return dict(n=int(len(rec)),
                mn_ls=sdiv(S["nm"], S["nn"]), mn_rms=ssqrt(sdiv(S["mm"], S["nn"])), mn_inv=sdiv(S["mm"], S["nm"]),
                mn_iv=sdiv(S["mf"], S["nf"]), mn_ls_med=float(np.median(g["nm"] / g["nn"])),
                mn_rms_med=float(np.median(np.sqrt(g["mm"] / g["nn"]))),
                mn_env_med=float(np.median(g["em"] / g["en"])), mn_env_pool=sdiv(S["em"], S["en"]),
                fn_ls=sdiv(S["nf"], S["nn"]), fn_rms=ssqrt(sdiv(S["ff"], S["nn"])), fn_inv=sdiv(S["ff"], S["nf"]),
                fn_iv=sdiv(S["mf"], S["nm"]), fn_ls_med=float(np.median(g["nf"] / g["nn"])),
                fn_rms_med=float(np.median(np.sqrt(g["ff"] / g["nn"]))),
                fn_env_med=float(np.median(g["ef"] / g["en"])), fn_env_pool=sdiv(S["ef"], S["en"]),
                c_nm=sdiv(S["nm"], np.sqrt(S["nn"] * S["mm"])), c_nf=sdiv(S["nf"], np.sqrt(S["nn"] * S["ff"])),
                c_mf=sdiv(S["mf"], np.sqrt(S["mm"] * S["ff"])),
                lag_mid_ms=float(np.median(g["Lm"]) * 2), lag_far_ms=float(np.median(g["Lf"]) * 2),
                lag_mid_mean_ms=float(np.mean(g["Lm"]) * 2), lag_far_mean_ms=float(np.mean(g["Lf"]) * 2),
                lag_mid_edge=float(np.mean(np.abs(g["Lm"]) == 4)), lag_far_edge=float(np.mean(np.abs(g["Lf"]) == 4)),
                pk_n=float(np.median(g["en"])), pk_m=float(np.median(g["em"])), pk_f=float(np.median(g["ef"])))


RATK = ("mn_ls", "mn_rms", "mn_inv", "mn_iv", "mn_ls_med", "mn_rms_med", "mn_env_med", "mn_env_pool",
        "fn_ls", "fn_rms", "fn_inv", "fn_iv", "fn_ls_med", "fn_rms_med", "fn_env_med", "fn_env_pool")


def zr(rec):
    t = rec[:, 2]
    z = est(rec[(t >= 1408) & (t <= 1492)])
    rf = est(rec[((t >= 1100) & (t < 1380)) | ((t > 1560) & (t <= 1900))])
    if z is None or rf is None:
        return None
    out = dict(zone=z, ref=rf, zr={k: sdiv(z[k], rf[k]) for k in RATK})
    out["zr"]["near_pk"] = sdiv(z["pk_n"], rf["pk_n"])
    out["zr"]["mid_pk"] = sdiv(z["pk_m"], rf["pk_m"])
    out["zr"]["far_pk"] = sdiv(z["pk_f"], rf["pk_f"])
    return out


def run(ver, keepB=None):
    t0 = time.time()
    pat, fulp = VERS[ver]
    paths = {k: pat.format(k) for k in KS}
    paths["ful"] = fulp
    for k in paths:
        assert np.load(paths[k], mmap_mode="r").shape == (NIL, NXL, NS), k
    names = ("ssz", "ssr", "abz", "abr", "enz", "enr", "coz", "cor")
    T = {k: {n: np.zeros((NIL, NXL), np.float64) for n in names} for k in KS}
    W100 = {k: np.zeros((NIL, NXL, 12), np.float32) for k in KS}
    E2t = {k: {a: np.zeros(NT) for a in AREA} for k in KS}
    B = {k: np.zeros((NBI, NBX, NT), np.float32) for k in KS + ("ful",)}
    for r in range(NBI):
        a, b = r * BI, min(NIL, (r + 1) * BI)
        for k in KS + ("ful",):
            mmk = np.load(paths[k], mmap_mode="r")          # reopened per chunk so mapped pages are released
            x = np.asarray(mmk[a:b, :, S0:S1], dtype=np.float64)
            del mmk
            B[k][r] = x.reshape(b - a, NBX, BX, NT).mean((0, 2))
            if k == "ful":
                continue
            e2 = x * x
            T[k]["ssz"][a:b] = e2[..., ZM].sum(-1)
            T[k]["ssr"][a:b] = e2[..., RM].sum(-1)
            ax = np.abs(x)
            T[k]["abz"][a:b] = ax[..., ZM].mean(-1)
            T[k]["abr"][a:b] = ax[..., RM].mean(-1)
            env = np.abs(hilbert(x, axis=-1))
            T[k]["enz"][a:b] = env[..., ZM].mean(-1)
            T[k]["enr"][a:b] = env[..., RM].mean(-1)
            W100[k][a:b] = e2.reshape(b - a, NXL, 12, 50).sum(-1)
            for an, am in AREA.items():
                E2t[k][an] += (e2 * am[a:b, :, None]).sum((0, 1))
            co = x[:, :-4, :] * x[:, 4:, :]
            T[k]["coz"][a:b, :-4] = co[..., ZM].sum(-1)
            T[k]["cor"][a:b, :-4] = co[..., RM].sum(-1)
            del x, e2, ax, env, co
    live = np.ones((NIL, NXL), bool)
    for k in KS:
        live &= (T[k]["ssr"] > 0) & (T[k]["ssz"] > 0)
    res = dict(live_fraction=float(live.mean()), n_zone_samples=NZ, n_ref_samples=NR)

    # 1. level of each stack in the zone relative to the reference; 2a. window ratios between the stacks
    res["nivaa"] = {}
    res["vindu_forhold"] = {}
    for an, am in AREA.items():
        m = am & live
        L = {k: levels(T, k, m) for k in KS}
        PT = {k: per_trace(T, k, m) for k in KS}
        res["nivaa"][an] = {k: dict(pooled=L[k], trace={e: pct(PT[k][e]) for e in ("rms", "abs", "env")}) for k in KS}
        res["nivaa"][an]["n_traces"] = int(m.sum())
        vf = {}
        for nm_, (k1, k2) in (("mid_near", ("mid", "nea")), ("far_near", ("far", "nea")), ("far_mid", ("far", "mid"))):
            vf[nm_] = dict(pooled={e: sdiv(L[k1][e], L[k2][e]) for e in ("rms", "abs", "env", "coh")},
                           trace={e: pct(PT[k1][e] / PT[k2][e]) for e in ("rms", "abs", "env")})
        res["vindu_forhold"][an] = vf
    subL = []
    for s in range(NSUB):
        m = (SUB == s) & live
        subL.append({k: levels(T, k, m) for k in KS} if m.sum() > 500 else None)
    res["nivaa"]["delomraader"] = {k: {e: st([d[k][e] for d in subL if d]) for e in ("rms", "abs", "env", "coh")} for k in KS}
    res["nivaa"]["delomraader_verdier"] = {k: {e: [d[k][e] if d else None for d in subL] for e in ("rms", "coh")} for k in KS}
    res["vindu_forhold"]["delomraader"] = {
        nm_: {e: st([sdiv(d[k1][e], d[k2][e]) for d in subL if d]) for e in ("rms", "abs", "env", "coh")}
        for nm_, (k1, k2) in (("mid_near", ("mid", "nea")), ("far_near", ("far", "nea")), ("far_mid", ("far", "mid")))}

    # 3. gain: rms against time in 100 ms windows, 800-2000 ms
    g = {}
    for an, am in AREA.items():
        m = am & live
        rms = {k: np.sqrt(W100[k][m].astype(np.float64).sum(0) / (m.sum() * 50)) for k in KS}
        g[an] = dict(t_start_ms=[800 + 100 * j for j in range(12)], **{k: rms[k] for k in KS},
                     mid_near=rms["mid"] / rms["nea"], far_near=rms["far"] / rms["nea"], far_mid=rms["far"] / rms["mid"])
    res["gain_100ms"] = g
    # the same per sub-area: far/near and mid/near of the window 1400-1500 relative to the median of the clean
    # reference windows (1100-1300 and 1600-1900)
    refw = [3, 4, 8, 9, 10]
    sv = {"mid_near": [], "far_near": [], "nea": [], "mid": [], "far": []}
    for s in range(NSUB):
        m = (SUB == s) & live
        if m.sum() <= 500:
            continue
        rms = {k: np.sqrt(W100[k][m].astype(np.float64).sum(0) / (m.sum() * 50)) for k in KS}
        for k in KS:
            sv[k].append(rms[k][6] / np.median(rms[k][refw]))
        for nm_, k1 in (("mid_near", "mid"), ("far_near", "far")):
            q = rms[k1] / rms["nea"]
            sv[nm_].append(q[6] / np.median(q[refw]))
    res["gain_100ms"]["vindu_1400_1500_mot_rene_referansevinduer_delomraader"] = {k: st(v) for k, v in sv.items()}

    # 2b. per reflector, on block-mean traces
    bok = np.zeros((NBI, NBX), bool)
    for r in range(NBI):
        a, b = r * BI, min(NIL, (r + 1) * BI)
        bok[r] = live[a:b].reshape(b - a, NBX, BX).mean((0, 2)) > 0.9
    bok[NBI - 1] &= (NIL - (NBI - 1) * BI) >= BI // 2
    E = {k: np.abs(hilbert(B[k].astype(np.float64), axis=-1)) for k in KS}
    res["reflektor"] = {}
    prof = {}
    for src in ("nea", "mid", "far", "ful"):
        rec = reflectors(B, E, src, bok)
        rr = rec[:, 0].astype(int)
        cc = rec[:, 1].astype(int)
        d = {}
        for an, bm in BAREA.items():
            d[an] = zr(rec[bm[rr, cc]])
        subs = [zr(rec[BSUB[rr, cc] == s]) for s in range(NSUB)]
        d["delomraader"] = {k: st([x["zr"][k] for x in subs if x]) for k in RATK + ("near_pk", "mid_pk", "far_pk")}
        res["reflektor"][src] = d
        # profile against time in 20 ms bins (for the figure)
        tb = np.arange(1100, 1900, 20)
        pr = {q: [] for q in ("mn_rms_med", "fn_rms_med", "mn_ls", "fn_ls", "mn_iv", "fn_iv", "mn_env_med", "fn_env_med", "pk_n", "pk_m", "pk_f", "c_nm", "c_nf", "n")}
        for t_ in tb:
            e = est(rec[(rec[:, 2] >= t_) & (rec[:, 2] < t_ + 20)])
            for q in pr:
                pr[q].append(float("nan") if e is None else e[q])
        prof[src] = {q: np.array(v, float) for q, v in pr.items()}
        prof[src]["t"] = tb + 10.0

    # block-mean traces: correlation between the stacks in the zone and the reference windows (zero lag)
    cz = {}
    for nm_, (k1, k2) in (("nm", ("nea", "mid")), ("nf", ("nea", "far")), ("mf", ("mid", "far"))):
        for wn_, wm_ in (("zone", ZM), ("ref", RM)):
            a_, b_ = B[k1][bok][:, wm_].astype(np.float64), B[k2][bok][:, wm_].astype(np.float64)
            cz[f"{nm_}_{wn_}"] = sdiv((a_ * b_).sum(), np.sqrt((a_ * a_).sum() * (b_ * b_).sum()))
    res["korrelasjon_blokktraser_nullforsinkelse"] = cz
    # instrument estimator without any alignment (zero lag), whole windows, block-mean traces:
    # mid/near = sum(mid*far)/sum(near*far), far/near = sum(far*mid)/sum(near*mid); zone divided by reference
    iv = {}
    for wn_, wm_ in (("zone", ZM), ("ref", RM)):
        n_, m_, f_ = (B[k][bok][:, wm_].astype(np.float64) for k in KS)
        iv[wn_] = dict(mn_iv=sdiv((m_ * f_).sum(), (n_ * f_).sum()), fn_iv=sdiv((f_ * m_).sum(), (n_ * m_).sum()),
                       mn_ls=sdiv((m_ * n_).sum(), (n_ * n_).sum()), fn_ls=sdiv((f_ * n_).sum(), (n_ * n_).sum()),
                       mn_rms=ssqrt(sdiv((m_ * m_).sum(), (n_ * n_).sum())), fn_rms=ssqrt(sdiv((f_ * f_).sum(), (n_ * n_).sum())),
                       near_rms=ssqrt((n_ * n_).mean()), mid_rms=ssqrt((m_ * m_).mean()), far_rms=ssqrt((f_ * f_).mean()))
    iv["zr"] = {q: sdiv(iv["zone"][q], iv["ref"][q]) for q in iv["zone"]}
    res["blokktraser_vindu_nullforsinkelse"] = iv

    # 4. the copy of the package 1256-1356 ms in the window 1366-1466 ms, stack by stack (block-mean traces)
    y0, y1 = (1366 - 800) // 2, (1466 - 800) // 2
    cp = {}
    for k in KS:
        Bk = B[k].astype(np.float64)
        y = Bk[..., y0:y1 + 1]
        Ls = np.arange(46, 68)
        syp = np.stack([(y * Bk[..., y0 - L:y1 + 1 - L]).sum(-1) for L in Ls], -1)
        spp = np.stack([(Bk[..., y0 - L:y1 + 1 - L] ** 2).sum(-1) for L in Ls], -1)
        syy = (y * y).sum(-1)
        rho = syp[bok].sum(0) / np.sqrt(spp[bok].sum(0) * syy[bok].sum())
        j = int(rho.argmin())
        subc = [sdiv(syp[bok & (BSUB == s)][:, j].sum(), spp[bok & (BSUB == s)][:, j].sum()) for s in range(NSUB)]
        subr = [sdiv(syp[bok & (BSUB == s)][:, j].sum(), np.sqrt(spp[bok & (BSUB == s)][:, j].sum() * syy[bok & (BSUB == s)].sum())) for s in range(NSUB)]
        cp[k] = dict(lag_ms=float(Ls[j] * 2), coef=sdiv(syp[bok][:, j].sum(), spp[bok][:, j].sum()), rho=float(rho[j]),
                     rho2=float(rho[j] ** 2), coef_delomraader=st(subc), rho_delomraader=st(subr),
                     rho_by_lag={str(int(L * 2)): float(v) for L, v in zip(Ls, rho)})
    res["kopi_1366_1466_av_1256_1356"] = cp
    res["seconds"] = time.time() - t0
    np.savez_compressed(os.path.join(OUT, f"profiles_{ver}.npz"),
                        **{f"e2t_{k}_{a}": E2t[k][a] for k in KS for a in AREA}, tms=tms,
                        **{f"prof_{s}_{q}": v for s, d in prof.items() for q, v in d.items()},
                        **{f"gmean_{k}": B[k][bok].astype(np.float64).mean(0) for k in KS + ("ful",)})
    if keepB is not None:
        keepB[ver] = {k: B[k][bok].copy() for k in ("nea", "far")} if bok.all() or True else None
        keepB[ver]["bok"] = bok
    return res


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    todo = sys.argv[1:] or list(VERS)
    rp = os.path.join(OUT, "results.json")
    RES = json.load(open(rp)) if os.path.exists(rp) else {}
    RES["_meta"] = dict(zone_ms=[1408, 1492], ref_ms="1100-1900 utenom 1380-1560", dt_ms=2.0, cube_shape=[NIL, NXL, NS],
                        trace_spacing_m=12.5, block_traces=[BI, BX], subarea_traces=[83, 78], n_subareas=NSUB,
                        well_index=[WELL_IL, WELL_XL], area_bronn="IL-indeks 0-137, XL-indeks 63-223 (1,7 x 2,0 km)",
                        area_kube="hele kuben, 249 x 468 traser (3,1 x 5,85 km)")
    KB = {}
    for ver in todo:
        RES[ver] = rnd(run(ver, KB))
        print(ver, "done", RES[ver]["seconds"], "s", flush=True)
        json.dump(RES, open(rp, "w"), indent=1)
    if "p07" in KB and "p10" in KB:
        ok = KB["p07"]["bok"] & KB["p10"]["bok"]
        out = {}
        for k in ("nea", "far"):
            a = np.zeros((NBI, NBX, NT), np.float32)
            b = np.zeros((NBI, NBX, NT), np.float32)
            a[KB["p07"]["bok"]] = KB["p07"][k]
            b[KB["p10"]["bok"]] = KB["p10"][k]
            a, b = a[ok][:, 150:551].astype(np.float64), b[ok][:, 150:551].astype(np.float64)
            cc = {}
            for L in range(-3, 4):
                u, v = (a[:, :a.shape[1] - L], b[:, L:]) if L >= 0 else (a[:, -L:], b[:, :b.shape[1] + L])
                cc[str(2 * L)] = sdiv((u * v).sum(), np.sqrt((u * u).sum() * (v * v).sum()))
            out[k] = cc
        RES["p07_mot_p10_korrelasjon_blokktraser_1100_1900_etter_forsinkelse_ms"] = rnd(out)
        json.dump(RES, open(rp, "w"), indent=1)
    print("ALL DONE", flush=True)
