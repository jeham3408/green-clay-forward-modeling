#!/usr/bin/env python3
"""post.py (kontroll): every reported number, computed from acc_<ver>.npz written by pass.py. Writes results.json."""
import json, os
import numpy as np
from scipy.signal import hilbert
from scipy.stats import spearmanr

OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
VERS = ("p07", "p10", "p07_dmE", "p10_dmE", "p07_dm", "p07_dmA", "p07_dmB", "p07_dmD")
KS = ("nea", "mid", "far")
NT = 600
t = 800.0 + 2.0 * np.arange(NT)
XLAGS = np.arange(-4, 5)
ALAGS = np.arange(48, 67)
PAIRS = (("nea", "mid"), ("nea", "far"), ("mid", "far"))


def win(a, b):
    return (t >= a) & (t <= b)


ZM = win(1408, 1492)
RM = ((t >= 1100) & (t < 1380)) | ((t > 1560) & (t <= 1900))
RCLEAN = ((t >= 1100) & (t < 1250)) | ((t >= 1600) & (t <= 1900))     # reference without the strong package and its tail
WINS = {"sone_1408_1492": ZM, "kopivindu_1366_1466": win(1366, 1466), "over_leira_1366_1404": win(1366, 1404),
        "leirlaget_1412_1446": win(1412, 1446), "under_leira_1450_1492": win(1450, 1492),
        "under_sonen_1494_1560": win(1494, 1560), "pakken_1256_1356": win(1256, 1356)}


def rnd(o):
    if isinstance(o, dict):
        return {str(k): rnd(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [rnd(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 4)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return rnd(o.tolist())
    return o


def stat(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return dict(n=len(v), median=np.median(v), mean=v.mean(), sd=v.std(ddof=1), se=v.std(ddof=1) / np.sqrt(len(v)),
                min=v.min(), max=v.max())


def level(P, wm, rm):
    """ratio of mean of P over window wm to mean over rm, P is (..., NT); returns array over leading axes."""
    return P[..., wm].mean(-1) / P[..., rm].mean(-1)


def sq(x):
    return np.sqrt(np.where(x > 0, x, np.nan))


def analyse(ver):
    d = np.load(os.path.join(OUT, f"acc_{ver}.npz"))
    R = {"live_fraction": float(d["live"].mean())}
    ss = {k: d[f"ss_{k}"] for k in KS}
    ab = {k: d[f"ab_{k}"] for k in KS}
    nx = {k: d[f"nx_{k}"] for k in KS}
    ni = {k: d[f"ni_{k}"] for k in KS}
    cx = {p: d[f"cx_{p[0]}_{p[1]}"] for p in PAIRS}
    ac = {k: d[f"ac_{k}"] for k in KS}
    B = {k: d[f"B_{k}"].astype(np.float64) for k in KS}

    # ---- 1. level of each stack, zone / reference, for several windows and estimators
    lev = {}
    for wn, wm in WINS.items():
        lev[wn] = {}
        for k in KS:
            tot = dict(rms=sq(level(ss[k].sum(0), wm, RM)), abs=level(ab[k].sum(0), wm, RM),
                       nabo_xl=sq(level(nx[k].sum(0), wm, RM)), nabo_il=sq(level(ni[k].sum(0), wm, RM)),
                       rms_ren_ref=sq(level(ss[k].sum(0), wm, RCLEAN)))
            Bk = B[k].reshape(-1, NT)
            tot["blokktrase_rms"] = sq(level((Bk ** 2).sum(0), wm, RM))
            env = np.abs(hilbert(Bk, axis=-1))
            tot["blokktrase_env_median"] = np.median(env[:, wm]) / np.median(env[:, RM])
            sub = sq(level(ss[k], wm, RM))
            subn = sq(level(nx[k], wm, RM))
            lev[wn][k] = dict(hele=tot, delomraader_rms=stat(sub), delomraader_nabo=stat(subn))
        for a, b in (("mid", "nea"), ("far", "nea"), ("far", "mid")):
            r = {e: lev[wn][a]["hele"][e] / lev[wn][b]["hele"][e] for e in lev[wn][a]["hele"]}
            sub = sq(level(ss[a], wm, RM)) / sq(level(ss[b], wm, RM))
            subn = sq(level(nx[a], wm, RM)) / sq(level(nx[b], wm, RM))
            lev[wn][f"{a}/{b}"] = dict(hele=r, delomraader_rms=stat(sub), delomraader_nabo=stat(subn))
    R["nivaa_og_forhold"] = lev

    # ---- 2. least squares / instrument estimators, trace level, lag of best overall correlation per pair
    ls = {}
    lagbest = {}
    for (a, b) in PAIRS:
        c = cx[(a, b)].sum(1)                       # (lags, NT): sum a(t) b(t+L)
        saa, sbb = ss[a].sum(0), ss[b].sum(0)
        o = {}
        for wn, wm in (("sone", ZM), ("ref", RM)):
            cc = c[:, wm].sum(1) / np.sqrt(saa[wm].sum() * sbb[wm].sum())
            j = int(np.argmax(cc))
            if 0 < j < len(cc) - 1:
                y0, y1, y2 = cc[j - 1], cc[j], cc[j + 1]
                frac = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            else:
                frac = 0.0
            o[wn] = dict(korr_null=cc[4], korr_best=cc[j], lag_ms=2 * (XLAGS[j] + frac), lag_sample=int(XLAGS[j]),
                         paa_kanten=bool(j in (0, len(cc) - 1)))
        lagbest[(a, b)] = o
        ls[f"{a}-{b}"] = o
    R["innretting_traser"] = ls

    def lsq(a, b, wm, j):
        """b on a: slope = sum(ab)/sum(aa); inverse = sum(bb)/sum(ab); rms = sqrt(bb/aa)."""
        c = cx[(a, b)].sum(1)[j][wm].sum()
        saa, sbb = ss[a].sum(0)[wm].sum(), ss[b].sum(0)[wm].sum()
        return dict(ls=c / saa, inv=sbb / c, rms=np.sqrt(sbb / saa))
    est = {}
    for (a, b) in (("nea", "mid"), ("nea", "far")):
        for mode in ("null", "best"):
            jz = 4 if mode == "null" else int(np.where(XLAGS == lagbest[(a, b)]["sone"]["lag_sample"])[0][0])
            jr = 4 if mode == "null" else int(np.where(XLAGS == lagbest[(a, b)]["ref"]["lag_sample"])[0][0])
            z, r = lsq(a, b, ZM, jz), lsq(a, b, RM, jr)
            est[f"{b}/{a}_{mode}"] = dict(sone=z, ref=r, sone_over_ref={q: z[q] / r[q] for q in z})
    # instrument (third stack) estimators at zero lag, trace level
    nm, nf, mf = (cx[p].sum(1)[4] for p in PAIRS)
    iv = {}
    for wn, wm in (("sone", ZM), ("ref", RM)):
        iv[wn] = dict(mid_near=mf[wm].sum() / nf[wm].sum(), far_near=mf[wm].sum() / nm[wm].sum())
    iv["sone_over_ref"] = {q: iv["sone"][q] / iv["ref"][q] for q in iv["sone"]}
    est["instrument_null_lag"] = iv
    R["minste_kvadrater_traser"] = est

    # block-mean traces: correlations at zero lag, instrument estimator
    Bn, Bm, Bf = (B[k].reshape(-1, NT) for k in KS)
    bc = {}
    for wn, wm in (("sone", ZM), ("ref", RM)):
        n_, m_, f_ = Bn[:, wm], Bm[:, wm], Bf[:, wm]
        c = lambda x, y: (x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum())
        bc[wn] = dict(korr_nm=c(n_, m_), korr_nf=c(n_, f_), korr_mf=c(m_, f_),
                      iv_mid_near=(m_ * f_).sum() / (n_ * f_).sum(), iv_far_near=(m_ * f_).sum() / (n_ * m_).sum(),
                      rms_mid_near=np.sqrt((m_ * m_).sum() / (n_ * n_).sum()), rms_far_near=np.sqrt((f_ * f_).sum() / (n_ * n_).sum()))
    bc["sone_over_ref"] = {q: bc["sone"][q] / bc["ref"][q] for q in bc["sone"] if not q.startswith("korr")}
    R["blokktraser"] = bc

    # per reflector, neutral picks: envelope peaks of the sum of the three block-mean stacks (each scaled by its own
    # reference rms), envelope maximum of each stack within +-8 ms of the pick
    E = {k: np.abs(hilbert(B[k].reshape(-1, NT), axis=-1)) for k in KS}
    S = sum(B[k].reshape(-1, NT) / np.sqrt((B[k].reshape(-1, NT)[:, RM] ** 2).mean()) for k in KS)
    ES = np.abs(hilbert(S, axis=-1))
    i = np.arange(150, 551)
    thr = 0.5 * np.median(ES[:, 150:551], axis=1, keepdims=True)
    pk = (ES[:, i] > ES[:, i - 1]) & (ES[:, i] >= ES[:, i + 1]) & (ES[:, i] > thr)
    rr, ii = np.nonzero(pk)
    ii = i[ii]
    pv = {k: np.max(np.stack([E[k][rr, ii + o] for o in range(-4, 5)]), 0) for k in KS}
    tp = t[ii]
    inz = (tp >= 1408) & (tp <= 1492)
    inr = ((tp >= 1100) & (tp < 1380)) | ((tp > 1560) & (tp <= 1900))
    rf = dict(n_sone=int(inz.sum()), n_ref=int(inr.sum()))
    for k in KS:
        rf[k] = np.median(pv[k][inz]) / np.median(pv[k][inr])
    for a, b in (("mid", "nea"), ("far", "nea"), ("far", "mid")):
        q = pv[a] / pv[b]
        rf[f"{a}/{b}_median_av_forhold"] = np.median(q[inz]) / np.median(q[inr])
        rf[f"{a}/{b}_forhold_av_sum"] = (pv[a][inz].sum() / pv[b][inz].sum()) / (pv[a][inr].sum() / pv[b][inr].sum())
    R["reflektor_envelopetopper_noytrale_plukk"] = rf

    # ---- 3. gain: rms in 100 ms windows 800-2000 ms
    g = {"t_start": list(range(800, 2000, 100))}
    for k in KS:
        g[k] = np.sqrt(ss[k].sum(0).reshape(12, 50).sum(1) / (d["cnt"].sum() * 50))
    g["mid/nea"] = g["mid"] / g["nea"]
    g["far/nea"] = g["far"] / g["nea"]
    g["far/mid"] = g["far"] / g["mid"]
    R["gain_100ms"] = g
    # fine profile in 20 ms bins (for the figure and for the location of the minimum)
    f20 = {"t_mid": (800 + 20 * np.arange(60) + 10).tolist()}
    for k in KS:
        f20[k] = np.sqrt(ss[k].sum(0).reshape(60, 10).sum(1) / (d["cnt"].sum() * 10))
        f20[k + "_nabo"] = sq(nx[k].sum(0).reshape(60, 10).sum(1) / (d["cnt_nx"].sum() * 10))
    R["profil_20ms"] = f20

    # ---- selection effect: slide an 84 ms window (43 samples) over 1100-1900 ms, ratios relative to the whole 1100-1900
    ALL = win(1100, 1900)
    cs = {k: np.concatenate([[0], np.cumsum(ss[k].sum(0))]) for k in KS}
    i0 = np.arange(150, 551 - 43 + 1)
    tw0 = t[i0]
    pw = {k: (cs[k][i0 + 43] - cs[k][i0]) / 43 for k in KS}
    al = {k: ss[k].sum(0)[ALL].mean() for k in KS}
    sl = {}
    for a, b in (("mid", "nea"), ("far", "nea"), ("far", "mid")):
        q = np.sqrt((pw[a] / al[a]) / (pw[b] / al[b]))
        j = int(np.argmin(q))
        jz = int(np.where(tw0 == 1408)[0][0])
        # non-overlapping windows that do not touch 1366-1560 (copy window and its tail): starts every 86 ms
        starts = [s for s in range(1100, 1900 - 84 + 1, 86) if (s + 84 < 1366 or s > 1560)]
        vals = np.array([q[int(np.where(tw0 == s)[0][0])] for s in starts])
        sl[f"{a}/{b}"] = dict(min_verdi=q[j], min_vindu_start_ms=tw0[j], sone_verdi=q[jz],
                              andel_vinduer_lavere_enn_sonen=float((q < q[jz]).mean()),
                              uavhengige_vinduer_start=starts, uavhengige_vinduer_verdier=vals,
                              sd_ln_uavhengige=np.log(vals).std(ddof=1),
                              z_sone=(np.log(q[jz]) - np.log(vals).mean()) / np.log(vals).std(ddof=1))
    for k in KS:
        q = np.sqrt(pw[k] / al[k])
        sl[k] = dict(max_verdi=q.max(), max_vindu_start_ms=tw0[int(np.argmax(q))], sone_verdi=q[int(np.where(tw0 == 1408)[0][0])])
    R["glidende_vindu_84ms"] = sl

    # ---- 4. the copy: regression of a(t+L) on a(t), source window 1256-1356, per stack, trace level
    cp = {}
    src = win(1256, 1356)
    for k in KS:
        a_ = ac[k].sum(1)                            # (lags, NT) sum a(t) a(t+L)
        p = ss[k].sum(0)
        coef, rho2 = [], []
        for j, L in enumerate(ALAGS):
            idx = np.nonzero(src)[0]
            num = a_[j, idx].sum()
            coef.append(num / p[idx].sum())
            rho2.append(num ** 2 / (p[idx].sum() * p[idx + L].sum()) * np.sign(num))
        coef, rho2 = np.array(coef), np.array(rho2)
        j = int(np.argmin(coef))
        cp[k] = dict(lag_ms=2 * int(ALAGS[j]), koef=coef[j], andel_energi=abs(rho2[j]), koef_alle_lag=coef,
                     lag_ms_alle=(2 * ALAGS).tolist())
        # null: the same measurement for every 100 ms source window between 1100 and 1900-100-lag, at this stack's lag
        Lb = int(ALAGS[j])
        starts = np.arange(150, 551 - 51 - Lb, 5)
        cn = []
        for s0 in starts:
            idx = np.arange(s0, s0 + 51)
            num = a_[j, idx].sum()
            cn.append([t[s0], num / p[idx].sum(), np.sign(num) * num ** 2 / (p[idx].sum() * p[idx + Lb].sum())])
        cn = np.array(cn)
        cp[k]["null_start_ms"] = cn[:, 0]
        cp[k]["null_koef"] = cn[:, 1]
        cp[k]["null_signert_andel"] = cn[:, 2]
        far_from = (np.abs(cn[:, 0] - 1256) > 100)
        cp[k]["null_andel_utenfor_pakken"] = stat(cn[far_from, 2])
        # copy-explained energy in the zone itself (target 1408-1492, source 110 ms above), and in 1366-1466
        for wn, wm in (("sone_1408_1492", ZM), ("kopivindu_1366_1466", win(1366, 1466))):
            idx = np.nonzero(wm)[0] - Lb
            num = a_[j, idx].sum()
            cp[k]["andel_" + wn] = num ** 2 / (p[idx].sum() * p[idx + Lb].sum())
            cp[k]["koef_" + wn] = num / p[idx].sum()
        # block-mean traces, same thing
        Bk = B[k].reshape(-1, NT)
        idx = np.nonzero(src)[0]
        num = (Bk[:, idx] * Bk[:, idx + Lb]).sum()
        cp[k]["blokktrase_koef"] = num / (Bk[:, idx] ** 2).sum()
        cp[k]["blokktrase_andel"] = num ** 2 / ((Bk[:, idx] ** 2).sum() * (Bk[:, idx + Lb] ** 2).sum())
    R["kopi"] = cp
    # levels after removing the copy-explained share of the energy in the zone (first-order copy only, one coefficient)
    kr = {}
    for k in KS:
        kr[k] = lev["sone_1408_1492"][k]["hele"]["rms"] * np.sqrt(1 - cp[k]["andel_sone_1408_1492"])
    kr["mid/nea"] = kr["mid"] / kr["nea"]
    kr["far/nea"] = kr["far"] / kr["nea"]
    R["sone_uten_kopiandel"] = kr

    # ---- 5. does mid/near follow the predicted copy load over the whole 1100-1900 ms? 20 ms bins
    pn = ss["nea"].sum(0)
    L = 55
    bins = [(a, a + 20) for a in range(1100, 1900, 20)]
    xs, ym, yf, tb = [], [], [], []
    for a, b in bins:
        m = (t >= a) & (t < b)
        idx = np.nonzero(m)[0]
        xs.append(pn[idx - L].sum() / pn[idx].sum())             # power 110 ms above relative to local power, near
        ym.append(0.5 * np.log(ss["mid"].sum(0)[idx].sum() / pn[idx].sum()))
        yf.append(0.5 * np.log(ss["far"].sum(0)[idx].sum() / pn[idx].sum()))
        tb.append(a + 10)
    xs, ym, yf, tb = map(np.array, (xs, ym, yf, tb))
    outz = (tb < 1360) | (tb > 1500)
    R["kopilast_mot_forhold_20ms"] = dict(
        t=tb, kopilast=xs, ln_mid_near=ym, ln_far_near=yf,
        spearman_mid_alle=spearmanr(xs, ym)[0], spearman_far_alle=spearmanr(xs, yf)[0],
        spearman_mid_utenfor_1360_1500=spearmanr(xs[outz], ym[outz])[0],
        spearman_far_utenfor_1360_1500=spearmanr(xs[outz], yf[outz])[0], n_alle=len(tb), n_utenfor=int(outz.sum()))
    return R


res = {}
for v in VERS:
    if os.path.exists(os.path.join(OUT, f"acc_{v}.npz")):
        res[v] = analyse(v)
# sub-area consistency between the processings and filter effect on absolute rms (zone and reference separately)
fe = {}
raw = np.load(os.path.join(OUT, "acc_p07.npz"))
for v in VERS:
    if v in ("p07", "p10") or v not in res:
        continue
    base = np.load(os.path.join(OUT, "acc_p10.npz")) if v.startswith("p10") else raw
    dm = np.load(os.path.join(OUT, f"acc_{v}.npz"))
    fe[v] = {k: dict(sone=float(np.sqrt(dm[f"ss_{k}"].sum(0)[ZM].sum() / base[f"ss_{k}"].sum(0)[ZM].sum())),
                     ref=float(np.sqrt(dm[f"ss_{k}"].sum(0)[RM].sum() / base[f"ss_{k}"].sum(0)[RM].sum()))) for k in KS}
res["filter_beholdt_amplitude"] = fe
json.dump(rnd(res), open(os.path.join(OUT, "results.json"), "w"), indent=1, ensure_ascii=False)

# ---- print the main table
def g(v, wn, key, e="rms"):
    return res[v]["nivaa_og_forhold"][wn][key]["hele"][e]
for v in res:
    if v == "filter_beholdt_amplitude":
        continue
    print("=====", v)
    for wn in WINS:
        print(f" {wn:24s}", " ".join(f"{key}:" + "/".join(f"{g(v, wn, key, e):.2f}" for e in ("rms", "abs", "nabo_xl", "nabo_il", "blokktrase_rms", "blokktrase_env_median", "rms_ren_ref")) for key in ("nea", "mid", "far", "mid/nea", "far/nea", "far/mid")))
    z = res[v]["nivaa_og_forhold"]["sone_1408_1492"]
    for key in ("nea", "mid", "far", "mid/nea", "far/nea"):
        s, sn = z[key]["delomraader_rms"], z[key]["delomraader_nabo"]
        print(f"  delomr {key}: rms median {s['median']:.3f} se {s['se']:.3f} min {s['min']:.2f} max {s['max']:.2f} | nabo median {sn['median']:.3f} se {sn['se']:.3f}")
    print("  innretting:", json.dumps(rnd(res[v]["innretting_traser"])))
    print("  LSQ:", json.dumps(rnd({k: x["sone_over_ref"] for k, x in res[v]["minste_kvadrater_traser"].items()})))
    print("  IV raw:", json.dumps(rnd(res[v]["minste_kvadrater_traser"]["instrument_null_lag"])))
    print("  blokk:", json.dumps(rnd(res[v]["blokktraser"])))
    print("  refl:", json.dumps(rnd(res[v]["reflektor_envelopetopper_noytrale_plukk"])))
    gg = res[v]["gain_100ms"]
    for key in ("nea", "mid", "far", "mid/nea", "far/nea", "far/mid"):
        print(f"  gain {key:8s}", " ".join(f"{x:.3f}" for x in gg[key]))
    s = res[v]["glidende_vindu_84ms"]
    print("  glidende:", json.dumps(rnd({k: {q: x[q] for q in x if not q.startswith("uavh")} for k, x in s.items()})))
    print("  glidende uavh mid/nea:", rnd(s["mid/nea"]["uavhengige_vinduer_verdier"]), "far/nea:", rnd(s["far/nea"]["uavhengige_vinduer_verdier"]))
    for k in KS:
        c = res[v]["kopi"][k]
        print(f"  kopi {k}: lag {c['lag_ms']} koef {c['koef']:.3f} andel {c['andel_energi']:.3f} | sone-andel {c['andel_sone_1408_1492']:.3f} koef {c['koef_sone_1408_1492']:.3f} | 1366-1466 andel {c['andel_kopivindu_1366_1466']:.3f} | blokk koef {c['blokktrase_koef']:.3f} andel {c['blokktrase_andel']:.3f} | null utenfor pakken: median {c['null_andel_utenfor_pakken']['median']:.3f} min {c['null_andel_utenfor_pakken']['min']:.3f} max {c['null_andel_utenfor_pakken']['max']:.3f}")
        print("     koef alle lag:", " ".join(f"{x:.2f}" for x in c["koef_alle_lag"]))
    print("  sone uten kopiandel:", json.dumps(rnd(res[v]["sone_uten_kopiandel"])))
    kl = res[v]["kopilast_mot_forhold_20ms"]
    print("  kopilast spearman:", {q: round(float(kl[q]), 3) for q in kl if q.startswith("spearman")}, kl["n_alle"], kl["n_utenfor"])
print("filter:", json.dumps(rnd(fe)))
