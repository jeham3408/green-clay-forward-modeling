#!/usr/bin/env python3
"""post2.py (kontroll): extra tests on top of post.py. Adds the key "ekstra" to results.json and draws two figures.
 - where in time the deficits sit, relative to the soft layers of the model (converted to two-way time with the
   time-depth curve of 15/9-13) and relative to where the copy of the strong package lands
 - copy share with the clay window as target; sub-area test; per-trace medians; filter bookkeeping
"""
import json, os
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
SL = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
VERS = ("p07", "p10", "p07_dmE", "p10_dmE", "p07_dm", "p07_dmA", "p07_dmB", "p07_dmD")
KS = ("nea", "mid", "far")
NT = 600
t = 800.0 + 2.0 * np.arange(NT)
ALAGS = np.arange(48, 67)
win = lambda a, b: (t >= a) & (t <= b)
ZM = win(1408, 1492)
RM = ((t >= 1100) & (t < 1380)) | ((t > 1560) & (t <= 1900))


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


sq = lambda x: np.sqrt(np.where(x > 0, x, np.nan))
lev = lambda P, wm, rm=RM: P[..., wm].mean(-1) / P[..., rm].mean(-1)

# soft layers of the model (metres below sea level, from prediksjoner.json) -> two-way time with the well's TD curve
w = json.load(open(SL + "/cubes/sleipner94/wells.json"))["wells"][0]
md = np.array(w["td"]["md"], float)
tw = np.array(w["td"]["twt"], float)
ok = np.isfinite(md) & np.isfinite(tw)
soft_m = [(1353, 1362), (1366, 1373), (1375, 1399), (1425, 1435)]
soft_ms = [[float(np.interp(a + w["kb"], md[ok], tw[ok])), float(np.interp(b + w["kb"], md[ok], tw[ok]))] for a, b in soft_m]
s_top, s_bot = soft_ms[0][0], soft_ms[-1][1]
HALF = 16.0                                   # half length of the wavelet, ms (one side lobe of a 30-35 Hz wavelet)
a_hi = 2 * np.floor((s_top - HALF) / 2)
b_lo = 2 * np.ceil((s_bot + HALF) / 2)
WIN2 = {f"over_myke_lag_1366_{a_hi:.0f}": win(1366, a_hi), f"myke_lag_{a_hi + 2:.0f}_{b_lo - 2:.0f}": win(a_hi + 2, b_lo - 2),
        f"under_myke_lag_{b_lo:.0f}_1560": win(b_lo, 1560), "lenger_ned_1562_1700": win(1562, 1700)}
EX = {"myke_lag_ms_fra_td_kurven": soft_ms, "vinduer": {k: [float(t[m][0]), float(t[m][-1])] for k, m in WIN2.items()}}

A = {v: np.load(os.path.join(OUT, f"acc_{v}.npz")) for v in VERS}
for v in VERS:
    d = A[v]
    ss = {k: d[f"ss_{k}"] for k in KS}
    nx = {k: d[f"nx_{k}"] for k in KS}
    E = {}
    # 1. windows relative to the soft layers
    W = {}
    for wn, wm in WIN2.items():
        o = {}
        for k in KS:
            o[k] = dict(rms=sq(lev(ss[k].sum(0), wm)), nabo=sq(lev(nx[k].sum(0), wm)))
        for a, b in (("mid", "nea"), ("far", "nea"), ("far", "mid")):
            o[f"{a}/{b}"] = {e: o[a][e] / o[b][e] for e in ("rms", "nabo")}
            sub = sq(lev(ss[a], wm)) / sq(lev(ss[b], wm))
            o[f"{a}/{b}"]["delomraader_rms_median"] = np.nanmedian(sub)
            o[f"{a}/{b}"]["delomraader_rms_se"] = np.nanstd(sub, ddof=1) / np.sqrt(len(sub))
        W[wn] = o
    E["vinduer_mot_myke_lag"] = W
    names = list(WIN2)
    zs = {}
    for q in ("mid/nea", "far/nea", "far/mid"):
        for e in ("rms", "nabo"):
            zs[f"{q}_{e}"] = W[names[1]][q][e] / np.sqrt(W[names[0]][q][e] * W[names[2]][q][e])
    for k in KS:
        for e in ("rms", "nabo"):
            zs[f"{k}_{e}"] = W[names[1]][k][e] / np.sqrt(W[names[0]][k][e] * W[names[2]][k][e])
    E["myke_lag_mot_naboene_over_og_under"] = zs
    # the same per sub-area, for an uncertainty that respects spatial correlation (24 sub-areas)
    sz = {}
    for q, (a, b) in (("mid/nea", ("mid", "nea")), ("far/nea", ("far", "nea")), ("far/mid", ("far", "mid"))):
        r = [sq(lev(ss[a], WIN2[n])) / sq(lev(ss[b], WIN2[n])) for n in names[:3]]
        x = r[1] / np.sqrt(r[0] * r[2])
        sz[q] = dict(median=np.median(x), se=x.std(ddof=1) / np.sqrt(len(x)), min=x.min(), max=x.max())
    E["myke_lag_mot_naboene_delomraader_rms"] = sz

    # 2. copy share of the near, mid and far energy with the clay window and the zone as targets (lags 96-132 ms)
    cp = {}
    for k in KS:
        a_ = d[f"ac_{k}"].sum(1)
        p = ss[k].sum(0)
        o = {}
        for wn, wm in (("leirlaget_1412_1446", win(1412, 1446)), ("sone_1408_1492", ZM)):
            best = None
            for j, L in enumerate(ALAGS):
                idx = np.nonzero(wm)[0] - L
                num = a_[j, idx].sum()
                c = num / p[idx].sum()
                r2 = num ** 2 / (p[idx].sum() * p[idx + L].sum())
                if best is None or c < best[1]:
                    best = (2 * int(L), c, r2)
            o[wn] = dict(lag_ms=best[0], koef=best[1], andel=best[2])
        # reference: the same at this lag for all 84 ms targets fully inside the reference windows
        L = o["sone_1408_1492"]["lag_ms"] // 2
        j = int(np.where(ALAGS == L)[0][0])
        vals = []
        for s0 in range(1100 + 2 * L, 1900 - 84 + 1, 28):
            if (s0 + 84 < 1380) or (s0 > 1560):
                idx = np.nonzero(win(s0, s0 + 84))[0] - L
                num = a_[j, idx].sum()
                vals.append([s0, num / p[idx].sum(), np.sign(num) * num ** 2 / (p[idx].sum() * p[idx + L].sum())])
        vals = np.array(vals)
        o["referanse_maalvinduer_start_ms"] = vals[:, 0]
        o["referanse_koef"] = vals[:, 1]
        o["referanse_signert_andel"] = vals[:, 2]
        cp[k] = o
    E["kopiandel_maalvindu"] = cp
    # near level in the zone if the copy share is taken out of both the zone and (median share) the reference
    k = "nea"
    zl = sq(lev(ss[k].sum(0), ZM))
    refshare = float(np.median(np.clip(-cp[k]["referanse_signert_andel"], 0, None)))
    E["near_sone_uten_kopi_i_sone_og_referanse"] = dict(
        near_raa=zl, andel_sone=cp[k]["sone_1408_1492"]["andel"], andel_referanse_median=refshare,
        near_korrigert=zl * np.sqrt(1 - cp[k]["sone_1408_1492"]["andel"]) / np.sqrt(1 - refshare))

    # 3. sub-areas: does the near excess follow the copy? (24 sub-areas)
    L = cp["nea"]["sone_1408_1492"]["lag_ms"] // 2
    j = int(np.where(ALAGS == L)[0][0])
    idx = np.nonzero(ZM)[0] - L
    acn = d["ac_nea"][j]
    share = np.sign(acn[:, idx].sum(1)) * acn[:, idx].sum(1) ** 2 / (ss["nea"][:, idx].sum(1) * ss["nea"][:, idx + L].sum(1))
    nlev = sq(lev(ss["nea"], ZM))
    mn = sq(lev(ss["mid"], ZM)) / nlev
    fn = sq(lev(ss["far"], ZM)) / nlev
    pk = sq(lev(ss["nea"], win(1256, 1356)))
    E["delomraader_kopi"] = dict(
        n=len(share), kopiandel_min=float((-share).min()), kopiandel_max=float((-share).max()),
        spearman_kopiandel_mot_near=spearmanr(-share, nlev)[0], spearman_kopiandel_mot_mid_near=spearmanr(-share, mn)[0],
        spearman_kopiandel_mot_far_near=spearmanr(-share, fn)[0], spearman_pakkestyrke_mot_near_sone=spearmanr(pk, nlev)[0],
        spearman_pakkestyrke_mot_mid_near=spearmanr(pk, mn)[0])

    # 4. per trace: window 1400-1500 ms against 1100-1400 and 1600-1900 ms (100 ms windows; approximates zone/reference)
    live = d["live"]
    pt = {}
    for k in KS:
        tr = d[f"tr_{k}"].astype(np.float64)[live]
        r = np.sqrt(tr[:, 6] / tr[:, [3, 4, 5, 8, 9, 10]].mean(1))
        pt[k] = r
    q = lambda x: dict(median=np.median(x), p16=np.percentile(x, 16), p84=np.percentile(x, 84))
    E["per_trase_1400_1500"] = {**{k: q(pt[k]) for k in KS}, "mid/nea": q(pt["mid"] / pt["nea"]), "far/nea": q(pt["far"] / pt["nea"]),
                                "n": int(live.sum())}
    # around the well: 161 x 161 traces centred on (57, 143), clipped to the cube
    il0, il1, xl0, xl1 = max(0, 57 - 80), 57 + 81, 143 - 80, 143 + 81
    wl = {}
    for k in KS:
        tr = d[f"tr_{k}"].astype(np.float64)[il0:il1, xl0:xl1].reshape(-1, 12).sum(0)
        wl[k] = np.sqrt(tr[6] / tr[[3, 4, 5, 8, 9, 10]].mean())
    wl["mid/nea"] = wl["mid"] / wl["nea"]
    wl["far/nea"] = wl["far"] / wl["nea"]
    wl["n_traser"] = int((il1 - il0) * (xl1 - xl0))
    E["rundt_bronnen_1400_1500"] = wl
    whole = {}
    for k in KS:
        tr = d[f"tr_{k}"].astype(np.float64).reshape(-1, 12).sum(0)
        whole[k] = np.sqrt(tr[6] / tr[[3, 4, 5, 8, 9, 10]].mean())
    whole["mid/nea"] = whole["mid"] / whole["nea"]
    whole["far/nea"] = whole["far"] / whole["nea"]
    E["hele_kuben_1400_1500"] = whole
    EX[v] = E

# 5. filter bookkeeping: how much of the change in near zone/reference comes from the reference going up
fb = {}
for v in VERS[2:]:
    base = A["p10"] if v.startswith("p10") else A["p07"]
    o = {}
    for k in KS:
        z = np.sqrt(A[v][f"ss_{k}"].sum(0)[ZM].sum() / base[f"ss_{k}"].sum(0)[ZM].sum())
        r = np.sqrt(A[v][f"ss_{k}"].sum(0)[RM].sum() / base[f"ss_{k}"].sum(0)[RM].sum())
        o[k] = dict(sone=z, ref=r, sone_over_ref=z / r, andel_av_endring_fra_referansen=np.log(r) / (np.log(r) - np.log(z)))
    rb = lambda a, b, d_: np.sqrt(lev(d_[f"ss_{a}"].sum(0), ZM) / lev(d_[f"ss_{b}"].sum(0), ZM))
    for q, (a, b) in (("mid/nea", ("mid", "nea")), ("far/nea", ("far", "nea"))):
        o[q + "_andel_av_ln_avvik_igjen"] = np.log(rb(a, b, A[v])) / np.log(rb(a, b, base))
    fb[v] = o
EX["filter"] = fb

res = json.load(open(os.path.join(OUT, "results.json")))
res["ekstra"] = rnd(EX)
json.dump(res, open(os.path.join(OUT, "results.json"), "w"), indent=1, ensure_ascii=False)

print("myke lag ms:", rnd(soft_ms), "vinduer:", EX["vinduer"])
for v in VERS:
    E = EX[v]
    print("=====", v)
    for wn, o in E["vinduer_mot_myke_lag"].items():
        print(f"  {wn:28s}", " ".join(f"{k}:{o[k]['rms']:.2f}/{o[k]['nabo']:.2f}" for k in ("nea", "mid", "far", "mid/nea", "far/nea", "far/mid")),
              f"| delomr far/mid {o['far/mid']['delomraader_rms_median']:.3f}±{o['far/mid']['delomraader_rms_se']:.3f} mid/nea {o['mid/nea']['delomraader_rms_median']:.3f}±{o['mid/nea']['delomraader_rms_se']:.3f}")
    print("  myke lag mot naboene:", rnd(E["myke_lag_mot_naboene_over_og_under"]))
    print("  delomr:", rnd(E["myke_lag_mot_naboene_delomraader_rms"]))
    for k in KS:
        c = E["kopiandel_maalvindu"][k]
        print(f"  kopi {k}: leir {rnd(c['leirlaget_1412_1446'])} sone {rnd(c['sone_1408_1492'])} ref andel {rnd(c['referanse_signert_andel'])}")
    print("  near uten kopi:", rnd(E["near_sone_uten_kopi_i_sone_og_referanse"]))
    print("  delomr kopi:", rnd(E["delomraader_kopi"]))
    print("  per trase:", rnd(E["per_trase_1400_1500"]))
    print("  bronn:", rnd(E["rundt_bronnen_1400_1500"]), "hele:", rnd(E["hele_kuben_1400_1500"]))
print("filter:", json.dumps(rnd(fb)))

# ---------------- figures
def prof(v, n=10):
    d = A[v]
    o = {}
    for k in KS:
        p = d[f"ss_{k}"].sum(0)
        o[k] = np.sqrt(p.reshape(NT // n, n).mean(1) / p[RM].mean())
    return 800 + 2 * n * np.arange(NT // n) + n, o


def shade(ax):
    ax.axvspan(1256, 1356, color="#f2d16b", alpha=0.35, lw=0)
    for a, b in soft_ms:
        ax.axvspan(a, b, color="0.55", alpha=0.45, lw=0)
    ax.axvspan(1366, 1466, ymin=0.0, ymax=0.04, color="#c0392b", lw=0)
    ax.set_xlim(1100, 1900)
    ax.grid(alpha=0.3)


C = {"nea": "#1f77b4", "mid": "#2ca02c", "far": "#d62728"}
N = {"nea": "near", "mid": "mid", "far": "far"}
fig, ax = plt.subplots(3, 1, figsize=(9, 9.2), sharex=True)
for i, (v, ttl) in enumerate((("p07", "a) Styrke mot tid, ufiltrert (94p07), delt på hver stakks referansenivå"),
                             ("p07_dmE", "b) Det samme etter multippelfilter E"))):
    tt, o = prof(v)
    for k in KS:
        ax[i].plot(tt, o[k], color=C[k], lw=1.8, label=N[k])
    shade(ax[i])
    ax[i].set_title(ttl, fontsize=10, loc="left")
    ax[i].set_ylabel("RMS / referanse")
    ax[i].legend(loc="upper right", ncol=3, fontsize=9)
for v, ls_ in (("p07", "-"), ("p10", ":"), ("p07_dmE", "--")):
    tt, o = prof(v)
    for (a, b), col in ((("mid", "nea"), "#2ca02c"), (("far", "nea"), "#d62728"), (("far", "mid"), "#7f3fbf")):
        ax[2].plot(tt, o[a] / o[b], color=col, ls=ls_, lw=1.6,
                   label=f"{N[a]}/{N[b]}" + {"-": " 94p07", ":": " 94p10", "--": " 94p07 etter filter E"}[ls_])
shade(ax[2])
ax[2].axhline(1, color="k", lw=0.6)
ax[2].set_ylim(0.45, 1.3)
ax[2].set_title("c) Forhold mellom stakkene (1 = som referansen). Gult: sterk pakke. Grått: myke lag i modellen. Rød strek: der kopien lander",
                fontsize=9, loc="left")
ax[2].set_ylabel("forhold / referanse")
ax[2].set_xlabel("toveis gangtid (ms), 20 ms-vinduer")
ax[2].legend(loc="upper right", ncol=3, fontsize=7.5)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_kontroll_tid.png"), dpi=95)
plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for v, ls_ in (("p07", "-"), ("p07_dmE", "--")):
    d = A[v]
    src = np.nonzero(win(1256, 1356))[0]
    for k in KS:
        a_ = d[f"ac_{k}"].sum(1)
        p = d[f"ss_{k}"].sum(0)
        c = [a_[j, src].sum() / p[src].sum() for j in range(len(ALAGS))]
        ax[0].plot(2 * ALAGS, c, color=C[k], ls=ls_, lw=1.8, marker="o", ms=3, label=N[k] + (" ufiltrert" if ls_ == "-" else " etter filter E"))
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xlabel("forsinkelse (ms)")
ax[0].set_ylabel("kopi av pakken 1256–1356 ms (koeffisient)")
ax[0].set_title("a) Kopien under den sterke pakken, enkelttraser", fontsize=10, loc="left")
ax[0].legend(fontsize=7.5, ncol=2)
ax[0].grid(alpha=0.3)
names = list(WIN2)[:3]
lab = ["over de\nmyke lagene", "ved de\nmyke lagene", "under de\nmyke lagene"]
x = np.arange(3)
for i, (q, col) in enumerate((("mid/nea", "#2ca02c"), ("far/nea", "#d62728"), ("far/mid", "#7f3fbf"))):
    y = [EX["p07"]["vinduer_mot_myke_lag"][n][q]["rms"] for n in names]
    y2 = [EX["p07_dmE"]["vinduer_mot_myke_lag"][n][q]["rms"] for n in names]
    ax[1].bar(x + (i - 1) * 0.27 - 0.06, y, 0.12, color=col, label=q.replace("nea", "near") + " ufiltrert")
    ax[1].bar(x + (i - 1) * 0.27 + 0.06, y2, 0.12, color=col, alpha=0.45, label=q.replace("nea", "near") + " etter filter E")
ax[1].axhline(1, color="k", lw=0.6)
ax[1].set_xticks(x)
ax[1].set_xticklabels([f"{l}\n{EX['vinduer'][n][0]:.0f}–{EX['vinduer'][n][1]:.0f} ms" for l, n in zip(lab, names)], fontsize=8)
ax[1].set_ylim(0.5, 1.1)
ax[1].set_ylabel("forhold / referanse (RMS)")
ax[1].set_title("b) Hvor sitter avvikene? (94p07)", fontsize=10, loc="left")
ax[1].legend(fontsize=7, ncol=3, loc="lower center")
ax[1].grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_kontroll_kopi_og_sted.png"), dpi=95)
