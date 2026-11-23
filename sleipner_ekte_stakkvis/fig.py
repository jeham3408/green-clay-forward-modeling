#!/usr/bin/env python3
"""fig.py: the two figures, from results.json and profiles_*.npz only."""
import json
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_stakkvis")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis"
R = json.load(open(D + "/results.json"))
S = R["oppsummering"]
KS = ("nea", "mid", "far")
COL = {"nea": "#1f5fa8", "mid": "#2e8b57", "far": "#c0392b"}
LAB = {"nea": "near", "mid": "mid", "far": "far"}
tms = 800.0 + 2.0 * np.arange(600)
RM = ((tms >= 1100) & (tms < 1380)) | ((tms > 1560) & (tms <= 1900))
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def slide(v):
    P = np.load(f"{D}/profiles_{v}.npz")
    out = {}
    for k in KS:
        e = P[f"e2t_{k}_kube"]
        out[k] = np.sqrt(np.convolve(e, np.ones(43) / 43, mode="same") / e[RM].mean())
    return out


def deco(ax, ymax):
    ax.axvspan(1408, 1492, color="0.85", zorder=0)
    ax.axvspan(1256, 1356, color="#f6e7c1", zorder=0)
    ax.plot([1412, 1446], [0.04 * ymax + ax.get_ylim()[0]] * 2, color="k", lw=4, solid_capstyle="butt")
    ax.axhline(1, color="0.5", lw=0.8)
    ax.set_xlim(1100, 1900)


fig, ax = plt.subplots(3, 1, figsize=(9, 9.6), sharex=True)
for a, (v1, v2, ttl) in zip(ax[:2], (("p07", "p10", "a) Ufiltrert: styrken til hver stakk (RMS i glidende 84 ms-vindu, delt på referansen)"),
                                     ("p07_dmE", "p10_dmE", "b) Det samme etter multippelfilter E"))):
    s1, s2 = slide(v1), slide(v2)
    for k in KS:
        a.plot(tms, s1[k], color=COL[k], lw=1.8, label=f"{LAB[k]}, 94p07")
        a.plot(tms, s2[k], color=COL[k], lw=1.2, ls="--", label=f"{LAB[k]}, 94p10")
    a.set_ylim(0.5, 1.75)
    deco(a, 1.25)
    a.set_title(ttl, loc="left", fontsize=10.5)
    a.set_ylabel("styrke / referanse")
ax[0].legend(ncol=3, fontsize=8.5, frameon=False, loc="upper right")
ax[0].text(1450, 1.66, "sonen", ha="center", fontsize=9)
ax[0].text(1306, 1.66, "sterk pakke", ha="center", fontsize=9)
a = ax[2]
for v, ls, lw, nm in (("p07", "-", 1.8, "ufiltrert"), ("p07_dmE", "--", 1.5, "etter filter E")):
    s = slide(v)
    a.plot(tms, s["mid"] / s["nea"], color=COL["mid"], ls=ls, lw=lw, label=f"mid/near, {nm}")
    a.plot(tms, s["far"] / s["nea"], color=COL["far"], ls=ls, lw=lw, label=f"far/near, {nm}")
a.set_ylim(0.5, 1.3)
deco(a, 0.8)
a.set_title("c) Forholdet mellom stakkene, 94p07 (delt på forholdet i referansen)", loc="left", fontsize=10.5)
a.set_ylabel("forhold / referanse")
a.set_xlabel("toveistid (ms). Svart strek: det mykeste leirlaget i 15/9-13 (1412–1446 ms)")
a.legend(ncol=2, fontsize=8.5, frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(D + "/fig_stakkvis_tid.png", dpi=110)
plt.close(fig)

# figure 2: the pattern against the two explanations
fig, a = plt.subplots(figsize=(9.2, 5.6))
cats = ["near", "mid", "far", "mid/near", "far/near"]
filt = ["p07_dmE", "p10_dmE", "p07_dm", "p07_dmD", "p07_dmA", "p07_dmB"]


def vals(v, i):
    if i < 3:
        n = S[v]["nivaa"][KS[i]]
        return n["rms"], [n[q] for q in ("rms", "abs", "env", "coh", "env_topp_reflektor")]
    d = S[v]["forhold"]["mid_near" if i == 3 else "far_near"]
    return d["vindu_rms"], [d["amplitudemaal_min"], d["amplitudemaal_max"]]


for i in range(5):
    for v, dx, mk, lab in (("p07", -0.30, "o", "målt, 94p07 ufiltrert"), ("p10", -0.20, "s", "målt, 94p10 ufiltrert")):
        c, rng = vals(v, i)
        a.plot([i + dx] * 2, [min(rng), max(rng)], color="k", lw=1.5)
        a.plot(i + dx, c, mk, color="k", ms=7, label=lab if i == 0 else None)
    cs, rs = zip(*[vals(v, i) for v in filt])
    allr = [x for r in rs for x in r]
    a.plot([i - 0.08] * 2, [min(allr), max(allr)], color="0.45", lw=1.5)
    a.plot(i - 0.08, np.median(cs), "o", mfc="w", mec="0.3", ms=7, label="målt, etter multippelfilter (seks varianter)" if i == 0 else None)
P = S["prediksjoner"]["sett"]
sets = list(P)
CASES = (("ingen_forklaring", "0.6", "modell: ingen av forklaringene"), ("b_multipler_skalert", "#d98200", "modell: (b) multipler, skalert til målt kopistyrke"),
         ("a_myke_lag_vs07", "#a98bd6", "modell: (a) Vs 7 % lavere i de myke lagene"), ("a_myke_lag_vs15", "#5b2a86", "modell: (a) Vs 15 % lavere i de myke lagene"))
QK = ("near", "mid", "far", "mid_near", "far_near")
for j, (case, c, lab) in enumerate(CASES):
    for i, q in enumerate(QK):
        ys = [P[s_]["forklaringer"][case][q] for s_ in sets]
        x = i + 0.10 + 0.11 * j
        a.plot([x, x], [min(ys), max(ys)], color=c, lw=5, solid_capstyle="butt", label=lab if i == 0 else None)
        a.plot([x] * len(ys), ys, "_", color=c, ms=11, mew=2.5)
a.axhline(1, color="0.5", lw=0.8)
a.axvline(2.68, color="0.8", lw=0.8)
a.set_xticks(range(5))
a.set_xticklabels(["near\n(styrke sone/ref.)", "mid\n(styrke sone/ref.)", "far\n(styrke sone/ref.)", "mid/near\n(sone/ref.)", "far/near\n(sone/ref.)"])
a.set_ylim(0.3, 1.65)
a.set_xlim(-0.6, 4.75)
a.set_ylabel("sonen 1408–1492 ms delt på referansen")
a.set_title("Mønsteret over near, mid og far: målt (hele kuben) mot de to forklaringene.\nFargede streker: modellen med to valg av offsetområder. Streker på målingene: spenn over amplitudemålene.", loc="left", fontsize=10)
a.legend(fontsize=8, frameon=False, loc="lower left", ncol=2)
fig.tight_layout()
fig.savefig(D + "/fig_stakkvis_monster.png", dpi=110)
print("figs ok", os.path.getsize(D + "/fig_stakkvis_tid.png"), os.path.getsize(D + "/fig_stakkvis_monster.png"))
