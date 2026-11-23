#!/usr/bin/env python3
"""Two figures from maps.npz and results.json (s2_analyse.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig"
M = np.load(OUT + "/maps.npz")
R = json.load(open(OUT + "/results.json"))
ext = [0, 26 * 0.225, 16 * 0.2, 0]
WELL = (143 * 0.0125, 57 * 0.0125)

fig, ax = plt.subplots(2, 3, figsize=(13, 6.2))
ok = M["p07_ok"]
panels = [("p07_t1316", "Sterk reflektor ~1316 ms: tid (ms)", "viridis", None),
          ("iso86", "Avstand 1316 -> 1386 ms (ms)", "viridis", None),
          ("p07_copy", "Kopi paa near: koeffisient ved ~110 ms", "magma", (-1.5, -0.3)),
          ("p07_D_fn", "far/near i sonen delt paa referansen", "RdBu", (0.2, 1.8)),
          ("p07_L_nea", "near: sone / referanse", "cividis", None),
          ("p07_L_far", "far: sone / referanse", "cividis", None)]
for a, (k, title, cmap, lim) in zip(ax.ravel(), panels):
    v = (M["p07_t1386"] - M["p07_t1316"]) if k == "iso86" else M[k]
    v = np.where(ok, v, np.nan)
    if lim is None:
        lim = (np.nanpercentile(v, 2), np.nanpercentile(v, 98))
    im = a.imshow(v, extent=ext, cmap=cmap, vmin=lim[0], vmax=lim[1], aspect="equal", interpolation="nearest")
    a.plot(*WELL, "k*", ms=9, mec="w")
    a.set_title(title, fontsize=10)
    a.set_xlabel("km (crossline)", fontsize=8)
    a.set_ylabel("km (inline)", fontsize=8)
    plt.colorbar(im, ax=a, shrink=0.75)
fig.suptitle("Sleipner 94p07, blokker 200 x 225 m. Stjerne: 15/9-13", fontsize=11)
fig.tight_layout()
fig.savefig(OUT + "/fig_romlig_kart.png", dpi=75)

fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
for p, col in (("p07", "C0"), ("p10", "C3")):
    okp = M[p + "_ok"]
    x, y = -M[p + "_copy"][okp], M[p + "_D_fn_refbelow"][okp]
    ax[0].loglog(x, y, ".", color=col, ms=4, alpha=0.5, label=f"{p}: r = {R[p]['test2']['far_near_refbelow_vs_copy']['pearson_r']:.2f}")
    qs = np.percentile(x, np.linspace(0, 100, 8))
    xm = [np.median(x[(x >= a) & (x <= b)]) for a, b in zip(qs[:-1], qs[1:])]
    ym = [np.exp(np.median(np.log(y[(x >= a) & (x <= b)]))) for a, b in zip(qs[:-1], qs[1:])]
    ax[0].plot(xm, ym, "-o", color=col, lw=2, ms=5, mec="k")
    t = M[p + "_prof_tau_ms"]
    ax[1].plot(t, np.exp(M[p + "_prof_all"]), color=col, lw=2, label=f"far/near {p}")
    ax[1].plot(t, np.exp(M[p + "_prof_mid_all"]), color=col, lw=1.2, ls="--", label=f"mid/near {p}")
ax[0].axhline(1, color="k", lw=0.8)
ax[0].set_xlim(0.25, 2.5)
ax[0].set_ylim(0.12, 3.0)
ax[0].set_xlabel("styrke paa kopien paa near (|koeffisient|)")
ax[0].set_ylabel("far/near i sonen / referansen under")
ax[0].set_title("Test 2: underskudd mot kopistyrke, blokk for blokk", fontsize=10)
ax[0].legend(fontsize=8)
ax[1].axvspan(84, 184, color="0.85")
for v in (-40, 0):
    ax[1].axvline(v, color="k", lw=0.8, ls=":")
ax[1].axvline(110, color="k", lw=0.8)
ax[1].set_xlabel("ms etter reflektoren ved 1316 ms (flatet)")
ax[1].set_ylabel("amplitudeforhold (vindu +-24 ms)")
ax[1].set_title("Forhold mot tid. Graatt: sonen. Strek: +110 ms", fontsize=10)
ax[1].legend(fontsize=8)
g = R["p07"]["test2"]["copy_strength_groups_with_ci"]
g10 = R["p10"]["test2"]["copy_strength_groups_with_ci"]
labs = ["svakeste\ntidel", "svak\ntredel", "sterk\ntredel"]
keys = ["weakest_decile", "weak_tercile", "strong_tercile"]
w = 0.13
for j, (nm, col) in enumerate((("level_near_refbelow", "C2"), ("level_mid_refbelow", "C1"), ("level_far_refbelow", "C4"))):
    for q, (gg, hatch) in enumerate(((g, None), (g10, "//"))):
        vals = [gg[k][nm]["value"] for k in keys]
        lo = [gg[k][nm]["value"] - gg[k][nm]["ci95"][0] for k in keys]
        hi = [gg[k][nm]["ci95"][1] - gg[k][nm]["value"] for k in keys]
        ax[2].bar(np.arange(3) + (2 * j + q - 2.5) * w, vals, w, yerr=[lo, hi], color=col, hatch=hatch, edgecolor="k", lw=0.5,
                  label=(nm.split("_")[1] + (" 2007" if q == 0 else " 2010")), error_kw=dict(lw=0.8))
ax[2].set_xticks(range(3))
ax[2].set_xticklabels(labs)
ax[2].set_ylabel("amplitude i sonen / referansen under")
ax[2].set_title("Near, mid og far hver for seg, etter kopistyrke", fontsize=10)
ax[2].set_ylim(0, 1.65)
ax[2].legend(fontsize=7, ncol=3, loc="upper left")
fig.tight_layout()
fig.savefig(OUT + "/fig_romlig_tester.png", dpi=75)
print("figs ok")
