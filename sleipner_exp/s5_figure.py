"""Figure: near / mid / far stacks with and without the sea surface, and far/near per reflector."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
S = np.load(os.path.join(ROOT, "stacks.npz"))
R = json.load(open(os.path.join(ROOT, "results.json")))
t = 1000 * S["t"]
w = (t >= 1100) & (t <= 1900)
fig, axs = plt.subplots(1, 4, figsize=(15, 8.5), gridspec_kw=dict(width_ratios=[1, 1, 1, 1.25]))
sc = None
for ax, (k, ttl) in zip(axs[:3], (("nea", "near (150–760 m)"), ("mid", "mid (760–1760 m)"), ("far", "far (1760–3150 m)"))):
    g = S["G_trend_main__" + k]
    sc = sc or -1.0 / np.abs(S["G_trend_main__nea"][w]).max()      # minus: the ghosts turn the pulse upside down; shown so that
    #                                                                 an increase in impedance is a positive peak, as in the real data
    ax.axhspan(1408, 1492, color="0.88", zorder=0)
    ax.axhspan(R["checks"]["t0_soft_layer_ms"][0], R["checks"]["t0_soft_layer_ms"][1], color="#f4c095", zorder=0)
    ax.plot(sc * g[w], t[w], "k", lw=1.0, label="uten havflate (bare primærer og spøkelser)")
    if "FS_trend_main__" + k in S.files:
        ax.plot(sc * S["FS_trend_main__" + k][w], t[w], "C3", lw=1.0, label="med havflate (multipler med)")
    ax.set_ylim(1900, 1100)
    ax.set_xlim(-1.6, 1.6)
    ax.set_title(ttl, fontsize=10)
    ax.set_xlabel("amplitude (near uten havflate = 1)")
    ax.grid(alpha=0.3)
axs[0].set_ylabel("toveis gangtid (ms)")
axs[0].legend(fontsize=7, loc="lower left")
ax = axs[3]
ax.axhspan(1408, 1492, color="0.88", zorder=0)
C = R["cases"]["main"]
for name, col, mk, lab in (("G_trend", "k", "o", "uten havflate, Vs på trend"), ("G_m15", "C0", "s", "uten havflate, Vs 15 % lavere"),
                           ("FS_trend", "C3", "o", "med havflate"), ("FS_trend_filt", "C2", "^", "med havflate, etter multippelfilter")):
    if name not in C:
        continue
    key = "events_FS" if name.startswith("FS") and "events_FS" in C[name] else "events_G"
    T = C[name][key]["table"]
    ax.plot([e["fn"] for e in T], [e["t_ms"] for e in T], mk, color=col, ms=5, label=lab, alpha=0.85)
ax.set_xlim(0, 2.0)
ax.set_ylim(1900, 1100)
ax.set_xlabel("far/near per reflektor")
ax.set_title("far/near (grått felt: sonen 1408–1492 ms)", fontsize=10)
ax.legend(fontsize=7, loc="lower right")
ax.grid(alpha=0.3)
fig.suptitle("Syntetiske near-, mid- og far-stakker fra 15/9-13 (elastisk FD, 25 Hz kilde, spøkelser med i begge). Oransje: det myke leirlaget", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "fig_sleipner_syntetisk.png"), dpi=110)
print("figure written")
