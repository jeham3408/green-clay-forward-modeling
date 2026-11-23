
# ---- fm_sl_fig2.py: figure for the report: (a) the strongest reflector and its sea-bed multiples on near, mid
# and far, (b) far/near per reflector against time, before and after the multiple filter ----
# Run:  cat sl_common.py sl_events_lib.py fm_sl_fig2.py | ssh linux-pc '<geoscience python> -I -'
# Panel (b) uses the numbers printed by fm_sl_robust.py (before) and fm_sl_demult.py (after, filter C).
os.environ["MPLCONFIGDIR"] = os.path.join(AVO, "mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = stacks()
PRE, POST = 30, 150
T0 = 1316
i_c = T0 // 2
tau_s, keep, ccb, lo = track(S["ful"], i_c)
lo2 = i_c - MAXLAG - PRE
mw = {}
for k in ("nea", "mid", "far"):
    W = np.asarray(S[k][:, :, lo2:i_c + MAXLAG + POST + 1], np.float32)
    idx = (MAXLAG + PRE + tau_s)[..., None] + np.arange(-PRE, POST + 1)
    F = np.take_along_axis(W, idx, axis=2).astype(np.float64)
    mw[k] = F[keep].mean(0)
    del W, F
tt = np.arange(-PRE, POST + 1) * 2.0
ref = {k: mw[k][PRE - 2:PRE + 3][np.argmin(mw[k][PRE - 2:PRE + 3])] for k in mw}     # parent is a trough

BEFORE = {1104: (0.54, 0.52, 0.56), 1120: (0.46, 0.44, 0.48), 1140: (0.63, 0.60, 0.65), 1156: (0.64, 0.62, 0.66),
          1168: (0.60, 0.58, 0.62), 1190: (0.48, 0.45, 0.50), 1206: (0.43, 0.41, 0.45), 1228: (0.50, 0.48, 0.52),
          1276: (0.59, 0.57, 0.62), 1316: (0.68, 0.66, 0.69), 1332: (0.64, 0.62, 0.65), 1352: (0.50, 0.48, 0.51),
          1386: (0.44, 0.41, 0.46), 1408: (0.29, 0.28, 0.31), 1422: (0.25, 0.24, 0.27), 1456: (0.27, 0.25, 0.28),
          1472: (0.28, 0.27, 0.30), 1492: (0.31, 0.29, 0.32), 1536: (0.29, 0.27, 0.31), 1584: (0.49, 0.47, 0.52),
          1608: (0.48, 0.45, 0.51), 1644: (0.38, 0.35, 0.41), 1722: (0.51, 0.48, 0.55), 1746: (0.64, 0.61, 0.67),
          1776: (0.61, 0.58, 0.64), 1808: (0.61, 0.58, 0.65), 1828: (0.51, 0.48, 0.54), 1862: (0.58, 0.56, 0.61),
          1900: (0.63, 0.61, 0.66)}
# after the Backus filter of fm_sl_demult3.py (filter D)
AFTER = {1104: (0.57, 0.55, 0.58), 1140: (0.62, 0.60, 0.63), 1156: (0.62, 0.60, 0.64), 1168: (0.60, 0.58, 0.62),
         1190: (0.49, 0.46, 0.52), 1206: (0.45, 0.43, 0.47), 1228: (0.54, 0.52, 0.57), 1276: (0.60, 0.58, 0.62),
         1316: (0.69, 0.67, 0.70), 1332: (0.65, 0.64, 0.67), 1352: (0.47, 0.46, 0.49), 1386: (0.52, 0.48, 0.55),
         1408: (0.35, 0.33, 0.36), 1422: (0.35, 0.33, 0.37), 1456: (0.36, 0.34, 0.37), 1472: (0.36, 0.35, 0.38),
         1492: (0.38, 0.36, 0.40), 1536: (0.50, 0.47, 0.53), 1584: (0.60, 0.57, 0.63), 1608: (0.48, 0.45, 0.51),
         1644: (0.47, 0.43, 0.51), 1722: (0.63, 0.60, 0.67), 1746: (0.72, 0.69, 0.75), 1776: (0.61, 0.58, 0.64),
         1808: (0.65, 0.61, 0.68), 1828: (0.52, 0.50, 0.55), 1862: (0.55, 0.53, 0.57), 1900: (0.65, 0.62, 0.67)}
REF_MED, REF_SD = 0.599, 0.134      # 21 events outside 1380-1560 ms after filter D: median, sd of log
# the 2010 processing of the same survey (94p10), after its own Backus filter (fm_sl10_events.py, filter E)
AFTER10 = {1104: 0.74, 1118: 0.68, 1138: 0.73, 1154: 0.72, 1168: 0.62, 1190: 0.49, 1204: 0.46, 1228: 0.56, 1276: 0.59,
           1316: 0.68, 1332: 0.66, 1352: 0.50, 1386: 0.53, 1408: 0.36, 1422: 0.36, 1454: 0.38, 1470: 0.39, 1492: 0.41,
           1536: 0.49, 1584: 0.62, 1608: 0.46, 1722: 0.64, 1744: 0.72, 1772: 0.59, 1808: 0.59, 1828: 0.55, 1860: 0.53}

INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
COL = {"nea": "#2a78d6", "mid": "#eb6834", "far": "#1baf7a"}
STY = {"nea": "-", "mid": "--", "far": ":"}
NAME = {"nea": "near (3-16°)", "mid": "mid (15-34°)", "far": "far (31-48°)"}
plt.rcParams.update({"font.size": 10, "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})

a = ax[0]
for x0 in (110,):
    a.axvline(x0, color=GRID, lw=6, zorder=0)
a.axhline(0, color=AXIS, lw=1)
for k in ("nea", "mid", "far"):
    a.plot(tt, mw[k] / abs(ref[k]), STY[k], color=COL[k], lw=2, label=NAME[k], solid_capstyle="round")
a.set_xlim(-60, 300)
a.set_ylim(-1.15, 1.05)
a.grid(axis="y", color=GRID, lw=1)
a.set_xlabel("tid etter reflektoren ved 1316 ms (ms)")
a.set_ylabel("amplitude i forhold til reflektoren")
a.annotate("reflektor ved 1316 ms", xy=(1, -0.99), xytext=(28, -0.95), ha="left", va="center", fontsize=9, color=INK2,
           arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
a.text(146, 0.95, "kopi rundt 110 ms under, motsatt\nfortegn. Passer med tiden gjennom\nvannet (81 m), men er bare delvis\nen multippel", ha="left", va="top", fontsize=9, color=INK2)
a.legend(frameon=False, loc="lower right", fontsize=9)
a.set_title("a) Den sterkeste reflektoren og kopien 110 ms under", loc="left", color=INK, fontsize=11)

b = ax[1]
b.axvspan(1395, 1545, color="#f0efec", zorder=0)
b.axhspan(REF_MED * np.exp(-2 * REF_SD), REF_MED * np.exp(2 * REF_SD), color="#e6f1fb", zorder=0)
b.axhline(REF_MED, color="#85b7eb", lw=1.2, ls="--", zorder=0)
b.text(1915, REF_MED * np.exp(2 * REF_SD) + 0.01, "vanlige reflektorer\n(median ± 2 sd)", ha="right", va="bottom", fontsize=8.5, color=INK2)
b.text(1470, 0.84, "sonen der kopiene\nav 1276-1332 ms havner", ha="center", va="top", fontsize=9, color=INK2)
b.plot([1412, 1446], [0.05, 0.05], color=INK2, lw=3, solid_capstyle="butt")
b.text(1429, 0.07, "mykeste lag\ni 15/9-13", ha="center", va="bottom", fontsize=8.5, color=INK2)
tb = np.array(sorted(BEFORE))
vb = np.array([BEFORE[t] for t in tb])
ta = np.array(sorted(AFTER))
va = np.array([AFTER[t] for t in ta])
b.errorbar(tb, vb[:, 0], yerr=[vb[:, 0] - vb[:, 1], vb[:, 2] - vb[:, 0]], fmt="o", ms=6, color=COL["nea"],
           mec="white", mew=1.2, elinewidth=1, capsize=0, label="før multippelfilter")
b.errorbar(ta + 4, va[:, 0], yerr=[va[:, 0] - va[:, 1], va[:, 2] - va[:, 0]], fmt="s", ms=5.5, mfc="white",
           mec=COL["mid"], mew=1.6, color=COL["mid"], elinewidth=1, capsize=0, label="etter multippelfilter")
t10 = np.array(sorted(AFTER10))
b.plot(t10 + 8, [AFTER10[t] for t in t10], "^", ms=6, mfc="white", mec=COL["far"], mew=1.6, ls="none",
       label="etter filter, 2010-prosessering")
b.set_xlim(1080, 1925)
b.set_ylim(0, 0.85)
b.grid(axis="y", color=GRID, lw=1)
b.set_xlabel("toveis gangtid (ms)")
b.set_ylabel("far/near (amplitude)")
hd, lb = b.get_legend_handles_labels()
order = [lb.index(s) for s in ("før multippelfilter", "etter multippelfilter", "etter filter, 2010-prosessering")]
b.legend([hd[i] for i in order], ["2007, før filter", "2007, etter filter", "2010, etter filter"],
         frameon=False, loc="lower left", fontsize=9)
b.set_title("b) Far/near for 29 reflektorer i Hordaland-leira", loc="left", color=INK, fontsize=11)
fig.text(0.01, 0.005, "Sleipner 1994, prosessert 2007 og 2010 (åpne data). Vinklene er utledet fra dataene ved 1,45 s. Strekene i b) er 95 % bootstrap-intervall over blokker på 200 x 225 m.",
         fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
out = os.path.join(AVO, "fig_sleipner_multipler.png")
fig.savefig(out, dpi=150, facecolor="white")
print("written", out)
print("copy ratios read off the normalised mean traces (value at +110 and +220 ms windows):")
for k in ("nea", "mid", "far"):
    y = mw[k] / abs(ref[k])
    i1 = PRE + 50 + int(np.argmin(y[PRE + 50:PRE + 66]))
    i2 = PRE + 104 + int(np.argmin(-y[PRE + 104:PRE + 122]))
    print(f"   {k}: trough at {tt[i1]:.0f} ms = {y[i1]:+.2f}; peak at {tt[i2]:.0f} ms = {y[i2]:+.2f}")
