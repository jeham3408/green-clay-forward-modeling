
# ---- fm_sl_demult.py: remove the water-layer multiples with a three-point filter per stack and measure again ----
# Run:  cat sl_common.py sl_events_lib.py fm_sl_demult.py | ssh linux-pc '<geoscience python> -I -'
#
# fm_sl_mult.py showed that every stack holds a copy of each strong reflector one water-layer period later (opposite
# sign) and two periods later (same sign). Copy amplitudes a1, a2 relative to the parent (mean of the parents at 1276
# and 1316 ms): near -0.68/+0.35 at 110/220 ms, mid -0.365/+0.28 at 112/224 ms, far -0.28/+0.12 at 119/234 ms.
# The filter y(t) = x(t) + f1*x(t - lag1) + f2*x(t - lag2) with f1 = -a1 and f2 = a1^2 - a2 removes both copies
# (for near this is close to the Backus water-layer operator with sea-bed reflection coefficient 0.34-0.375).
# This is a rough robustness check, not production demultiple.
import sys

S = stacks()
A = {"nea": (55, -0.68, 110, 0.35), "mid": (56, -0.365, 112, 0.28), "far": (59, -0.28, 117, 0.12)}
DM = {}
for k, (l1, a1, l2, a2) in A.items():
    x = np.asarray(S[k], np.float32)
    y = x.copy()
    y[:, :, l1:] += np.float32(-a1) * x[:, :, :-l1]
    y[:, :, l2:] += np.float32(a1 * a1 - a2) * x[:, :, :-l2]
    DM[k] = y
    del x
    print(f"filter {k}: + {-a1:.3f} x(t-{l1*2} ms) + {a1*a1-a2:+.3f} x(t-{l2*2} ms)")
PRE, POST = 30, 150


def flat_long(arr, i_c, tau_s):
    lo = i_c - MAXLAG - PRE
    W = np.asarray(arr[:, :, lo:i_c + MAXLAG + POST + 1], np.float32)
    idx = (MAXLAG + PRE + tau_s)[..., None] + np.arange(-PRE, POST + 1)
    return np.take_along_axis(W, idx, axis=2).astype(np.float64)


print("\n== copies of the parent before and after the filter (amplitude relative to the parent, at the lag found before)")
for T0 in (1316, 1276, 1140, 1746):
    i_c = T0 // 2
    if i_c + MAXLAG + POST + 1 > NS:
        continue
    tau_s, keep, ccb, lo = track(S["ful"], i_c)
    row = []
    for k in ("nea", "mid", "far"):
        l1, a1, l2, a2 = A[k]
        res = []
        for arr in (S[k], DM[k]):
            F = flat_long(arr, i_c, tau_s)
            mw = F[keep].mean(0)
            par0 = flat_long(S[k], i_c, tau_s)[keep].mean(0)[PRE - 10:PRE + 21]
            pn = (par0 ** 2).sum()
            c1 = max((((mw[PRE - 10 + L:PRE + 21 + L] * par0).sum() / pn) for L in range(l1 - 3, l1 + 4)), key=abs)
            c2 = max((((mw[PRE - 10 + L:PRE + 21 + L] * par0).sum() / pn) for L in range(l2 - 4, l2 + 5)), key=abs)
            res.append((c1, c2))
        row.append(f"{k}: first copy {res[0][0]:+.2f} -> {res[1][0]:+.2f}, second copy {res[0][1]:+.2f} -> {res[1][1]:+.2f}")
    print(f"  parent {T0}: " + " | ".join(row))

EV = [1104, 1140, 1190, 1228, 1276, 1316, 1332, 1386, 1408, 1422, 1456, 1472, 1492, 1536, 1584, 1608, 1722, 1746, 1776, 1808, 1828, 1900]
print("\n== events before and after the multiple filter (same horizons, tracked on the original full stack)")
print("  time | before: far/near  mid/near  near amp | after: far/near (95 %)      mid/near (95 %)      near amp | near amp kept")
S3 = {k: S[k] for k in ("nea", "mid", "far")}
OUT = {}
for T0 in EV:
    hz = track(S["ful"], T0 // 2)
    m0 = measure(S3, S["ful"], T0, horizon=hz)
    m1 = measure(DM, S["ful"], T0, horizon=hz)
    c = H
    a0 = m0["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m0["mw"]["nea"][c - 2:c + 3]))]
    a1_ = m1["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m1["mw"]["nea"][c - 2:c + 3]))]
    OUT[T0] = (m0, m1, a0, a1_)
    print(f"  {T0:4d} | {m0['fn'][0]:14.2f} {m0['mn'][0]:9.2f} {a0:+9.3f} | {m1['fn'][0]:12.2f} ({m1['fn'][1]:.2f}-{m1['fn'][2]:.2f})   {m1['mn'][0]:5.2f} ({m1['mn'][1]:.2f}-{m1['mn'][2]:.2f})  {a1_:+8.3f} | {a1_/a0:6.2f}")
    sys.stdout.flush()

clean = [T for T in EV if abs(OUT[T][3] / OUT[T][2]) > 0.7 and abs(OUT[T][3]) > 0.08]
print("\n== events that keep more than 70 % of their near amplitude after the filter (taken as primaries):", clean)
fnv = np.array([OUT[T][1]["fn"][0] for T in clean])
mnv = np.array([OUT[T][1]["mn"][0] for T in clean])
print(f"   far/near: median {np.median(fnv):.2f}, min {fnv.min():.2f}, max {fnv.max():.2f}, spread (sd of log) {np.std(np.log(fnv)):.2f}")
print(f"   mid/near: median {np.median(mnv):.2f}, min {mnv.min():.2f}, max {mnv.max():.2f}, spread (sd of log) {np.std(np.log(mnv)):.2f}")

# well tie before and after
w, t, L = well()
SYN = np.load(os.path.join(AVO, "synth_trend.npy"))
print("\n== well tie at 15/9-13 before and after the filter (trend-Vs synthetic, lags +-40 ms)")
for win in ((1270, 1560), (1560, 1960)):
    i0, i1 = win[0] // 2, win[1] // 2
    for j, k in enumerate(("nea", "mid", "far")):
        res = []
        for arr in (S[k], DM[k]):
            real = np.asarray(arr[WELL_IL - 2:WELL_IL + 3, WELL_XL - 2:WELL_XL + 3, :], np.float64).mean((0, 1))
            res.append(xc(SYN[j][i0:i1], real[i0:i1]))
        print(f"  {win[0]}-{win[1]} {k}: before cc0 {res[0][0]:+.2f}, best {res[0][1]:+.2f} at {res[0][2]:+.0f} ms | after cc0 {res[1][0]:+.2f}, best {res[1][1]:+.2f} at {res[1][2]:+.0f} ms")

# share of the energy removed by the filter per window
print("\n== rms after / before the filter, per stack")
for win in ((1100, 1380), (1380, 1560), (1560, 1960)):
    i0, i1 = win[0] // 2, win[1] // 2
    print(f"  {win[0]}-{win[1]}: " + ", ".join(f"{k} {np.sqrt(np.mean(DM[k][::3, ::3, i0:i1].astype(np.float64)**2) / np.mean(np.asarray(S[k][::3, ::3, i0:i1], np.float64)**2)):.2f}" for k in ("nea", "mid", "far")))
for k in DM:
    np.save(os.path.join(AVO, f"dm_{k}.npy"), DM[k])
