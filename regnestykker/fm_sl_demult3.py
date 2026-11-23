
# ---- fm_sl_demult3.py: multiple filter with physically consistent coefficients, and the reference spread
# computed without the zone ----
# Run:  cat sl_common.py sl_events_lib.py fm_sl_demult3.py | ssh linux-pc '<geoscience python> -I -u -'
# Water-layer reverberation on both sides of the spread: 1/(1 + R z)^2 = 1 - 2R z + 3R^2 z^2 - ...
# The exact inverse (Backus) is (1 + R z)^2 = 1 + 2R z + R^2 z^2, so f1 = 2R = -a1 and f2 = R^2 = a1^2 / 4.
# fm_sl_demult.py used f2 = a1^2 - a2 with the measured second copy a2. On mid and far that copy also holds
# primary energy (verification), which made f2 negative and removed primaries 224-234 ms below strong events.
import sys

S = stacks()
A1 = {"nea": (55, -0.68), "mid": (56, -0.365), "far": (59, -0.28)}      # lag (samples), first copy (mean of two parents)
DM = {}
for k, (l1, a1) in A1.items():
    l2 = 2 * l1
    f1, f2 = -a1, a1 * a1 / 4
    p = os.path.join(AVO, f"dmD_{k}.npy")
    if not os.path.exists(p):
        out = np.lib.format.open_memmap(p, mode="w+", dtype=np.float32, shape=S[k].shape)
        for il in range(S[k].shape[0]):
            x = np.asarray(S[k][il], np.float32)
            y = x.copy()
            y[:, l1:] += np.float32(f1) * x[:, :-l1]
            y[:, l2:] += np.float32(f2) * x[:, :-l2]
            out[il] = y
        out.flush()
        del out
    DM[k] = np.load(p, mmap_mode="r")
    print(f"filter D {k}: + {f1:.3f} x(t-{l1*2} ms) + {f2:.4f} x(t-{l2*2} ms)")

EV = [1104, 1140, 1156, 1168, 1190, 1206, 1228, 1276, 1316, 1332, 1352, 1386, 1408, 1422, 1456, 1472, 1492, 1536,
      1584, 1608, 1644, 1722, 1746, 1776, 1808, 1828, 1862, 1900]
ZONE = (1400, 1500)
S3 = {k: S[k] for k in ("nea", "mid", "far")}
R = {}
print("\n  time | unfiltered: far/near mid/near near | filter D: far/near (95 %)      mid/near (95 %)     near | kept")
for T0 in EV:
    hz = track(S["ful"], T0 // 2)
    m0 = measure(S3, S["ful"], T0, horizon=hz)
    m1 = measure(DM, S["ful"], T0, horizon=hz)
    c = H
    a0 = m0["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m0["mw"]["nea"][c - 2:c + 3]))]
    a1_ = m1["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m1["mw"]["nea"][c - 2:c + 3]))]
    R[T0] = (m0["fn"][0], m0["mn"][0], m1["fn"], m1["mn"], a0, a1_)
    print(f"  {T0:4d} | {m0['fn'][0]:9.2f} {m0['mn'][0]:8.2f} {a0:+.3f} | {m1['fn'][0]:8.2f} ({m1['fn'][1]:.2f}-{m1['fn'][2]:.2f})   {m1['mn'][0]:5.2f} ({m1['mn'][1]:.2f}-{m1['mn'][2]:.2f})  {a1_:+.3f} | {a1_/a0:5.2f}")
    sys.stdout.flush()

inz = [t for t in EV if ZONE[0] <= t <= ZONE[1]]
out = [t for t in EV if not (1380 <= t <= 1560)]
for lab, idx in (("unfiltered", 0), ("filter D", 2)):
    fo = np.array([R[t][idx][0] if idx else R[t][0] for t in out])
    fz = np.array([R[t][idx][0] if idx else R[t][0] for t in inz])
    mo = np.array([R[t][3][0] if idx else R[t][1] for t in out])
    mz = np.array([R[t][3][0] if idx else R[t][1] for t in inz])
    sd = np.std(np.log(fo), ddof=1)
    zsc = (np.log(np.median(fz)) - np.mean(np.log(fo))) / sd
    print(f"\n== {lab}: reference = {len(out)} events outside 1380-1560 ms; zone = {inz}")
    print(f"   far/near reference median {np.median(fo):.3f} (range {fo.min():.2f}-{fo.max():.2f}, sd of log {sd:.3f}); zone median {np.median(fz):.3f} -> zone/reference {np.median(fz)/np.median(fo):.2f}, {zsc:+.1f} sd")
    print(f"   mid/near reference median {np.median(mo):.3f} (sd of log {np.std(np.log(mo), ddof=1):.3f}); zone median {np.median(mz):.3f} -> zone/reference {np.median(mz)/np.median(mo):.2f}")
