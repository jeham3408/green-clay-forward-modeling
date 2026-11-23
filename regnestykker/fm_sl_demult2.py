
# ---- fm_sl_demult2.py: does what is left after the multiple filter depend on the filter? ----
# Run:  cat sl_common.py sl_events_lib.py fm_sl_demult2.py | ssh linux-pc '<geoscience python> -I -u -'
# Three versions of the three-point filter: coefficients from the parent at 1316 ms (A), from the parent at
# 1276 ms (B), and their mean (C, made by fm_sl_demult.py). Low memory: filtered cubes are written as memmaps.
import sys

S = stacks()
SETS = {
    "A": {"nea": (55, -0.75, 110, 0.38), "mid": (56, -0.36, 112, 0.32), "far": (60, -0.27, 117, 0.16)},
    "B": {"nea": (55, -0.60, 110, 0.31), "mid": (56, -0.37, 112, 0.24), "far": (59, -0.29, 118, 0.07)},
}
DMS = {"C": {k: np.load(os.path.join(AVO, f"dm_{k}.npy"), mmap_mode="r") for k in ("nea", "mid", "far")}}
for name, A in SETS.items():
    DMS[name] = {}
    for k, (l1, a1, l2, a2) in A.items():
        p = os.path.join(AVO, f"dm{name}_{k}.npy")
        if not os.path.exists(p):
            out = np.lib.format.open_memmap(p, mode="w+", dtype=np.float32, shape=S[k].shape)
            for il in range(S[k].shape[0]):
                x = np.asarray(S[k][il], np.float32)
                y = x.copy()
                y[:, l1:] += np.float32(-a1) * x[:, :-l1]
                y[:, l2:] += np.float32(a1 * a1 - a2) * x[:, :-l2]
                out[il] = y
            out.flush()
            del out
        DMS[name][k] = np.load(p, mmap_mode="r")
print("filters:  A from the 1316 parent, B from the 1276 parent, C their mean")

EV = [1276, 1316, 1332, 1386, 1408, 1422, 1456, 1472, 1492, 1536, 1584, 1746, 1828, 1900]
S3 = {k: S[k] for k in ("nea", "mid", "far")}
HOR = {}
print("\n  time | unfiltered far/near mid/near near | A: far/near mid/near near | B: far/near mid/near near | C: far/near mid/near near")
R = {}
for T0 in EV:
    hz = track(S["ful"], T0 // 2)
    HOR[T0] = hz
    row = []
    R[T0] = {}
    for name, src in (("0", S3), ("A", DMS["A"]), ("B", DMS["B"]), ("C", DMS["C"])):
        m = measure(src, S["ful"], T0, horizon=hz)
        c = H
        a = m["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m["mw"]["nea"][c - 2:c + 3]))]
        R[T0][name] = (m["fn"][0], m["mn"][0], a)
        row.append(f"{m['fn'][0]:5.2f} {m['mn'][0]:5.2f} {a:+.3f}")
    print(f"  {T0:4d} | " + " | ".join(row))
    sys.stdout.flush()

ref = [1276, 1316, 1332, 1746, 1828, 1900]
zone = [1408, 1422, 1456, 1472, 1492]
print("\n== zone 1408-1492 ms against the reference reflectors (1276, 1316, 1332, 1746, 1828, 1900): ratio of medians")
for name in ("0", "A", "B", "C"):
    fz, fr = np.median([R[t][name][0] for t in zone]), np.median([R[t][name][0] for t in ref])
    mz, mr = np.median([R[t][name][1] for t in zone]), np.median([R[t][name][1] for t in ref])
    print(f"   filter {name}: far/near zone {fz:.2f}, reference {fr:.2f}, zone/reference {fz/fr:.2f}   |   mid/near zone {mz:.2f}, reference {mr:.2f}, zone/reference {mz/mr:.2f}")

print("\n== do the zone events follow the structure of the reflector one water-layer period above? Correlation of horizon-time maps (kept traces, block medians)")


def bmap(T0):
    tau_s, keep, ccb, lo = HOR[T0]
    s = block_sum(tau_s * keep.astype(float))
    n = block_sum(keep.astype(float))
    return np.where(n > 0.3 * BS[0] * BS[1], s / np.maximum(n, 1), np.nan)


pairs = [(1386, 1276), (1422, 1316), (1456, 1316), (1456, 1332), (1536, 1316), (1492, 1386), (1746, 1316), (1828, 1316), (1900, 1316), (1332, 1316), (1276, 1316)]
for a, b in pairs:
    x, y = bmap(a), bmap(b)
    ok = np.isfinite(x) & np.isfinite(y)
    cc = np.corrcoef(x[ok], y[ok])[0, 1]
    sl = np.polyfit(y[ok], x[ok], 1)[0]
    print(f"   {a} against {b}: correlation {cc:+.2f}, slope {sl:+.2f}, relief (sd) {np.std(x[ok])*2:.1f} ms against {np.std(y[ok])*2:.1f} ms, blocks {ok.sum()}")
