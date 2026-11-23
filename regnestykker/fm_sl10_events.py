
# ---- fm_sl10_events.py: the same reflectors in two processings of the same 1994 survey (2007: 94p07, 2010: 94p10).
# If the far deficit at the soft clay layer is in the rock, it should be there in both. ----
# Run:  cat sl_common.py sl_events_lib.py fm_sl10_events.py | ssh linux-pc '<geoscience python> -I -u -'
import sys

DIRS = {"07": R + "/sleipner/avo", "10": R + "/sleipner/avo10"}
ST = {tag: {k: np.load(os.path.join(d, k + ".npy"), mmap_mode="r") for k in KEYS} for tag, d in DIRS.items()}
COV = np.load(os.path.join(DIRS["10"], "coverage.npy"))
print("94p10 covers", int(COV.sum()), "of", COV.size, "94p07 nodes")
PRE, POST = 30, 150

# ---------- 1. fold weights in 94p10
print("\n== full ~ a*near + b*mid + c*far, 94p10 (100 ms windows)")
sub10 = {k: np.asarray(ST["10"][k][::3, ::3, :1000], np.float64)[COV[::3, ::3]] for k in KEYS}
for c0 in (600, 900, 1100, 1300, 1500, 1700, 1900):
    i0, i1 = (c0 - 50) // 2, (c0 + 50) // 2
    X = np.stack([sub10[k][:, i0:i1].ravel() for k in ("nea", "mid", "far")], 1)
    y = sub10["ful"][:, i0:i1].ravel()
    co, *_ = np.linalg.lstsq(X, y, rcond=None)
    r2 = 1 - np.sum((y - X @ co) ** 2) / np.sum(y ** 2)
    print(f"   {c0} ms: a {co[0]:.3f} b {co[1]:.3f} c {co[2]:.3f} (sum {co.sum():.3f}), R2 {r2:.3f}")
del sub10

# ---------- 2. events from the cube-wide mean full trace of each processing
EV = {}
for tag in ("07", "10"):
    mt = np.asarray(ST[tag]["ful"][::2, ::2, :1000], np.float64)
    if tag == "10":
        mt = mt[COV[::2, ::2]]
    else:
        mt = mt.reshape(-1, mt.shape[-1])
    mt = mt.mean(0)
    r1 = np.sqrt(np.mean(mt[540:985] ** 2))
    ev = []
    for i in range(541, 951):
        v = mt[i]
        if abs(v) > 0.8 * r1 and abs(v) >= abs(mt[i - 1]) and abs(v) >= abs(mt[i + 1]):
            if not ev or i * 2 - ev[-1] >= 12:
                ev.append(i * 2)
    EV[tag] = ev
    print(f"events {tag}: {ev}")
pairs = []
for t7 in EV["07"]:
    t10 = min(EV["10"], key=lambda e: abs(e - t7))
    if abs(t10 - t7) <= 6:
        pairs.append((t7, t10))
print("matched events (94p07 time, 94p10 time):", pairs)


# ---------- 3. copies 110 ms below the strongest reflectors, in each processing
def flat_long(arr, i_c, tau_s):
    lo = i_c - MAXLAG - PRE
    W = np.asarray(arr[:, :, lo:i_c + MAXLAG + POST + 1], np.float32)
    idx = (MAXLAG + PRE + tau_s)[..., None] + np.arange(-PRE, POST + 1)
    return np.take_along_axis(W, idx, axis=2).astype(np.float64)


A1, LAG = {}, {}
print("\n== copy one water-layer period below the parent (LS amplitude of the parent wavelet, best lag 100-125 ms)")
for tag in ("07", "10"):
    for T0 in (1316, 1276):
        Tp = min(EV[tag], key=lambda e: abs(e - T0))
        tau_s, keep, ccb, lo = track(ST[tag]["ful"], Tp // 2)
        row = []
        for k in ("nea", "mid", "far"):
            F = flat_long(ST[tag][k], Tp // 2, tau_s)
            mw = F[keep].mean(0)
            par = mw[PRE - 7:PRE + 8]
            pn = (par ** 2).sum()
            best = min(((mw[PRE - 7 + L:PRE + 8 + L] * par).sum() / pn, L) for L in range(50, 63))
            A1.setdefault(tag, {}).setdefault(k, []).append(best[0])
            LAG.setdefault(tag, {}).setdefault(k, []).append(best[1])
            row.append(f"{k} {best[0]:+.2f} at {best[1]*2} ms")
        print(f"   {tag} parent {Tp}: " + "; ".join(row))
for tag in A1:
    for k in A1[tag]:
        A1[tag][k] = float(np.mean(A1[tag][k]))
        LAG[tag][k] = int(round(np.mean(LAG[tag][k])))
print("mean first-copy ratios and lags (samples) used for the filter:", {t: {k: (round(A1[t][k], 3), LAG[t][k]) for k in A1[t]} for t in A1})


# ---------- 4. Backus filter per processing (f1 = -a1, f2 = a1^2/4 at twice the lag)
def backus(tag):
    out = {}
    for k in ("nea", "mid", "far"):
        p = os.path.join(DIRS[tag], f"dmE_{k}.npy")
        l1 = LAG[tag][k]
        a1 = A1[tag][k]
        if not os.path.exists(p):
            src = ST[tag][k]
            o = np.lib.format.open_memmap(p, mode="w+", dtype=np.float32, shape=src.shape)
            for il in range(src.shape[0]):
                x = np.asarray(src[il], np.float32)
                y = x.copy()
                y[:, l1:] += np.float32(-a1) * x[:, :-l1]
                y[:, 2 * l1:] += np.float32(a1 * a1 / 4) * x[:, :-2 * l1]
                o[il] = y
            o.flush()
            del o
        out[k] = np.load(p, mmap_mode="r")
    return out


DM = {tag: backus(tag) for tag in ("07", "10")}

# ---------- 5. ratios per event, both processings, unfiltered and filtered, two estimators
print("\n== far/near and mid/near per event. IV = third stack as instrument; PK = block-mean peak ratio (no alignment between stacks)")
print(" t07  t10 | 94p07 IV f/n m/n  PK f/n | 94p07 filt IV f/n m/n | 94p10 IV f/n m/n  PK f/n | 94p10 filt IV f/n m/n | kept 07 10")
RES = {}
for t7, t10 in pairs:
    r = {}
    for tag, T in (("07", t7), ("10", t10)):
        hz = track(ST[tag]["ful"], T // 2)
        S3 = {k: ST[tag][k] for k in ("nea", "mid", "far")}
        m0 = measure(S3, ST[tag]["ful"], T, horizon=hz)
        m1 = measure(DM[tag], ST[tag]["ful"], T, horizon=hz)
        c = H
        a0 = m0["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m0["mw"]["nea"][c - 2:c + 3]))]
        a1_ = m1["mw"]["nea"][c - 2:c + 3][np.argmax(np.abs(m1["mw"]["nea"][c - 2:c + 3]))]
        r[tag] = (m0["fn"][0], m0["mn"][0], m0["pk_fn"][0], m1["fn"][0], m1["mn"][0], a1_ / a0, m0["kept"])
    RES[t7] = r
    a, b = r["07"], r["10"]
    print(f" {t7} {t10} | {a[0]:9.2f} {a[1]:4.2f} {a[2]:6.2f} | {a[3]:15.2f} {a[4]:4.2f} | {b[0]:9.2f} {b[1]:4.2f} {b[2]:6.2f} | {b[3]:15.2f} {b[4]:4.2f} | {a[5]:4.2f} {b[5]:4.2f}")
    sys.stdout.flush()


def zone_stats(tag, idx, label):
    zone = [t for t in RES if 1400 <= t <= 1500]
    ref = [t for t in RES if not (1380 <= t <= 1560)]
    fz = np.array([RES[t][tag][idx] for t in zone])
    fr = np.array([RES[t][tag][idx] for t in ref])
    ok = (fr > 0) & np.isfinite(fr)
    sd = np.std(np.log(fr[ok]), ddof=1)
    zs = (np.log(np.median(fz)) - np.mean(np.log(fr[ok]))) / sd
    print(f"   {tag} {label}: zone {zone} median {np.median(fz):.3f}; reference ({len(ref)} events) median {np.median(fr):.3f}, sd(log) {sd:.3f}; "
          f"zone/reference {np.median(fz)/np.median(fr):.2f} ({zs:+.1f} sd)")


print("\n== zone 1400-1500 ms against reflectors outside 1380-1560 ms")
for tag in ("07", "10"):
    zone_stats(tag, 0, "far/near IV, unfiltered")
    zone_stats(tag, 2, "far/near PK, unfiltered")
    zone_stats(tag, 3, "far/near IV, Backus filter")
    zone_stats(tag, 1, "mid/near IV, unfiltered")
    zone_stats(tag, 4, "mid/near IV, Backus filter")

print("\n== per-event agreement between the processings (far/near IV unfiltered and filtered): correlation over events")
t = sorted(RES)
for idx, lab in ((0, "unfiltered"), (3, "Backus filter")):
    x = np.log([RES[k]["07"][idx] for k in t])
    y = np.log([RES[k]["10"][idx] for k in t])
    print(f"   {lab}: r = {np.corrcoef(x, y)[0, 1]:+.2f}; mean log difference 94p10-94p07 {np.mean(y - x):+.3f}")
