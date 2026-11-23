#!/usr/bin/env python3
"""Sleipner 94p07 near/mid/far/full stacks: SEG-Y (IBM float) -> float32 arrays, plus basic controls.

Reads the raw open files only. Writes derived arrays to OUT, outside the app folder.
Run:  ssh linux-pc '<geoscience python> -I -' < fm_sl_read.py
"""
import os
import numpy as np

RAW = "/mnt/ssd/Data/force-hackathon-2026/sleipner/raw/seismic/2007 processing/data"
OUT = "/mnt/ssd/Data/force-hackathon-2026/sleipner/avo"
NS, NIL, NXL, DT = 1001, 249, 468, 0.002
NAMES = {"nea": "94p07nea.sgy", "mid": "94p07mid.sgy", "far": "94p07far.sgy", "ful": "94p07ful.sgy"}


def ibm2ieee(u):
    u = u.astype(np.uint32)
    sign = np.where(u >> 31, -1.0, 1.0)
    exp = ((u >> 24) & 0x7F).astype(np.int32) - 64
    mant = (u & 0x00FFFFFF).astype(np.float64) / 16777216.0
    return (sign * mant * np.power(16.0, exp)).astype(np.float32)


def read(key):
    out = os.path.join(OUT, key + ".npy")
    if os.path.exists(out):
        return np.load(out, mmap_mode="r")
    p = os.path.join(RAW, NAMES[key])
    ntr = (os.path.getsize(p) - 3600) // (240 + 4 * NS)
    assert ntr == NIL * NXL, ntr
    dt = np.dtype([("h", ">i4", 60), ("d", ">u4", NS)])
    mm = np.memmap(p, dtype=dt, mode="r", offset=3600, shape=(ntr,))
    h = np.asarray(mm["h"][[0, NXL - 1, NXL, ntr - 1]])
    assert (h[:, 47].tolist(), h[:, 48].tolist()) == ([1716, 1716, 1717, 1964], [956, 1423, 956, 1423]), h[:, 47:49]
    arr = np.empty((ntr, NS), np.float32)
    for a in range(0, ntr, 8000):
        arr[a:a + 8000] = ibm2ieee(np.asarray(mm["d"][a:a + 8000]))
    arr = arr.reshape(NIL, NXL, NS)
    np.save(out, arr)
    return arr


def spectrum(x, nfft=256):
    w = np.hanning(x.shape[-1])
    X = np.abs(np.fft.rfft(x * w, nfft, axis=-1)) ** 2
    P = X.reshape(-1, X.shape[-1]).mean(0)
    f = np.fft.rfftfreq(nfft, DT)
    A = np.sqrt(P)
    cen = (f * A).sum() / A.sum()
    cum = np.cumsum(P) / P.sum()
    med = f[np.searchsorted(cum, 0.5)]
    return cen, f[np.argmax(A)], med


def xcorr_lag(a, b, maxlag=6):
    """Mean normalised cross-correlation of b against a over traces; positive lag = b later than a."""
    n = a.shape[-1]
    lags = np.arange(-maxlag, maxlag + 1)
    cc = []
    den = np.sqrt((a ** 2).sum(-1) * (b ** 2).sum(-1)) + 1e-30
    for L in lags:
        if L >= 0:
            num = (a[..., :n - L] * b[..., L:]).sum(-1)
        else:
            num = (a[..., -L:] * b[..., :n + L]).sum(-1)
        cc.append(np.mean(num / den))
    cc = np.array(cc)
    k = int(np.argmax(cc))
    if 0 < k < len(cc) - 1:
        y0, y1, y2 = cc[k - 1], cc[k], cc[k + 1]
        d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    else:
        d = 0.0
    return cc[maxlag], cc[k], (lags[k] + d) * DT * 1000


def main():
    os.makedirs(OUT, exist_ok=True)
    S = {k: read(k) for k in NAMES}
    print("arrays:", {k: v.shape for k, v in S.items()})
    t = np.arange(NS) * DT * 1000
    sub = {k: np.asarray(v[::3, ::3, :], dtype=np.float64) for k, v in S.items()}
    ntr = sub["nea"].shape[0] * sub["nea"].shape[1]
    print("subset traces:", ntr)

    print("\n== dead traces and first live sample")
    for k, v in sub.items():
        live = np.abs(v) > 0
        dead = (~live.any(-1)).sum()
        first = np.where(live.any(-1), live.argmax(-1), NS)
        fl = first[first < NS]
        print(f"{k}: dead traces {dead}; first live sample time ms: min {fl.min()*2} median {np.median(fl)*2:.0f} p95 {np.percentile(fl,95)*2:.0f} max {fl.max()*2}; "
              f"fraction of exact zeros inside 300-1900 ms {np.mean(v[..., 150:950] == 0):.5f}; abs max {np.abs(v).max():.4g}")

    print("\n== RMS per 100 ms window (all subset traces) and ratios")
    print(" t0-t1 ms |   near     mid     far    full | mid/nea far/nea ful/nea | live% nea mid far")
    for a in range(0, 2000, 100):
        i0, i1 = a // 2, (a + 100) // 2
        r = {k: np.sqrt(np.mean(v[..., i0:i1] ** 2)) for k, v in sub.items()}
        lv = {k: 100 * np.mean(np.abs(v[..., i0:i1]).sum(-1) > 0) for k, v in sub.items()}
        n = r["nea"] + 1e-30
        print(f"{a:5d}-{a+100:4d} | {r['nea']:7.3g} {r['mid']:7.3g} {r['far']:7.3g} {r['ful']:7.3g} | {r['mid']/n:6.3f} {r['far']/n:7.3f} {r['ful']/n:7.3f} | {lv['nea']:5.0f} {lv['mid']:4.0f} {lv['far']:4.0f}")

    print("\n== spectra per 200 ms window: centroid / peak / median frequency (Hz)")
    print(" t0-t1 ms |  near c/p/m      mid c/p/m      far c/p/m      full c/p/m   | cen far/nea  mid/nea")
    for a in range(100, 1900, 200):
        i0, i1 = a // 2, (a + 200) // 2
        sp = {k: spectrum(v[..., i0:i1]) for k, v in sub.items()}
        print(f"{a:5d}-{a+200:4d} | " + "  ".join(f"{sp[k][0]:5.1f}/{sp[k][1]:4.1f}/{sp[k][2]:4.1f}" for k in ("nea", "mid", "far", "ful"))
              + f" | {sp['far'][0]/sp['nea'][0]:6.3f} {sp['mid'][0]/sp['nea'][0]:8.3f}")

    print("\n== correlation between stacks per 200 ms window: cc at zero lag, best cc, best lag (ms; positive = second stack later)")
    print(" t0-t1 ms |   near-mid            near-far            mid-far             near-full")
    for a in range(100, 1900, 200):
        i0, i1 = a // 2, (a + 200) // 2
        row = []
        for p, q in (("nea", "mid"), ("nea", "far"), ("mid", "far"), ("nea", "ful")):
            c0, cb, lag = xcorr_lag(sub[p][..., i0:i1], sub[q][..., i0:i1])
            row.append(f"{c0:5.2f} {cb:5.2f} {lag:+5.2f}")
        print(f"{a:5d}-{a+200:4d} | " + "   ".join(row))

    print("\n== full stack as a combination of the three: ful ~ a*nea + b*mid + c*far (least squares per window)")
    print(" t0-t1 ms |     a       b       c    sum   | R2")
    for a in range(100, 1900, 200):
        i0, i1 = a // 2, (a + 200) // 2
        X = np.stack([sub[k][..., i0:i1].ravel() for k in ("nea", "mid", "far")], 1)
        y = sub["ful"][..., i0:i1].ravel()
        co, *_ = np.linalg.lstsq(X, y, rcond=None)
        r2 = 1 - np.sum((y - X @ co) ** 2) / np.sum(y ** 2)
        print(f"{a:5d}-{a+200:4d} | {co[0]:7.3f} {co[1]:7.3f} {co[2]:7.3f} {co.sum():6.3f} | {r2:.4f}")

    print("\n== amplitude ratio estimators per 200 ms window")
    print(" t0-t1 ms | rms far/nea | slope cov(f,n)/var(n) | IV cov(f,m)/cov(n,m) | rms mid/nea | IV cov(m,f)/cov(n,f)")
    for a in range(100, 1900, 200):
        i0, i1 = a // 2, (a + 200) // 2
        n, m, f = (sub[k][..., i0:i1].ravel() for k in ("nea", "mid", "far"))
        print(f"{a:5d}-{a+200:4d} | {np.sqrt(np.mean(f**2)/np.mean(n**2)):11.3f} | {np.mean(f*n)/np.mean(n*n):21.3f} | {np.mean(f*m)/np.mean(n*m):20.3f} | "
              f"{np.sqrt(np.mean(m**2)/np.mean(n**2)):11.3f} | {np.mean(m*f)/np.mean(n*f):20.3f}")


if __name__ == "__main__":
    main()
