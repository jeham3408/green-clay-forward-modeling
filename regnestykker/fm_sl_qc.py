#!/usr/bin/env python3
"""Sleipner near/mid/far: why do the partial stacks correlate so poorly?

Cross-spectra between the stacks (coherence, phase, transfer function), common-signal share per stack
(three-cornered hat), lateral continuity, autocorrelation (multiples), per-trace time shifts, and a picture.
Reads the arrays made by fm_sl_read.py. Writes one PNG to OUT.
"""
import json
import os
import numpy as np

OUT = "/mnt/ssd/Data/force-hackathon-2026/sleipner/avo"
CUBE = "/mnt/ssd/Data/force-hackathon-2026/sleipner/cubes/sleipner94"
DT = 0.002
S = {k: np.load(os.path.join(OUT, k + ".npy"), mmap_mode="r") for k in ("nea", "mid", "far", "ful")}
NIL, NXL, NS = S["nea"].shape


def cross_spectra(win, step=2, nfft=256):
    i0, i1 = int(win[0] / 2), int(win[1] / 2)
    w = np.hanning(i1 - i0)
    F = {k: np.fft.rfft(np.asarray(v[::step, ::step, i0:i1], np.float64) * w, nfft, axis=-1).reshape(-1, nfft // 2 + 1) for k, v in S.items() if k != "ful"}
    f = np.fft.rfftfreq(nfft, DT)
    auto = {k: np.mean(np.abs(F[k]) ** 2, 0) for k in F}
    cross = {(a, b): np.mean(F[a] * np.conj(F[b]), 0) for a, b in (("nea", "mid"), ("nea", "far"), ("mid", "far"))}
    return f, auto, cross


def band(f, x, lo, hi):
    m = (f >= lo) & (f < hi)
    return x[m].sum()


print("== cube meta and well file (structure only)")
meta = json.load(open(os.path.join(CUBE, "meta.json")))
print({k: (v if not isinstance(v, (list, dict)) else type(v).__name__ + str(len(v))) for k, v in meta.items()})
wells = json.load(open(os.path.join(CUBE, "wells.json")))
wl = wells if isinstance(wells, list) else wells.get("wells", wells)
print(type(wells).__name__, (list(wells.keys())[:20] if isinstance(wells, dict) else len(wells)))
w0 = wl[0] if isinstance(wl, list) else list(wl.values())[0]
print({k: (v if not isinstance(v, (list, dict)) else type(v).__name__ + str(len(v))) for k, v in w0.items()})
for k, v in w0.items():
    if isinstance(v, dict):
        print("  ", k, {kk: (vv if not isinstance(vv, (list, dict)) else type(vv).__name__ + str(len(vv))) for kk, vv in v.items()})

for win in ((860, 1060), (1100, 1500), (1500, 1900)):
    f, auto, cross = cross_spectra(win)
    print(f"\n== window {win[0]}-{win[1]} ms: per frequency band")
    print("  band Hz | coh n-m  n-f  m-f | phase n-f deg  n-m deg | |H f/n| IV   |H m/n| IV | common share nea  mid  far | power share nea")
    tot = band(f, auto["nea"], 4, 90)
    for lo in range(5, 85, 5):
        hi = lo + 5
        P = {k: band(f, auto[k], lo, hi) for k in auto}
        C = {k: band(f, v, lo, hi) for k, v in cross.items()}
        coh = {k: abs(C[k]) / np.sqrt(P[k[0]] * P[k[1]]) for k in C}
        ph_nf = np.degrees(np.angle(C[("nea", "far")]))
        ph_nm = np.degrees(np.angle(C[("nea", "mid")]))
        Hfn = abs(C[("mid", "far")]) / abs(C[("nea", "mid")])      # far/near seen through mid
        Hmn = abs(C[("mid", "far")]) / abs(C[("nea", "far")])      # mid/near seen through far
        sn = abs(C[("nea", "mid")]) * abs(C[("nea", "far")]) / abs(C[("mid", "far")])
        sm = abs(C[("nea", "mid")]) * abs(C[("mid", "far")]) / abs(C[("nea", "far")])
        sf = abs(C[("nea", "far")]) * abs(C[("mid", "far")]) / abs(C[("nea", "mid")])
        print(f"  {lo:3d}-{hi:3d} | {coh[('nea','mid')]:6.2f} {coh[('nea','far')]:5.2f} {coh[('mid','far')]:5.2f} | {ph_nf:+10.0f} {ph_nm:+10.0f} | {Hfn:9.3f} {Hmn:12.3f} | "
              f"{sn/P['nea']:14.2f} {sm/P['mid']:5.2f} {sf/P['far']:5.2f} | {P['nea']/tot:8.3f}")

print("\n== lateral continuity: mean correlation between neighbouring traces (12.5 m apart, crossline direction)")
print(" window ms |  near   mid   far  full")
for win in ((860, 1060), (1100, 1500), (1500, 1900)):
    i0, i1 = int(win[0] / 2), int(win[1] / 2)
    row = []
    for k in ("nea", "mid", "far", "ful"):
        v = np.asarray(S[k][::2, :, i0:i1], np.float64)
        a, b = v[:, :-1], v[:, 1:]
        cc = (a * b).sum(-1) / np.sqrt((a * a).sum(-1) * (b * b).sum(-1) + 1e-30)
        row.append(f"{cc.mean():5.2f}")
    print(f" {win[0]:4d}-{win[1]:4d} | " + " ".join(row))

print("\n== autocorrelation 1100-1900 ms (mean over traces), value at lag ms")
lags = [0, 10, 20, 30, 50, 80, 100, 104, 108, 112, 116, 120, 150, 200, 210, 216, 220, 230, 300]
print(" lag ms  " + " ".join(f"{L:6d}" for L in lags))
for k in ("nea", "mid", "far", "ful"):
    v = np.asarray(S[k][::3, ::3, 550:950], np.float64).reshape(-1, 400)
    F = np.fft.rfft(v, 1024, axis=-1)
    ac = np.fft.irfft(np.mean(np.abs(F) ** 2, 0))
    ac /= ac[0]
    print(f" {k}     " + " ".join(f"{ac[L // 2]:6.2f}" for L in lags))

print("\n== per-trace best lag of far against near, and of mid against near (parabolic peak, search +-10 ms)")
for win in ((1100, 1500), (1500, 1900)):
    i0, i1 = int(win[0] / 2), int(win[1] / 2)
    n = np.asarray(S["nea"][::2, ::2, i0:i1], np.float64).reshape(-1, i1 - i0)
    for other in ("mid", "far"):
        o = np.asarray(S[other][::2, ::2, i0:i1], np.float64).reshape(-1, i1 - i0)
        L = np.arange(-5, 6)
        cc = np.empty((len(L), n.shape[0]))
        den = np.sqrt((n * n).sum(-1) * (o * o).sum(-1))
        N = n.shape[1]
        for j, l in enumerate(L):
            cc[j] = ((n[:, :N - l] * o[:, l:]).sum(-1) if l >= 0 else (n[:, -l:] * o[:, :N + l]).sum(-1)) / den
        k = cc.argmax(0)
        best = cc.max(0)
        kk = np.clip(k, 1, len(L) - 2)
        idx = np.arange(n.shape[0])
        y0, y1, y2 = cc[kk - 1, idx], cc[kk, idx], cc[kk + 1, idx]
        d = np.where(np.abs(y0 - 2 * y1 + y2) > 1e-12, 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2), 0)
        lag = (L[kk] + np.clip(d, -1, 1)) * 2
        q = np.percentile(lag, [10, 25, 50, 75, 90])
        print(f" {win[0]}-{win[1]} near vs {other}: cc at zero lag median {np.median(cc[5]):.2f}; best cc median {np.median(best):.2f} (p10 {np.percentile(best,10):.2f}, p90 {np.percentile(best,90):.2f}); "
              f"lag ms p10/25/50/75/90 = {q[0]:+.1f}/{q[1]:+.1f}/{q[2]:+.1f}/{q[3]:+.1f}/{q[4]:+.1f}")

# picture: one inline through the well area, four panels
try:
    from PIL import Image
    il = 57
    t0, t1 = 300, 1000
    cols = []
    for k in ("nea", "mid", "far", "ful"):
        v = np.asarray(S[k][il, :, t0:t1], np.float64).T
        s = np.percentile(np.abs(v[250:]), 99)
        x = np.clip(v / s, -1, 1)
        rgb = np.stack([np.where(x > 0, 1, 1 + x), 1 - np.abs(x), np.where(x < 0, 1, 1 - x)], -1)
        cols.append((rgb * 255).astype(np.uint8))
        cols.append(np.full((t1 - t0, 6, 3), 0, np.uint8))
    img = np.concatenate(cols[:-1], 1)
    Image.fromarray(img).save(os.path.join(OUT, "qc_inline_near_mid_far_full.png"))
    print("\npicture written:", img.shape, "inline index", il, "time", t0 * 2, "-", t1 * 2, "ms; panels near | mid | far | full, crossline increases to the right")
except Exception as e:
    print("picture failed:", e)
