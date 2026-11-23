#!/usr/bin/env python3
"""Sleipner 1994 survey, 2010 processing (94p10): check the bag, inspect the SEG-Y headers, and put the four
stacks on the grid of the 2007 processing (94p07) so that reflectors can be compared trace by trace.

Read-only on the raw files. Writes float32 arrays (249 x 468 x NS) to OUT.
Run:  ssh linux-pc '<geoscience python> -I -' < fm_sl10_read.py
"""
import hashlib
import os
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
BAG = R + "/raw/seismic/2010 processing"
RAW = BAG + "/data"
RAW7 = R + "/raw/seismic/2007 processing/data/94p07nea.sgy"
OUT = R + "/avo10"
NAMES = {"nea": "94p10nea", "mid": "94p10mid", "far": "94p10far", "ful": "94p10ful"}


def ibm2ieee(u):
    u = u.astype(np.uint32)
    sign = np.where(u >> 31, -1.0, 1.0)
    exp = ((u >> 24) & 0x7F).astype(np.int32) - 64
    mant = (u & 0x00FFFFFF).astype(np.float64) / 16777216.0
    return (sign * mant * np.power(16.0, exp)).astype(np.float32)


def binhdr(p):
    raw = open(p, "rb").read(3600)
    txt = raw[:3200]
    ebc = txt.count(b"\x40") > txt.count(b"\x20")
    t = txt.decode("cp037" if ebc else "ascii", "replace")
    b = raw[3200:]
    u2 = lambda byte: int.from_bytes(b[byte - 3201:byte - 3199], "big")
    return t, u2(3217), u2(3221), u2(3225)


def headers(p, ns, fmt_bytes=4):
    size = os.path.getsize(p)
    tl = 240 + fmt_bytes * ns
    ntr = (size - 3600) // tl
    assert 3600 + ntr * tl == size, (size, ntr, tl)
    dt = np.dtype({"names": ["h"], "formats": [("V240")], "offsets": [0], "itemsize": tl})
    mm = np.memmap(p, dtype=dt, mode="r", offset=3600, shape=(ntr,))
    H = np.frombuffer(np.asarray(mm["h"]).tobytes(), dtype=np.uint8).reshape(ntr, 240)
    i4 = lambda byte: H[:, byte - 1:byte + 3].copy().view(">i4").ravel().astype(np.int64)
    i2 = lambda byte: H[:, byte - 1:byte + 1].copy().view(">i2").ravel().astype(np.int64)
    return ntr, tl, i4, i2


# ---------- 1. bag checksums
print("== bag check (sha256 against manifest)")
man = {}
for line in open(os.path.join(BAG, "manifest-sha256.txt")):
    h, f = line.strip().split(None, 1)
    man[f] = h
for f, h in sorted(man.items()):
    p = os.path.join(BAG, f)
    d = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            d.update(chunk)
    print(f"   {f}: {'OK' if d.hexdigest() == h else 'MISMATCH'} ({os.path.getsize(p)} bytes)")
print("bag-info:", open(os.path.join(BAG, "bag-info.txt")).read().strip().replace("\n", " | "))

# ---------- 2. headers
t, dtus, ns, fmt = binhdr(os.path.join(RAW, NAMES["nea"]))
print(f"\n== 94p10nea binary header: dt {dtus} us, ns {ns}, format {fmt}")
lines = [t[i * 80:(i + 1) * 80].rstrip() for i in range(40)]
print("text header lines with content beyond the template:")
for l in lines:
    s = l[3:].strip() if len(l) > 3 else ""
    if s:
        print("   ", l[:80])
ntr, tl, i4, i2 = headers(os.path.join(RAW, NAMES["nea"]), ns)
print(f"traces {ntr}, trace length {tl} bytes")
for nm, byte in (("tr seq line 1", 1), ("ffid 9", 9), ("tr in field 13", 13), ("cdp 21", 21), ("offset 37", 37),
                 ("cdp x 181", 181), ("cdp y 185", 185), ("il 189", 189), ("xl 193", 193), ("b17", 17), ("b25", 25)):
    v = i4(byte)
    print(f"   {nm:14s}: min {v.min():>12d} max {v.max():>12d} unique {len(np.unique(v)):>7d} first {v[:3].tolist()}")
for nm, byte in (("scalco 71", 71), ("delay 109", 109), ("ns 115", 115), ("dt 117", 117)):
    v = i2(byte)
    print(f"   {nm:14s}: min {v.min():>6d} max {v.max():>6d} first {v[:3].tolist()}")
il10, xl10, x10, y10, sc10 = i4(189), i4(193), i4(181), i4(185), i2(71)

# 94p07 grid for reference
t7, dt7, ns7, fmt7 = binhdr(RAW7)
n7, tl7, j4, j2 = headers(RAW7, ns7)
il7, xl7, x7, y7, sc7 = j4(189), j4(193), j4(181), j4(185), j2(71)
print(f"\n94p07: traces {n7}, il {il7.min()}-{il7.max()}, xl {xl7.min()}-{xl7.max()}, scalco {sc7[0]}")


def scale(v, s):
    s = int(s)
    return v / -s if s < 0 else v * (s if s > 0 else 1)


X7, Y7 = scale(x7.astype(float), sc7[0]), scale(y7.astype(float), sc7[0])
X10, Y10 = scale(x10.astype(float), sc10[0]), scale(y10.astype(float), sc10[0])
A = np.c_[il7, xl7, np.ones(n7)]
cx, *_ = np.linalg.lstsq(A, X7, rcond=None)
cy, *_ = np.linalg.lstsq(A, Y7, rcond=None)
print(f"94p07 grid: x = {cx[0]:.3f} il + {cx[1]:.3f} xl + {cx[2]:.1f}; y = {cy[0]:.3f} il + {cy[1]:.3f} xl + {cy[2]:.1f}; "
      f"fit residual {np.sqrt(np.mean((A @ cx - X7) ** 2 + (A @ cy - Y7) ** 2)):.3f} m")
M = np.array([[cx[0], cx[1]], [cy[0], cy[1]]])
f = np.linalg.solve(M, np.vstack([X10 - cx[2], Y10 - cy[2]]))
fil, fxl = f[0], f[1]
ril, rxl = np.rint(fil).astype(int), np.rint(fxl).astype(int)
off = np.hypot(*(M @ np.vstack([fil - ril, fxl - rxl])))
print(f"94p10 traces mapped onto the 94p07 grid: distance to nearest node median {np.median(off):.2f} m, p99 {np.percentile(off, 99):.2f} m")
same_num = np.mean((ril == il10) & (rxl == xl10))
print(f"share of 94p10 traces whose own il/xl numbers equal the mapped 94p07 node: {same_num:.3f}")
inside = (ril >= 1716) & (ril <= 1964) & (rxl >= 956) & (rxl <= 1423) & (off < 3.0)
print(f"94p10 traces inside the 94p07 box: {inside.sum()} of {ntr}; 94p07 nodes covered: {len(np.unique(ril[inside] * 10000 + rxl[inside]))} of {n7}")

# ---------- 3. convert, on the 94p07 grid
os.makedirs(OUT, exist_ok=True)
gi = ril[inside] - 1716
gx = rxl[inside] - 956
src = np.nonzero(inside)[0]
cov = np.zeros((249, 468), bool)
cov[gi, gx] = True
np.save(os.path.join(OUT, "coverage.npy"), cov)
for k, nm in NAMES.items():
    p = os.path.join(RAW, nm)
    out = os.path.join(OUT, k + ".npy")
    if os.path.exists(out):
        print("exists:", out)
        continue
    dtp = np.dtype([("h", "V240"), ("d", ">u4", ns)])
    mm = np.memmap(p, dtype=dtp, mode="r", offset=3600, shape=(ntr,))
    arr = np.lib.format.open_memmap(out, mode="w+", dtype=np.float32, shape=(249, 468, ns))
    arr[:] = 0
    order = np.argsort(src)
    for a in range(0, len(order), 6000):
        sel = order[a:a + 6000]
        d = ibm2ieee(np.asarray(mm["d"][src[sel]]))
        arr[gi[sel], gx[sel]] = d
    arr.flush()
    del arr
    print("written:", out)

# ---------- 4. quick comparison with 94p07 (full and near stacks)
a7 = {k: np.load(R + f"/avo/{k}.npy", mmap_mode="r") for k in ("nea", "far", "ful")}
a10 = {k: np.load(os.path.join(OUT, k + ".npy"), mmap_mode="r") for k in ("nea", "far", "ful")}
n = min(a7["ful"].shape[2], a10["ful"].shape[2])
print(f"\n== 94p10 against 94p07 on common traces (every 3rd), samples {n}")
print(" window ms | stack | rms 94p10/94p07 | cc at zero lag | best cc | lag ms (94p10 later +)")
for win in ((900, 1100), (1100, 1380), (1380, 1560), (1560, 1900)):
    i0, i1 = win[0] // 2, win[1] // 2
    m = cov[::3, ::3]
    for k in ("ful", "nea", "far"):
        u = np.asarray(a7[k][::3, ::3, i0:i1], np.float64)[m]
        v = np.asarray(a10[k][::3, ::3, i0:i1], np.float64)[m]
        best = (-9, 0)
        for L in range(-10, 11):
            if L >= 0:
                c = (u[:, :u.shape[1] - L] * v[:, L:]).sum() / np.sqrt((u ** 2).sum() * (v ** 2).sum())
            else:
                c = (u[:, -L:] * v[:, :v.shape[1] + L]).sum() / np.sqrt((u ** 2).sum() * (v ** 2).sum())
            if c > best[0]:
                best = (c, L)
        c0 = (u * v).sum() / np.sqrt((u ** 2).sum() * (v ** 2).sum())
        print(f" {win[0]}-{win[1]} | {k} | {np.sqrt((v**2).mean()/(u**2).mean()):8.3f} | {c0:+.2f} | {best[0]:+.2f} | {best[1]*2:+d}")
