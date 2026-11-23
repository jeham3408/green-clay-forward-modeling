"""Scratch: how much does the centroid frequency of a window scatter, with no attenuation at all? Synthetic and the four cubes (read-only)."""
import sys, json, os, numpy as np
R = sys.argv[1]
rng = np.random.default_rng(7)
def centroid(x, dt_s, fmax=110.0):
    n = x.shape[-1]; w = np.hanning(n); nfft = 1 << int(np.ceil(np.log2(n * 4)))
    a = np.abs(np.fft.rfft(x * w, nfft, axis=-1)); f = np.fft.rfftfreq(nfft, dt_s); m = (f >= 3) & (f <= fmax)
    return (a[..., m] * f[m]).sum(-1) / a[..., m].sum(-1)
def ricker(fp, dt_s):
    half = int(np.ceil(1.5 / fp / dt_s)); t = np.arange(-half, half + 1) * dt_s; a = (np.pi * fp * t) ** 2
    return (1 - 2 * a) * np.exp(-a)
print("=== A. synthetic: white random reflectivity, Ricker 30 Hz, no attenuation. Difference in centroid between two windows of the same trace")
dt = 0.002
for win_ms in (200, 300, 500):
    n = int(win_ms / 1000 / dt); N = 4000
    r = rng.normal(size=(N, 2 * n + 400)); wv = ricker(30.0, dt)
    s = np.array([np.convolve(q, wv, mode="same") for q in r])
    c1 = centroid(s[:, 100:100 + n], dt); c2 = centroid(s[:, 300 + n:300 + 2 * n], dt)
    d = c2 - c1
    print("  window %3d ms: centroid %.1f Hz, single-window std %.2f Hz, std of the difference %.2f Hz -> traces with INDEPENDENT geology needed to see 0.7 / 1.3 / 2.5 Hz at 2 sigma: %d / %d / %d"
          % (win_ms, c1.mean(), c1.std(), d.std(), int(np.ceil((2 * d.std() / 0.7) ** 2)), int(np.ceil((2 * d.std() / 1.3) ** 2)), int(np.ceil((2 * d.std() / 2.5) ** 2))))
print("\n=== B. the reflectivity of Nini-4 itself (impedance log on the 2 ms grid), Ricker 30 Hz synthetic, centroid in sliding 300 ms windows")
d = json.load(open(R + "/nini/cubes/nini97/wells.json")); w = d["wells"][0]; L = w["logs"]
ai = np.array([np.nan if v is None else v for v in L["ai"]], float); ok = np.isfinite(ai); k0, k1 = np.flatnonzero(ok)[[0, -1]]
a = np.interp(np.arange(k0, k1 + 1), np.flatnonzero(ok), ai[ok]); rc = np.zeros(a.size); rc[1:] = (a[1:] - a[:-1]) / (a[1:] + a[:-1])
syn = np.convolve(rc, ricker(30.0, L["dt"] / 1000.0), mode="same"); t = L["t0"] + L["dt"] * np.arange(k0, k1 + 1)
n = int(300 / L["dt"]); vals = []
for s0 in range(0, syn.size - n + 1, n // 2):
    c = float(centroid(syn[s0:s0 + n], L["dt"] / 1000.0)); vals.append(c); print("  %4.0f-%4.0f ms: %.1f Hz" % (t[s0], t[s0 + n - 1], c))
print("  range %.1f Hz, std %.1f Hz (same wavelet everywhere; the differences are geology only)" % (max(vals) - min(vals), np.std(vals)))
print("\n=== C. the cubes: centroid frequency of consecutive windows in the overburden mudstone, 3000 random traces and 25x25-trace block means")
for name, p, wins in (("nini97", R + "/nini/cubes/nini97", [(900, 1200), (1200, 1500), (1280, 1580), (1680, 1980)]),
                      ("sleipner94", R + "/sleipner/cubes/sleipner94", [(1200, 1500), (1500, 1800), (1800, 2000)]),
                      ("tne01", R + "/smeaheia/cubes/tne01", [(600, 900), (900, 1200), (1200, 1500)]),
                      ("gn1101", R + "/smeaheia/cubes/gn1101", [(500, 800), (800, 1100), (1100, 1400)])):
    meta = json.load(open(p + "/meta.json")); cube = np.load(p + "/cube_ixt.npy", mmap_mode="r"); dtm, t0 = meta["dt_ms"], meta["t0_ms"]
    nI, nX, nT = cube.shape
    ii = rng.integers(0, nI, 3000); jj = rng.integers(0, nX, 3000); order = np.lexsort((jj, ii)); ii, jj = ii[order], jj[order]
    tr = np.stack([np.asarray(cube[i, j, :], np.float32) for i, j in zip(ii, jj)]); live = np.any(tr != 0, axis=1) & (np.abs(tr).max(axis=1) > 0)
    tr = tr[live]
    cs = []
    for a_, b_ in wins:
        ka, kb = int(round((a_ - t0) / dtm)), int(round((b_ - t0) / dtm))
        seg = tr[:, ka:kb]; good = np.abs(seg).sum(1) > 0
        c = np.full(tr.shape[0], np.nan); c[good] = centroid(seg[good].astype(np.float64), dtm / 1000.0, fmax=min(110.0, 0.9 * 500.0 / dtm)); cs.append(c)
        print("  %-10s %4d-%4d ms: centroid mean %.1f Hz, std between traces %.1f Hz (n=%d)" % (name, a_, b_, np.nanmean(c), np.nanstd(c), np.isfinite(c).sum()))
    for q in range(len(wins) - 1):
        dd = cs[q + 1] - cs[q]; print("  %-10s difference %s minus %s: mean %+.1f Hz, std between traces %.1f Hz" % (name, wins[q + 1], wins[q], np.nanmean(dd), np.nanstd(dd)))
    # block means: does averaging neighbours remove the scatter?
    a_, b_ = wins[0]; a2, b2 = wins[1]; ka, kb, kc, kd = [int(round((v - t0) / dtm)) for v in (a_, b_, a2, b2)]
    bm = []
    for _ in range(60):
        i0 = int(rng.integers(0, nI - 25)); j0 = int(rng.integers(0, nX - 25))
        blk = np.asarray(cube[i0:i0 + 25:2, j0:j0 + 25:2, :], np.float64).reshape(-1, nT); blk = blk[np.any(blk != 0, axis=1)]
        if blk.shape[0] < 50: continue
        bm.append(float(np.mean(centroid(blk[:, kc:kd], dtm / 1000.0, fmax=min(110.0, 0.9 * 500.0 / dtm)) - centroid(blk[:, ka:kb], dtm / 1000.0, fmax=min(110.0, 0.9 * 500.0 / dtm)))))
    print("  %-10s the same difference averaged over blocks of 25x25 traces (%d blocks): std between blocks %.1f Hz -> averaging neighbours leaves this much" % (name, len(bm), np.std(bm)))
print("\n=== D. the VSP file with the downgoing wavefield (headers only)")
vp = sys.argv[2]
raw = open(vp, "rb").read(3600); b = raw[3200:]
u2 = lambda byte: int.from_bytes(b[byte - 3201:byte - 3199], "big")
dtu, ns, fmt = u2(3217), u2(3221), u2(3225); size = os.path.getsize(vp); bps = {1: 4, 2: 4, 3: 2, 5: 4, 8: 1}.get(fmt, 4); ntr = (size - 3600) // (240 + bps * ns)
print("  %s: sample interval %d us, %d samples (%.0f ms), format code %d, %d traces" % (os.path.basename(vp), dtu, ns, ns * dtu / 1000.0, fmt, ntr))
