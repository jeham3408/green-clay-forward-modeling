
# ---- fm_sl_events.py: amplitude against offset, reflector by reflector, over the whole cube ----
# Run:  cat sl_common.py fm_sl_events.py | ssh linux-pc '<geoscience python> -I -'
#
# For each flat event: track it on the full stack, flatten near/mid/far on it, align the stacks block by block,
# and estimate the amplitude ratios with the third stack as instrument:
#     far/near = sum(far*mid) / sum(near*mid),   mid/near = sum(mid*far) / sum(near*far),   far/mid = sum(far*near) / sum(mid*near)
# Noise that is independent between the stacks does not bias these ratios. Uncertainty: bootstrap over blocks.
import sys

S = stacks()
NIL, NXL, NS = S["nea"].shape
BS = (16, 18)            # block size in traces: 200 m x 225 m
HALF, MAXLAG, H = 15, 20, 24
EVENTS = [(1104, -1), (1140, +1), (1228, -1), (1276, +1), (1316, -1), (1332, +1), (1386, -1), (1456, +1), (1536, -1), (1746, +1), (1776, -1)]
rng = np.random.default_rng(7)


def median_filter2(a, k=5):
    r = k // 2
    p = np.pad(a, r, mode="edge")
    st = np.stack([p[i:i + a.shape[0], j:j + a.shape[1]] for i in range(k) for j in range(k)], 0)
    return np.median(st, 0)


def blocks_of(shape):
    bi = np.arange(shape[0]) // BS[0]
    bx = np.arange(shape[1]) // BS[1]
    return bi[:, None] * 1000 + bx[None, :]


BID = blocks_of((NIL, NXL))
UB = np.unique(BID)


def block_sum(x):
    """Sum a (NIL, NXL, ...) array over blocks -> (nblocks, ...)."""
    flat = x.reshape(NIL * NXL, *x.shape[2:])
    idx = np.searchsorted(UB, BID.ravel())
    out = np.zeros((len(UB),) + x.shape[2:], np.float64)
    np.add.at(out, idx, flat)
    return out


def track(i_c):
    lo = i_c - HALF - MAXLAG - H
    hi = i_c + HALF + MAXLAG + H + 1
    Wf = np.asarray(S["ful"][:, :, lo:hi], np.float32)
    off = MAXLAG + H                     # index of (i_c - HALF) inside Wf
    pilot = Wf[:, :, off:off + 2 * HALF + 1].mean((0, 1)).astype(np.float64)
    tau_s, keep, ccbest = None, None, None
    for it in range(2):
        pn = np.sqrt((pilot ** 2).sum())
        nl = 2 * MAXLAG + 1
        cc = np.empty((nl, NIL, NXL), np.float32)
        for j in range(nl):
            seg = Wf[:, :, off - MAXLAG + j: off - MAXLAG + j + 2 * HALF + 1]
            cc[j] = (seg * pilot).sum(-1) / (np.sqrt((seg.astype(np.float64) ** 2).sum(-1)) * pn + 1e-20)
        # coarse guide from block-averaged correlation
        cb = block_sum(np.moveaxis(cc, 0, -1))              # (nblocks, nl)
        gl = cb.argmax(-1) - MAXLAG
        nbi, nbx = len(np.unique(BID[:, 0] // 1000)), len(np.unique(BID[0, :] % 1000))
        g = median_filter2(gl.reshape(nbi, nbx).astype(float), 3)
        guide = g[(np.arange(NIL) // BS[0])[:, None], (np.arange(NXL) // BS[1])[None, :]]
        # fine: best lag within guide +- 5 samples
        lag = np.arange(-MAXLAG, MAXLAG + 1)[:, None, None]
        masked = np.where(np.abs(lag - guide[None]) <= 5, cc, -9)
        tau = masked.argmax(0) - MAXLAG
        ccbest = masked.max(0)
        tau_s = np.rint(median_filter2(tau.astype(float), 5)).astype(int)
        keep = (ccbest > 0.4) & (np.abs(tau - tau_s) <= 2)
        # new pilot from the flattened, kept traces
        idx = (off + tau_s)[..., None] + np.arange(2 * HALF + 1)
        fl = np.take_along_axis(Wf, idx, axis=2)
        pilot = fl[keep].mean(0).astype(np.float64)
    return tau_s, keep, ccbest, lo


def flatten(key, i_c, tau_s, lo):
    W = np.asarray(S[key][:, :, lo:lo + 2 * (HALF + MAXLAG + H) + 1], np.float32)
    c = HALF + MAXLAG + H                # index of i_c inside W
    idx = (c + tau_s)[..., None] + np.arange(-H, H + 1)
    return np.take_along_axis(W, idx, axis=2).astype(np.float64)


def frac_shift(x, s):
    """Shift the last axis by s samples (positive = later), s broadcastable to x[..., 0]."""
    n = x.shape[-1]
    F = np.fft.rfft(x, axis=-1)
    f = np.fft.rfftfreq(n)
    return np.fft.irfft(F * np.exp(-2j * np.pi * f * np.asarray(s)[..., None]), n, axis=-1)


def block_lag(a, m, keep, maxlag=6, cw=10):
    """Lag (samples) of a relative to m per block: positive = a later than m. Block-summed cross-correlation."""
    c = a.shape[-1] // 2
    lags = np.arange(-maxlag, maxlag + 1)
    cf = np.empty((len(UB), len(lags)))
    k = keep[..., None]
    for j, L in enumerate(lags):
        prod = (m[:, :, c - cw:c + cw + 1] * a[:, :, c - cw + L:c + cw + 1 + L]) * k
        cf[:, j] = block_sum(prod.sum(-1))
    j = cf.argmax(1)
    jj = np.clip(j, 1, len(lags) - 2)
    r = np.arange(len(UB))
    y0, y1, y2 = cf[r, jj - 1], cf[r, jj], cf[r, jj + 1]
    den = y0 - 2 * y1 + y2
    d = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
    return lags[jj] + np.clip(d, -1, 1)


def boot(num, den, n=1000):
    """Ratio of block sums with a block bootstrap interval."""
    ok = np.isfinite(num) & np.isfinite(den)
    num, den = num[ok], den[ok]
    r = num.sum() / den.sum()
    idx = rng.integers(0, len(num), (n, len(num)))
    rs = num[idx].sum(1) / den[idx].sum(1)
    return r, np.percentile(rs, 2.5), np.percentile(rs, 97.5)


print("event  pol | kept % | time spread ms p5/p50/p95 | lag near-mid, far-mid ms (median, iqr) | far/near (95 %)        mid/near (95 %)        far/mid (95 %)       | mean-wavelet far/near mid/near | common share near mid far | block scatter far/near, mid/near (robust sd of log)")
RES = {}
for T0, pol in EVENTS:
    i_c = T0 // 2
    tau_s, keep, ccb, lo = track(i_c)
    seg = {k: flatten(k, i_c, tau_s, lo) for k in ("nea", "mid", "far")}
    ln = block_lag(seg["nea"], seg["mid"], keep)
    lf = block_lag(seg["far"], seg["mid"], keep)
    idxb = np.searchsorted(UB, BID)
    seg["nea"] = frac_shift(seg["nea"], -ln[idxb])
    seg["far"] = frac_shift(seg["far"], -lf[idxb])
    c, cw = H, 8
    k3 = keep[..., None]
    win = {k: seg[k][:, :, c - cw:c + cw + 1] * k3 for k in seg}
    P = {}
    for a, b in (("nea", "mid"), ("far", "mid"), ("nea", "far"), ("nea", "nea"), ("mid", "mid"), ("far", "far")):
        P[(a, b)] = block_sum((win[a] * win[b]).sum(-1))
    nk = block_sum(keep.astype(float))
    good = nk >= 0.3 * BS[0] * BS[1]
    for key in P:
        P[key] = np.where(good, P[key], np.nan)
    fn = boot(P[("far", "mid")], P[("nea", "mid")])
    mn = boot(P[("far", "mid")], P[("nea", "far")])
    fm = boot(P[("nea", "far")], P[("nea", "mid")])
    # mean wavelet after flattening and alignment
    mw = {k: seg[k][keep].mean(0) for k in seg}
    amp = {k: mw[k][c - 2:c + 3][np.argmax(np.abs(mw[k][c - 2:c + 3]))] for k in mw}
    tot = {k: np.nansum(P[(k, k)]) for k in ("nea", "mid", "far")}
    snm, sfm, snf = np.nansum(P[("nea", "mid")]), np.nansum(P[("far", "mid")]), np.nansum(P[("nea", "far")])
    share = (snm * snf / sfm / tot["nea"], snm * sfm / snf / tot["mid"], snf * sfm / snm / tot["far"])
    # block-level scatter of the ratios
    with np.errstate(divide="ignore", invalid="ignore"):
        bfn = P[("far", "mid")] / P[("nea", "mid")]
        bmn = P[("far", "mid")] / P[("nea", "far")]
    def rsd(x):
        x = x[np.isfinite(x) & (x > 0)]
        q = np.percentile(np.log(x), [25, 75])
        return (q[1] - q[0]) / 1.349, len(x)
    s1, n1 = rsd(bfn)
    s2, n2 = rsd(bmn)
    ts = (tau_s[keep] * 2.0)
    RES[T0] = dict(fn=fn, mn=mn, fm=fm, amp=amp, share=share, kept=keep.mean(), sfn=s1, smn=s2, nb=int(good.sum()))
    print(f"{T0:5d}  {pol:+d} | {100*keep.mean():5.1f} | {np.percentile(ts,5):+5.0f}/{np.percentile(ts,50):+4.0f}/{np.percentile(ts,95):+4.0f}       | "
          f"{np.median(ln[good])*2:+5.1f} ({(np.percentile(ln[good],75)-np.percentile(ln[good],25))*2:4.1f}), {np.median(lf[good])*2:+5.1f} ({(np.percentile(lf[good],75)-np.percentile(lf[good],25))*2:4.1f}) | "
          f"{fn[0]:5.2f} ({fn[1]:.2f}-{fn[2]:.2f})   {mn[0]:5.2f} ({mn[1]:.2f}-{mn[2]:.2f})   {fm[0]:5.2f} ({fm[1]:.2f}-{fm[2]:.2f}) | "
          f"{amp['far']/amp['nea']:6.2f} {amp['mid']/amp['nea']:6.2f}  [{amp['nea']:+.3f} {amp['mid']:+.3f} {amp['far']:+.3f}] | "
          f"{share[0]:.2f} {share[1]:.2f} {share[2]:.2f} | {s1:.2f} (n={n1}), {s2:.2f} (n={n2})")
    sys.stdout.flush()
    np.save(os.path.join(AVO, f"event_{T0}_blocks.npy"), np.stack([P[("far", "mid")], P[("nea", "mid")], P[("nea", "far")], nk]))

print("\n== the same ratios relative to the median of all events (removes the unknown scaling between the stacks)")
mfn = np.median([RES[t]["fn"][0] for t in RES])
mmn = np.median([RES[t]["mn"][0] for t in RES])
print(f"   median over events: far/near {mfn:.2f}, mid/near {mmn:.2f}")
for T0, pol in EVENTS:
    r = RES[T0]
    print(f"   {T0}: far/near {r['fn'][0]/mfn:5.2f} ({r['fn'][1]/mfn:.2f}-{r['fn'][2]/mfn:.2f})   mid/near {r['mn'][0]/mmn:5.2f} ({r['mn'][1]/mmn:.2f}-{r['mn'][2]/mmn:.2f})")
