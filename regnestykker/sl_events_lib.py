
# ---- sl_events_lib.py: event machinery shared by the fm_sl_events*.py scripts (concatenated after sl_common.py) ----
# Track a flat event on the full stack, flatten near/mid/far on it, align the stacks block by block and
# estimate amplitude ratios with the third stack as instrument:
#     far/near = sum(far*mid) / sum(near*mid),   mid/near = sum(mid*far) / sum(near*far)
# Noise that is independent between the stacks does not bias these ratios. Uncertainty: bootstrap over blocks.
BS = (16, 18)            # block size in traces: 200 m x 225 m
HALF, MAXLAG, H = 15, 20, 24
NIL, NXL, NS = 249, 468, 1001
rng = np.random.default_rng(7)
BID = (np.arange(NIL) // BS[0])[:, None] * 1000 + (np.arange(NXL) // BS[1])[None, :]
UB = np.unique(BID)
IDXB = np.searchsorted(UB, BID)
NBI, NBX = (NIL - 1) // BS[0] + 1, (NXL - 1) // BS[1] + 1


def median_filter2(a, k=5):
    r = k // 2
    p = np.pad(a, r, mode="edge")
    st = np.stack([p[i:i + a.shape[0], j:j + a.shape[1]] for i in range(k) for j in range(k)], 0)
    return np.median(st, 0)


def block_sum(x):
    flat = x.reshape(NIL * NXL, *x.shape[2:])
    out = np.zeros((len(UB),) + x.shape[2:], np.float64)
    np.add.at(out, IDXB.ravel(), flat)
    return out


def track(full, i_c):
    lo = i_c - HALF - MAXLAG - H
    hi = i_c + HALF + MAXLAG + H + 1
    Wf = np.asarray(full[:, :, lo:hi], np.float32)
    off = MAXLAG + H
    pilot = Wf[:, :, off:off + 2 * HALF + 1].mean((0, 1)).astype(np.float64)
    tau_s = keep = ccbest = None
    for it in range(2):
        pn = np.sqrt((pilot ** 2).sum())
        nl = 2 * MAXLAG + 1
        cc = np.empty((nl, NIL, NXL), np.float32)
        for j in range(nl):
            seg = Wf[:, :, off - MAXLAG + j: off - MAXLAG + j + 2 * HALF + 1]
            cc[j] = (seg * pilot).sum(-1) / (np.sqrt((seg.astype(np.float64) ** 2).sum(-1)) * pn + 1e-20)
        cb = block_sum(np.moveaxis(cc, 0, -1))
        g = median_filter2((cb.argmax(-1) - MAXLAG).reshape(NBI, NBX).astype(float), 3)
        guide = g[(np.arange(NIL) // BS[0])[:, None], (np.arange(NXL) // BS[1])[None, :]]
        lag = np.arange(-MAXLAG, MAXLAG + 1)[:, None, None]
        masked = np.where(np.abs(lag - guide[None]) <= 5, cc, -9)
        tau = masked.argmax(0) - MAXLAG
        ccbest = masked.max(0)
        tau_s = np.rint(median_filter2(tau.astype(float), 5)).astype(int)
        keep = (ccbest > 0.4) & (np.abs(tau - tau_s) <= 2)
        idx = (off + tau_s)[..., None] + np.arange(2 * HALF + 1)
        pilot = np.take_along_axis(Wf, idx, axis=2)[keep].mean(0).astype(np.float64)
    return tau_s, keep, ccbest, lo


def flatten(arr, i_c, tau_s, lo):
    W = np.asarray(arr[:, :, lo:lo + 2 * (HALF + MAXLAG + H) + 1], np.float32)
    c = HALF + MAXLAG + H
    idx = (c + tau_s)[..., None] + np.arange(-H, H + 1)
    return np.take_along_axis(W, idx, axis=2).astype(np.float64)


def frac_shift(x, s):
    n = x.shape[-1]
    F = np.fft.rfft(x, axis=-1)
    f = np.fft.rfftfreq(n)
    return np.fft.irfft(F * np.exp(-2j * np.pi * f * np.asarray(s)[..., None]), n, axis=-1)


def block_lag(a, m, keep, maxlag=6, cw=10):
    c = a.shape[-1] // 2
    lags = np.arange(-maxlag, maxlag + 1)
    cf = np.empty((len(UB), len(lags)))
    k = keep[..., None]
    for j, L in enumerate(lags):
        cf[:, j] = block_sum(((m[:, :, c - cw:c + cw + 1] * a[:, :, c - cw + L:c + cw + 1 + L]) * k).sum(-1))
    jj = np.clip(cf.argmax(1), 1, len(lags) - 2)
    r = np.arange(len(UB))
    y0, y1, y2 = cf[r, jj - 1], cf[r, jj], cf[r, jj + 1]
    den = y0 - 2 * y1 + y2
    d = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
    return lags[jj] + np.clip(d, -1, 1)


def boot(num, den, n=1000, sel=None):
    ok = np.isfinite(num) & np.isfinite(den)
    if sel is not None:
        ok &= sel
    num, den = num[ok], den[ok]
    r = num.sum() / den.sum()
    idx = rng.integers(0, len(num), (n, len(num)))
    rs = num[idx].sum(1) / den[idx].sum(1)
    return r, np.percentile(rs, 2.5), np.percentile(rs, 97.5)


def bandpass(x, lo, hi):
    """Zero-phase band-pass of the last axis with cosine edges (Hz)."""
    n = x.shape[-1]
    f = np.fft.rfftfreq(n, DT)
    w = np.clip((f - 0.7 * lo) / (0.3 * lo), 0, 1) * np.clip((1.3 * hi - f) / (0.3 * hi), 0, 1)
    w = 0.5 - 0.5 * np.cos(np.pi * w)
    return np.fft.irfft(np.fft.rfft(x, axis=-1) * w, n, axis=-1)


def measure(S3, full, T0, band=None, cw=8, horizon=None):
    """IV amplitude ratios for the event at T0 ms. S3: dict with near, mid, far arrays. Returns a dict."""
    i_c = T0 // 2
    tau_s, keep, ccb, lo = horizon if horizon is not None else track(full, i_c)
    seg = {k: flatten(S3[k], i_c, tau_s, lo) for k in ("nea", "mid", "far")}
    if band:
        seg = {k: bandpass(v, *band) for k, v in seg.items()}
    raw = {k: v.copy() for k, v in seg.items()}
    ln = block_lag(seg["nea"], seg["mid"], keep)
    lf = block_lag(seg["far"], seg["mid"], keep)
    seg["nea"] = frac_shift(seg["nea"], -ln[IDXB])
    seg["far"] = frac_shift(seg["far"], -lf[IDXB])
    c = H
    k3 = keep[..., None]
    win = {k: seg[k][:, :, c - cw:c + cw + 1] * k3 for k in seg}
    P = {}
    for a, b in (("nea", "mid"), ("far", "mid"), ("nea", "far"), ("nea", "nea"), ("mid", "mid"), ("far", "far")):
        P[(a, b)] = block_sum((win[a] * win[b]).sum(-1))
    nk = block_sum(keep.astype(float))
    good = nk >= 0.3 * BS[0] * BS[1]
    for key in P:
        P[key] = np.where(good, P[key], np.nan)
    out = dict(T0=T0, kept=float(keep.mean()), good=good, P=P, nk=nk, tau=tau_s, keep=keep, ln=ln, lf=lf,
               fn=boot(P[("far", "mid")], P[("nea", "mid")]), mn=boot(P[("far", "mid")], P[("nea", "far")]),
               fm=boot(P[("nea", "far")], P[("nea", "mid")]))
    # mean wavelets: aligned (mw) and not aligned between stacks (mw_raw)
    out["mw"] = {k: seg[k][keep].mean(0) for k in seg}
    out["mw_raw"] = {k: raw[k][keep].mean(0) for k in raw}
    # alignment-free block estimate: extremum of the block-mean wavelet near the centre
    bm = {k: block_sum(raw[k] * k3) / np.maximum(nk, 1)[:, None] for k in raw}
    pol = np.sign(out["mw_raw"]["mid"][c - 2:c + 3][np.argmax(np.abs(out["mw_raw"]["mid"][c - 2:c + 3]))])
    pk = {}
    for k, hw in (("nea", 5), ("mid", 5), ("far", 8)):
        pk[k] = np.where(good, (pol * bm[k][:, c - hw:c + hw + 1]).max(1), np.nan)
    out["pk"] = pk
    out["pol"] = pol
    out["pk_fn"] = boot(pk["far"], pk["nea"])
    out["pk_mn"] = boot(pk["mid"], pk["nea"])
    return out
