
# ---- fm_sl_mult.py: are the events below the strongest reflector water-layer multiples of it? ----
# Run:  cat sl_common.py sl_events_lib.py fm_sl_mult.py | ssh linux-pc '<geoscience python> -I -'
# Flatten every stack on a parent reflector and look for a copy of it one and two water-layer periods later
# (sea bed at 109 ms two-way time): opposite sign after one period, same sign after two. A multiple follows the
# parent in structure and in amplitude, block by block; an independent primary does not.
S = stacks()
PRE, POST = 30, 150          # samples before and after the horizon: -60 ms to +300 ms


def flat_long(arr, i_c, tau_s):
    lo = i_c - MAXLAG - PRE
    W = np.asarray(arr[:, :, lo:i_c + MAXLAG + POST + 1], np.float32)
    idx = (MAXLAG + PRE + tau_s)[..., None] + np.arange(-PRE, POST + 1)
    return np.take_along_axis(W, idx, axis=2).astype(np.float64)


for T0 in (1316, 1276):
    i_c = T0 // 2
    tau_s, keep, ccb, lo = track(S["ful"], i_c)
    print(f"\n==== parent event at {T0} ms (kept {100*keep.mean():.0f} % of the traces)")
    for k in ("nea", "mid", "far", "ful"):
        F = flat_long(S[k], i_c, tau_s)
        mw = F[keep].mean(0)
        par = mw[PRE - 10:PRE + 21]                      # parent wavelet, -20 to +40 ms
        pn = (par ** 2).sum()
        out = []
        for name, (a, b), sign in (("one period", (40, 75), -1), ("two periods", (95, 130), +1)):
            best = (-9, 0, 0)
            for L in range(a, b):
                seg = mw[PRE - 10 + L:PRE + 21 + L]
                cc = sign * (seg * par).sum() / np.sqrt(pn * (seg ** 2).sum())
                if cc > best[0]:
                    best = (cc, L, (seg * par).sum() / pn)
            cc, L, sc = best
            # block by block: amplitude of the copy against amplitude of the parent
            pb = block_sum((F[:, :, PRE - 10:PRE + 21] * par).sum(-1) * keep) / pn
            cb = block_sum((F[:, :, PRE - 10 + L:PRE + 21 + L] * par).sum(-1) * keep) / pn
            nk = block_sum(keep.astype(float))
            ok = nk >= 0.3 * BS[0] * BS[1]
            pbm, cbm = pb[ok] / nk[ok], cb[ok] / nk[ok]
            r = np.corrcoef(pbm, cbm)[0, 1]
            slope = (pbm * cbm).sum() / (pbm * pbm).sum()
            out.append(f"{name}: best lag {L*2} ms, wavelet match {cc:+.2f}, amplitude {sc:+.2f} of the parent, block correlation {r:+.2f}, block slope {slope:+.2f}")
        print(f"  {k}: " + " | ".join(out))
        if k in ("nea", "far"):
            print(f"     mean flattened trace, {k}, every 6 ms from -30 to +282 ms:")
            print("     " + " ".join(f"{v:+.2f}" for v in mw[PRE - 15:PRE + 142:3]))

    # control: the same search at lags where no multiple is expected (60-78 ms and 150-180 ms), near stack
    F = flat_long(S["nea"], i_c, tau_s)
    mw = F[keep].mean(0)
    par = mw[PRE - 10:PRE + 21]
    pn = (par ** 2).sum()
    for name, (a, b) in (("control lags 60-78 ms", (30, 40)), ("control lags 150-184 ms", (75, 93))):
        best = max(((abs((mw[PRE - 10 + L:PRE + 21 + L] * par).sum()) / np.sqrt(pn * (mw[PRE - 10 + L:PRE + 21 + L] ** 2).sum()), L, (mw[PRE - 10 + L:PRE + 21 + L] * par).sum() / pn) for L in range(a, b)))
        L = best[1]
        pb = block_sum((F[:, :, PRE - 10:PRE + 21] * par).sum(-1) * keep) / pn
        cb = block_sum((F[:, :, PRE - 10 + L:PRE + 21 + L] * par).sum(-1) * keep) / pn
        nk = block_sum(keep.astype(float))
        ok = nk >= 0.3 * BS[0] * BS[1]
        r = np.corrcoef(pb[ok] / nk[ok], cb[ok] / nk[ok])[0, 1]
        print(f"  near, {name}: best |match| {best[0]:.2f} at {L*2} ms, amplitude {best[2]:+.2f} of the parent, block correlation {r:+.2f}")
