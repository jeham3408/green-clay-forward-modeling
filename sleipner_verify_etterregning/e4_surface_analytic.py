"""Check of the pieces with an exact solution: the water-only runs against the analytic 2D line-source Green's
function with a perfect pressure-release surface (Lloyd mirror).  FS water run: p = g(r1) - g(r2), r1 = source (6 m)
to receiver (8 m), r2 = image source (-6 m) to receiver.  G water run: 2 [g(r1) - g(r2)].
g = H(t - r/c) / (2 pi sqrt(t^2 - r^2/c^2)), convolved with a unit Ricker(35 Hz) centred on t = 0.
Tests: time zero (lag), amplitude against offset (a surface coefficient different from -1 shows up at far offsets,
where the direct wave and its ghost nearly cancel), and the 35 -> 25 Hz conversion."""
import numpy as np
from scipy.signal import butter, filtfilt

import e1_lib as L

C = 1480.0


def ricker(f0, t):
    a = (np.pi * f0 * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


def analytic(offs, terms, f0, tgrid):
    dtf = 1e-4
    tf = np.arange(-0.1, tgrid[-1] + 0.1, dtf)
    tw = np.arange(-800, 801) * dtf
    w = ricker(f0, tw)
    out = np.zeros((len(offs), len(tgrid)))
    for i, x in enumerate(offs):
        acc = np.zeros(len(tf))
        for sgn, dz in terms:
            tr = np.hypot(x, dz) / C
            F = lambda tt: np.where(tt > tr, np.arccosh(np.maximum(tt, tr) / tr), 0.0) / (2 * np.pi)
            acc += sgn * (F(tf + dtf / 2) - F(tf - dtf / 2)) / dtf
        out[i] = np.interp(tgrid, tf, np.convolve(acc, w, mode="same") * dtf)
    return out


def cmp(fd, an, t, x, pre=0.06, post=0.10):
    w = (t > x / C - pre) & (t < x / C + post)
    a, b = fd[w], an[w]
    sc = float(np.sum(a * b) / np.sum(b * b))
    ncc = float(np.sum(a * b) / np.sqrt(np.sum(a * a) * np.sum(b * b)))
    lags = np.arange(-8, 9)
    i0, i1 = np.where(w)[0][[0, -1]]
    cc = np.array([np.sum(fd[i0 + l:i1 + 1 + l] * an[i0:i1 + 1]) for l in lags])
    j = int(np.clip(np.argmax(cc * np.sign(sc)), 1, len(lags) - 2))
    y0, y1, y2 = cc[j - 1], cc[j], cc[j + 1]
    lag = lags[j] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    res = float(np.sqrt(np.sum((a - sc * b) ** 2) / np.sum((sc * b) ** 2)))
    return sc, ncc, float(lag), res


R = {}
offs_sel = [150, 300, 500, 750, 1000, 1500, 2000, 2500, 3000, 3150]
bb, aa = butter(2, 3.0 / 500.0, "high")
for cfg, terms, fac in (("FS", [(1, 2.0), (-1, 14.0)], 1.0), ("G", [(1, 2.0), (-1, 14.0)], 2.0)):
    d, t, off = L.load_gather(cfg, "water")
    ix = [int(np.argmin(np.abs(off - x))) for x in offs_sel]
    an35 = fac * analytic(off[ix], terms, 35.0, t)
    an35f = filtfilt(bb, aa, an35, axis=-1)
    rows = []
    for j, i in enumerate(ix):
        sc, ncc, lag, res = cmp(d[i], an35f[j], t, off[i])
        sc0, ncc0, lag0, res0 = cmp(d[i], an35[j], t, off[i])
        rows.append(dict(offset=float(off[i]), scale=sc, ncc=ncc, lag_ms=lag, resid_rel=res, scale_no_lowcut=sc0, ncc_no_lowcut=ncc0))
    s0 = rows[0]["scale"]
    for r in rows:
        r["scale_rel_to_150m"] = r["scale"] / s0
    R[cfg + "_35Hz"] = rows
    print("\n%s water run against the analytic Lloyd-mirror solution (35 Hz, as modelled):" % cfg)
    for r in rows:
        print("  x=%5.0f  scale/scale(150)=%.4f  ncc=%.5f  lag=%+.3f ms  residual=%.4f   (without the 3 Hz low-cut on the analytic: ncc=%.5f)" % (
            r["offset"], r["scale_rel_to_150m"], r["ncc"], r["lag_ms"], r["resid_rel"], r["ncc_no_lowcut"]))
    # sensitivity: what a surface coefficient of -0.99 would do (1 % of the image-source term left over)
    if cfg == "FS":
        sens = []
        for j, i in enumerate(ix):
            leak = 0.01 * analytic(off[[i]], [(1, 14.0)], 35.0, t)[0]
            w = (t > off[i] / C - 0.06) & (t < off[i] / C + 0.10)
            sens.append(dict(offset=float(off[i]), rms_leak_over_rms_signal=float(np.sqrt(np.mean(leak[w] ** 2) / np.mean(an35[j][w] ** 2)))))
        R["sensitivity_surface_coefficient_minus_0.99"] = sens
        print("  a surface coefficient of -0.99 instead of -1 would leave (rms, relative to the signal):", [(int(s["offset"]), round(s["rms_leak_over_rms_signal"], 3)) for s in sens])
    # conversion 35 -> 25 Hz: converted FD against the analytic solution for a 25 Hz Ricker, both with sqrt(i omega)
    conv = L.to_wavelet(d[ix], 35.0, 25.0, point=True)
    an25 = L.to_wavelet(filtfilt(bb, aa, fac * analytic(off[ix], terms, 25.0, t), axis=-1), 25.0, 25.0, point=True)
    ref35 = L.to_wavelet(d[ix], 35.0, 35.0, point=True)
    an35p = L.to_wavelet(an35f, 35.0, 35.0, point=True)
    rows = []
    for j, i in enumerate(ix):
        sc, ncc, lag, res = cmp(conv[j], an25[j], t, off[i])
        sc35, ncc35, lag35, res35 = cmp(ref35[j], an35p[j], t, off[i])
        rows.append(dict(offset=float(off[i]), scale_25=sc, scale_35=sc35, scale_25_over_35=sc / sc35, ncc_25=ncc, lag_ms_25=lag, resid_rel_25=res, resid_rel_35=res35))
    R[cfg + "_conversion_25Hz"] = rows
    print("  35 -> 25 Hz conversion (own filter, unit-peak Rickers): scale(25)/scale(35), ncc, lag, residual 25 / residual 35")
    for r in rows:
        print("  x=%5.0f  %.4f  ncc=%.5f  lag=%+.3f ms  resid25=%.4f resid35=%.4f" % (r["offset"], r["scale_25_over_35"], r["ncc_25"], r["lag_ms_25"], r["resid_rel_25"], r["resid_rel_35"]))
# their filter uses (f_from/f_to)^2, mine ^3: constant factor
R["their_conversion_constant_factor_vs_unit_peak"] = 25.0 / 35.0
L.jdump(R, "e4_results.json")
print("saved")
