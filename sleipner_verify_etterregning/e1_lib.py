"""Etterregning (independent recomputation) of the Sleipner synthetic experiment: own processing library.
Does NOT import slp.py or s3_stack_measure.py. Only numpy / scipy.
"""
import json
import os

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import hilbert

SRC = os.environ.get("SRC", "/workspace/sleipner_exp")
OUT = os.environ.get("OUT", "/workspace/sleipner_verify_etterregning")
DT = 0.001
T0 = np.arange(0.0, 2.2005, DT)
RANGES = {"near": (150.0, 760.0), "mid": (760.0, 1760.0), "far": (1760.0, 3150.0)}
ZONE = (1.408, 1.492)
EXCL = (1.380, 1.560)
TLIM = (1.100, 1.900)
VW, HW = 1480.0, 81.0


def load_gather(cfg, case, src=None):
    z = np.load(os.path.join(src or SRC, "gather_%s_%s.npz" % (cfg, case)))
    return z["d"].astype(float), z["t"].astype(float), z["offsets"].astype(float)


def load_model(src=None):
    M = np.load(os.path.join(src or SRC, "model_1d.npz"))
    return {k: M[k].copy() for k in ("z", "ze", "vp", "vs", "rho", "twt", "soft")}


def to_wavelet(d, f_from=35.0, f_to=25.0, point=True):
    """Ricker(f_from) -> Ricker(f_to), both unit peak: ratio of the analytic spectra 2/sqrt(pi) f^2/f0^3 exp(-f^2/f0^2).
    point=True: multiply by sqrt(i omega) (undo the far-field half-integration of a 2D line source)."""
    n = d.shape[-1]
    nfft = 4 * 2 ** int(np.ceil(np.log2(n)))
    f = np.fft.rfftfreq(nfft, DT)
    H = (f_from / f_to) ** 3 * np.exp(-f ** 2 * (1.0 / f_to ** 2 - 1.0 / f_from ** 2)) if abs(f_from - f_to) > 1e-9 else np.ones(len(f))
    H = H.astype(complex)
    if point:
        H = H * np.sqrt(2j * np.pi * f)
    return np.fft.irfft(np.fft.rfft(d, nfft, axis=-1) * H, nfft, axis=-1)[..., :n]


def box_slow(v, n):
    if n <= 1:
        return v.copy()
    s = np.convolve(np.pad(1.0 / v, n // 2, mode="edge"), np.ones(n) / n, mode="valid")
    return 1.0 / s


def raytab(v, offsets, n_p=3000, n_water=0, pmax=1.0 / 1700.0):
    """Two-way reflection time / ray parameter / dX/dp at the base of every 1 m cell, surface to surface, for the
    given offsets. n_water: extra round trips in the water layer (peg-leg of that order). NaN where no ray."""
    n = len(v)
    p = np.linspace(0.0, pmax, n_p)
    pv2 = (p[None, :] * v[:, None]) ** 2
    with np.errstate(invalid="ignore"):
        c = np.sqrt(np.where(pv2 < 1.0, 1.0 - pv2, np.nan))
    X = 2.0 * np.cumsum(p[None, :] * v[:, None] / c, axis=0)
    T = 2.0 * np.cumsum(1.0 / (v[:, None] * c), axis=0)
    if n_water:
        cw = np.sqrt(1.0 - (p * VW) ** 2)
        X = X + n_water * 2.0 * HW * p * VW / cw
        T = T + n_water * 2.0 * HW / (VW * cw)
    nx = len(offsets)
    tx = np.full((n, nx), np.nan)
    px = np.full((n, nx), np.nan)
    dxdp = np.full((n, nx), np.nan)
    for k in range(n):
        g = np.isfinite(X[k])
        if g.sum() < 5:
            continue
        Xk = X[k, g]
        tx[k] = np.interp(offsets, Xk, T[k, g], right=np.nan)
        px[k] = np.interp(offsets, Xk, p[g], right=np.nan)
        dxdp[k] = np.interp(offsets, Xk, np.gradient(Xk, p[g]), right=np.nan)
    return tx, px, dxdp


class Kin:
    def __init__(self, m, offsets, smooth=11):
        self.off = np.asarray(offsets, float)
        v = m["vp"].copy()
        v[m["z"] < HW] = VW
        self.v = box_slow(v, smooth)
        self.v[m["z"] < HW] = VW
        self.t0 = 2.0 * np.cumsum(1.0 / self.v)            # two-way vertical time to the base of each cell
        self.tx, self.px, self.dxdp = raytab(self.v, self.off)
        self.vbig = box_slow(v, 201)                       # for a smooth angle mute
        with np.errstate(invalid="ignore"):
            self.ang = np.degrees(np.arcsin(np.clip(self.px * self.vbig[:, None], 0, 1)))


def nmo(d, t, kin, mute=("angle", 48.0)):
    """Exact-time NMO (own ray table), cubic interpolation along time. Returns corrected gather and live mask."""
    nx = d.shape[0]
    out = np.zeros((nx, len(T0)))
    live = np.zeros((nx, len(T0)), bool)
    for ix in range(nx):
        tk = kin.tx[:, ix]
        ok = np.isfinite(tk)
        if ok.sum() < 10:
            continue
        tmap = np.interp(T0, kin.t0[ok], tk[ok], left=np.nan, right=np.nan)
        if mute[0] == "angle":
            a = np.interp(T0, kin.t0[ok], kin.ang[ok, ix], left=np.nan, right=np.nan)
            lv = np.isfinite(tmap) & (a <= mute[1])
        elif mute[0] == "stretch":
            g = np.gradient(tmap, T0)                      # dt/dt0 ; stretch = 1/g - 1
            with np.errstate(invalid="ignore", divide="ignore"):
                lv = np.isfinite(tmap) & (g > 0) & (1.0 / g - 1.0 <= mute[1])
            # keep only the part below the first (deepest-going) violation: one mute time per trace
            bad = np.where(~lv)[0]
            if len(bad):
                lv[:bad.max() + 1] = False
        else:
            lv = np.isfinite(tmap)
        # make the angle mute one time per trace as well (everything above the last dead sample is dead)
        bad = np.where(~lv & (T0 < 1.95))[0]
        if len(bad):
            lv[:bad.max() + 1] = False
        cs = CubicSpline(t, d[ix])
        tm = np.where(lv, tmap, 0.0)
        out[ix] = np.where(lv & (tm <= t[-1]), cs(tm), 0.0)
        live[ix] = lv & (tm <= t[-1])
    return out, live


def stack3(c, live, off):
    S = {}
    for k, (a, b) in RANGES.items():
        m = (off >= a - 0.5) & (off < b - 0.5) if k != "far" else (off >= a - 0.5) & (off <= b + 0.5)
        n = live[m].sum(0)
        S[k] = np.where(n > 0, c[m].sum(0) / np.maximum(n, 1), 0.0)
        S["n_" + k] = n
    return S


def win(a, b):
    return (T0 >= a) & (T0 <= b)


def rms(x, m):
    return float(np.sqrt(np.mean(x[m] ** 2)))


ZW = win(*ZONE)
RW = (win(TLIM[0], EXCL[0]) | win(EXCL[1], TLIM[1]))


def est_rms(S, num="far", den="near"):
    return (rms(S[num], ZW) / rms(S[den], ZW)) / (rms(S[num], RW) / rms(S[den], RW))


def picks(x, rel=0.25, min_sep=0.016, tlim=TLIM):
    i0, i1 = int(round(tlim[0] / DT)), int(round(tlim[1] / DT))
    idx = [i for i in range(i0, i1 + 1) if abs(x[i]) >= abs(x[i - 1]) and abs(x[i]) > abs(x[i + 1]) and x[i] * x[i - 1] > 0]
    if not idx:
        return []
    amax = max(abs(x[i]) for i in idx)
    idx = sorted([i for i in idx if abs(x[i]) >= rel * amax], key=lambda i: -abs(x[i]))
    keep = []
    for i in idx:
        if all(abs(i - j) * DT >= min_sep for j in keep):
            keep.append(i)
    return sorted(keep)


def per_reflector(S, pk, num="far", den="near", hw_den=4, hw_num=12):
    """Two own estimators per picked reflector (pick = sample index on a reference trace):
    signed peak: extremum of the pick's polarity within +-hw on each stack (far may be shifted / stretched);
    envelope peak: maximum of the Hilbert envelope within the same windows."""
    env = {k: np.abs(hilbert(S[k])) for k in (num, den)}
    rows = []
    for i in pk:
        pol = np.sign(S[den][i]) or 1.0
        a_den = float((pol * S[den][i - hw_den:i + hw_den + 1]).max())
        a_num = float((pol * S[num][i - hw_num:i + hw_num + 1]).max())
        e_den = float(env[den][i - hw_den:i + hw_den + 1].max())
        e_num = float(env[num][i - hw_num:i + hw_num + 1].max())
        rows.append(dict(t_ms=round(1000 * T0[i], 1), pk=a_num / a_den if a_den > 0 else np.nan, env=e_num / e_den))
    return rows


def zone_ref(rows, key):
    t = np.array([r["t_ms"] for r in rows]) / 1000.0
    v = np.array([r[key] for r in rows], float)
    z = (t >= ZONE[0]) & (t <= ZONE[1])
    r = (t < EXCL[0]) | (t > EXCL[1])
    return dict(zone=[round(float(a), 3) for a in v[z]], zone_med=float(np.median(v[z])) if z.any() else None,
                ref_med=float(np.median(v[r])) if r.any() else None, n=[int(z.sum()), int(r.sum())],
                ratio=float(np.median(v[z]) / np.median(v[r])) if z.any() and r.any() else None,
                ref_minmax=[round(float(v[r].min()), 3), round(float(v[r].max()), 3)] if r.any() else None)


def copy_fit(par_tr, tr, ic, lag_lo, lag_hi, sign):
    """Copy of the event at ic of par_tr inside tr: lag (s) of the best normalised correlation of the requested
    sign for the parent window -20..+40 ms, and the least-squares amplitude there."""
    a, b = 20, 40
    par = par_tr[ic - a:ic + b + 1]
    pn = float(np.sum(par ** 2))
    best = (-9.0, 0, 0.0)
    for L in range(int(round(lag_lo / DT)), int(round(lag_hi / DT)) + 1):
        if ic + b + 1 + L >= len(tr):
            break
        seg = tr[ic - a + L:ic + b + 1 + L]
        dot = float(np.sum(seg * par))
        cc = dot / np.sqrt(pn * float(np.sum(seg ** 2)) + 1e-300)
        if sign * cc > best[0]:
            best = (sign * cc, L, dot / pn)
    return dict(coef=round(best[2], 3), lag_ms=round(1000 * best[1] * DT, 1), ncc=round(best[0], 3))


def jdump(obj, name):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(obj, f, indent=1, default=conv)
