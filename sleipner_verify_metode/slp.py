"""Shared functions for the Sleipner 1D synthetic experiment: ray times, NMO, stacks, Zoeppritz synthetic, measures.

All times in seconds unless a name ends in _ms. The model is model_1d.npz from s1_model.py (1 m cells).
"""
import json
import os

import numpy as np

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
X0, X1, X2, X3 = 150.0, 761.0, 1763.0, 3150.0       # offset limits of near / mid / far (as the parent session's sl_common.py)
MUTE_DEG = 48.0                                      # outer mute as incidence angle at the reflector
W_FULL = dict(nea=0.20, mid=0.32, far=0.45)          # the full stack as weighted sum (parent session, section 2.10)
ZD = 0.0                                             # datum of the ray times: the sea surface. The pulse with its source (6 m)
#                                                      and receiver (8 m) ghosts is centred exactly on the surface-to-surface
#                                                      ray time (arrivals at -14, -2, +2, +14 m of vertical water path around it)
VW = 1480.0
DT = 0.001
T_OUT = np.arange(0.0, 2.2001, DT)
ZONE = (1.408, 1.492)
REF_EXCL = (1.380, 1.560)
T_EV = (1.100, 1.900)
VSCEN = {"trend": 0.0, "m07": 0.07, "m15": 0.15, "m20": 0.20, "m30": 0.30, "m40": 0.40}


def model(case="trend"):
    M = np.load(os.path.join(ROOT, "model_1d.npz"))
    m = {k: M[k].copy() for k in ("z", "ze", "vp", "vs", "rho", "twt", "soft")}
    m["vs"][m["soft"]] *= (1.0 - VSCEN.get(case, 0.0))
    m["info"] = json.loads(str(M["info"]))
    return m


def smooth_slowness(v, n=5):
    k = np.ones(n) / n
    s = np.convolve(np.pad(1.0 / v, n // 2, mode="edge"), k, mode="valid")
    return 1.0 / s


class Rays:
    """Exact reflection times in the flat model for all cell bottoms (reflector depth k+1 m) and offsets.

    t[k, ix]   data time of the primary P reflection (source and receiver at the datum ZD)
    p[k, ix]   ray parameter;  ang[k, ix] incidence angle at the reflector (degrees);  t0[k] two-way time from z = 0
    J[k, ix]   2D spreading term (as the Troll project's 05_raytracing.py)
    """

    def __init__(self, m, offsets, n_p=1600, n_water=0, h_water=81.0):
        """n_water > 0: the same reflection with n_water extra round trips in the water layer (water-layer multiple)."""
        v = smooth_slowness(m["vp"], 5)
        v[m["z"] < 81.0] = VW
        self.v = v
        n = len(v)
        i0 = int(round(ZD))
        pmax = 0.9995 / v[i0:].min()
        self.i0 = i0
        p = pmax * np.sin(np.linspace(0.0, np.pi / 2, n_p))
        pv2 = (p[:, None] * v[None, :]) ** 2
        with np.errstate(invalid="ignore"):
            c = np.sqrt(np.where(pv2 < 1.0, 1.0 - pv2, np.nan))
        dx = p[:, None] * v[None, :] / c
        dt = 1.0 / (v[None, :] * c)
        dx[:, :i0] = 0.0
        dt[:, :i0] = 0.0
        X = np.cumsum(dx, axis=1)
        T = np.cumsum(dt, axis=1)
        h = np.asarray(offsets, float)
        self.t = np.full((n, len(h)), np.nan)
        self.p = np.full((n, len(h)), np.nan)
        self.J = np.full((n, len(h)), np.nan)
        cw = np.sqrt(np.clip(1 - (p * VW) ** 2, 1e-12, None))
        Xw = n_water * 2 * h_water * p * VW / cw
        Tw = n_water * 2 * h_water / (VW * cw)
        for k in range(i0 + 9, n):
            ok = np.isfinite(X[:, k])
            if ok.sum() < 5:
                continue
            Xk, Tk, Pk = 2 * X[ok, k] + Xw[ok], 2 * T[ok, k] + Tw[ok], p[ok]
            self.t[k] = np.interp(h, Xk, Tk, right=np.nan)
            self.p[k] = np.interp(h, Xk, Pk, right=np.nan)
            dXdp = np.gradient(Xk, Pk)
            self.J[k] = np.interp(h, Xk, dXdp, right=np.nan)
        cs = np.sqrt(np.clip(1 - (self.p * VW) ** 2, 0, 1))
        self.J = self.J * cs * cs / VW
        with np.errstate(invalid="ignore"):
            self.ang = np.degrees(np.arcsin(np.clip(self.p * v[:, None], 0, 1)))
        self.t0 = m["twt"][1:]
        self.offsets = np.asarray(offsets, float)
        self.vrms = np.sqrt(np.cumsum(m["vp"] ** 2 * (2.0 / m["vp"])) / np.cumsum(2.0 / m["vp"]))


def nmo(d, t, rays, mode="exact", mute_deg=MUTE_DEG, stretch_max=None):
    """NMO-correct a gather d[ix, it] (time axis t, zero = source pulse peak). Returns (corrected[ix, it0], live[ix, it0]).

    mode "exact": the ray-traced reflection time of each depth (primaries are flat by construction).
    mode "hyper": t = sqrt(t0^2 + x^2 / Vrms(t0)^2), minus the static of the source and receiver depths.
    Mute: incidence angle above mute_deg, and / or wavelet stretch above stretch_max (0.5 = 50 %).
    """
    out = np.zeros((d.shape[0], len(T_OUT)), np.float32)
    live = np.zeros((d.shape[0], len(T_OUT)), bool)
    t0k = rays.t0
    for ix in range(d.shape[0]):
        if mode == "exact":
            tk = rays.t[:, ix]
        else:
            tk = np.sqrt(t0k ** 2 + (rays.offsets[ix] / rays.vrms) ** 2)
        ok = np.isfinite(tk) & np.isfinite(rays.ang[:, ix])
        if ok.sum() < 10:
            continue
        g = np.gradient(np.where(ok, tk, np.nan), t0k)
        with np.errstate(invalid="ignore"):
            if mute_deg is not None:
                ok &= rays.ang[:, ix] <= mute_deg
            if stretch_max is not None:
                ok &= (1.0 / np.maximum(g, 1e-6)) <= 1.0 + stretch_max
            ok &= np.isfinite(g) & (g > 0)
        if ok.sum() < 10:
            continue
        vals = np.interp(tk[ok], t, d[ix], left=0.0, right=0.0)
        out[ix] = np.interp(T_OUT, t0k[ok], vals, left=0.0, right=0.0)
        lv = np.interp(T_OUT, t0k, ok.astype(float), left=0.0, right=0.0) > 0.99
        live[ix] = lv
        out[ix] *= lv
    return out, live


def stacks(c, live, offsets):
    o = np.asarray(offsets)
    S = {}
    for k, (a, b) in (("nea", (X0 - 1, X1)), ("mid", (X1, X2)), ("far", (X2, X3 + 1))):
        m = (o > a) & (o <= b)
        n = live[m].sum(0)
        S[k] = np.where(n > 0, c[m].sum(0) / np.maximum(n, 1), 0.0)
        S["n_" + k] = n
    S["ful"] = sum(W_FULL[k] * S[k] for k in W_FULL)
    return S


def reshape_wavelet(d, f_from, f_to, line_to_point=True):
    """Turn data modelled with a Ricker of peak f_from into data for a Ricker of peak f_to (exact, linear), and remove
    the phase of the 2D line source (the far field of a line source is the half-integral of the pulse; multiplying the
    spectrum by sqrt(i omega) gives the pulse a point source would give; the 2D amplitude decay is left as it is)."""
    n = d.shape[-1]
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfftfreq(nfft, DT)
    H = (f_from / f_to) ** 2 * np.exp(-f ** 2 * (1.0 / f_to ** 2 - 1.0 / f_from ** 2)) if abs(f_from - f_to) > 1e-6 else np.ones(len(f))
    H = H.astype(complex)
    if line_to_point:
        H = H * np.sqrt(1j * 2 * np.pi * f)
    return np.fft.irfft(np.fft.rfft(d, nfft, axis=-1) * H, nfft, axis=-1)[..., :n]


# ---------------------------------------------------------------------------------------- primaries-only synthetic
def zoep_pp(a1, b1, r1, a2, b2, r2, p):
    """Exact P-P reflection coefficient (Aki & Richards) for ray parameter p; arrays broadcast. Real part returned.
    Fluid over solid (b1 = 0) and fluid over fluid are handled separately."""
    a1, b1, r1, a2, b2, r2, p = np.broadcast_arrays(a1, b1, r1, a2, b2, r2, p)
    p = p.astype(complex)
    c1 = np.sqrt(1 - (p * a1) ** 2)
    c2 = np.sqrt(1 - (p * a2) ** 2)
    out = np.zeros(p.shape)
    fluid = b1 < 1.0
    sol = ~fluid
    if sol.any():
        A1, B1, R1, A2, B2, R2, P, C1, C2 = (x[sol] for x in (a1, b1, r1, a2, b2, r2, p, c1, c2))
        d1 = np.sqrt(1 - (P * B1) ** 2)
        d2 = np.sqrt(1 - (P * B2) ** 2)
        a = R2 * (1 - 2 * B2 ** 2 * P ** 2) - R1 * (1 - 2 * B1 ** 2 * P ** 2)
        b = R2 * (1 - 2 * B2 ** 2 * P ** 2) + 2 * R1 * B1 ** 2 * P ** 2
        c = R1 * (1 - 2 * B1 ** 2 * P ** 2) + 2 * R2 * B2 ** 2 * P ** 2
        d = 2 * (R2 * B2 ** 2 - R1 * B1 ** 2)
        E = b * C1 / A1 + c * C2 / A2
        F = b * d1 / B1 + c * d2 / B2
        G = a - d * C1 / A1 * d2 / B2
        H = a - d * C2 / A2 * d1 / B1
        D = E * F + G * H * P ** 2
        out[sol] = (((b * C1 / A1 - c * C2 / A2) * F - (a + d * C1 / A1 * d2 / B2) * H * P ** 2) / D).real
    if fluid.any():
        A1, R1, A2, B2, R2, P, C1, C2 = (x[fluid] for x in (a1, r1, a2, b2, r2, p, c1, c2))
        z1 = R1 * A1 / C1
        cj = np.sqrt(1 - (P * B2) ** 2)
        k = 1 - 2 * (P * B2) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            z2 = R2 * A2 / C2 * k ** 2 + np.where(B2 > 1.0, R2 * B2 / cj * (2 * P * B2 * cj) ** 2, 0.0)
            out[fluid] = ((z2 - z1) / (z2 + z1)).real
    return np.nan_to_num(out)


def ricker(f, n=301):
    t = (np.arange(n) - n // 2) * DT
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


def conv_gather(m, rays, f0, ghosts=True, t_len=2.502, mult=()):
    """Primaries-only gather: exact Zoeppritz P-P at every 1 m interface, exact times, 2D spreading, optional source
    and receiver ghosts (6 m and 8 m). No transmission loss, no multiples, no converted waves."""
    n = len(m["vp"])
    nt = int(round(t_len / DT)) + 1
    d = np.zeros((len(rays.offsets), nt))
    a1, b1, r1 = m["vp"][:-1, None], m["vs"][:-1, None], m["rho"][:-1, None]
    a2, b2, r2 = m["vp"][1:, None], m["vs"][1:, None], m["rho"][1:, None]
    isb = int(np.argmax(m["vs"] > 0))                         # first sediment cell: sea-bed properties (mean of 4 m)
    sb = [float(np.mean(m[k][isb:isb + 4])) for k in ("vp", "vs", "rho")]
    parts = [(1.0, -14.0)] if not ghosts else [(1.0, -14.0), (-1.0, -2.0), (-1.0, 2.0), (1.0, 14.0)]
    # primaries (order 0) and, optionally, water-layer multiples of order nw: amplitude (nw+1) (-R_seabed)^nw (Backus)
    for ry, nw in [(rays, 0)] + list(mult):
        P = ry.p[:n - 1]
        Pz = np.nan_to_num(P)
        R = zoep_pp(a1, b1, r1, a2, b2, r2, Pz)
        if nw:
            Rsb = zoep_pp(VW, 0.0, 1.03, sb[0], sb[1], sb[2], Pz)
            R = R * (nw + 1) * (-Rsb) ** nw
        with np.errstate(invalid="ignore", divide="ignore"):
            A = R / np.sqrt(np.abs(ry.J[:n - 1]))
        T = ry.t[:n - 1]
        ok = np.isfinite(T) & np.isfinite(A) & (ry.ang[:n - 1] < 70.0)
        ok[:isb + 1] = ok[:isb + 1] & (nw == 0)                # no "multiple of the sea bed itself" here (it is muted anyway)
        qw = np.sqrt(np.clip(1.0 / VW ** 2 - Pz ** 2, 0, None))
        for ix in range(d.shape[0]):
            m_ = ok[:, ix]
            for sgn, dz in parts:
                tt = (T[m_, ix] + dz * qw[m_, ix]) / DT
                i0 = np.floor(tt).astype(int)
                fr = tt - i0
                g = (i0 >= 0) & (i0 + 1 < nt)
                np.add.at(d[ix], i0[g], sgn * A[m_, ix][g] * (1 - fr[g]))
                np.add.at(d[ix], i0[g] + 1, sgn * A[m_, ix][g] * fr[g])
    # wavelet: Ricker with the far-field response of a 2D line source (half integration, causal)
    w = ricker(f0)
    nfft = 1 << int(np.ceil(np.log2(nt + len(w))))
    om = 2 * np.pi * np.fft.rfftfreq(nfft, DT)
    Wf = np.fft.rfft(w, nfft) * np.exp(1j * om * (len(w) // 2) * DT)
    out = np.fft.irfft(np.fft.rfft(d, nfft, axis=-1) * Wf[None, :], nfft, axis=-1)[:, :nt]
    return out.astype(np.float32), np.arange(nt) * DT


# ---------------------------------------------------------------------------------------- measures on stacked traces
def pick_events(full, tlim=T_EV, rel=0.25, min_sep=0.016):
    """Local extrema of the full stack inside tlim, at least rel x the largest, at least min_sep apart."""
    i0, i1 = int(tlim[0] / DT), int(tlim[1] / DT)
    x = full
    idx = [i for i in range(i0 + 1, i1 - 1) if (x[i] - x[i - 1]) * (x[i + 1] - x[i]) <= 0 and abs(x[i]) > 0]
    if not idx:
        return []
    amax = max(abs(x[i]) for i in idx)
    idx = [i for i in idx if abs(x[i]) >= rel * amax]
    idx.sort(key=lambda i: -abs(x[i]))
    keep = []
    for i in idx:
        if all(abs(i - j) * DT >= min_sep for j in keep):
            keep.append(i)
    return sorted(keep)


def frac_shift(x, s):
    n = len(x)
    f = np.fft.rfftfreq(n)
    return np.fft.irfft(np.fft.rfft(x) * np.exp(-2j * np.pi * f * s), n)


def lag_to(a, m, c, cw, maxlag):
    """Lag (samples, positive = a later than m) maximising the cross-correlation of windows of half-width cw at c."""
    lags = np.arange(-maxlag, maxlag + 1)
    cf = np.array([np.sum(m[c - cw:c + cw + 1] * a[c - cw + L:c + cw + 1 + L]) for L in lags])
    j = int(np.clip(np.argmax(cf), 1, len(lags) - 2))
    y0, y1, y2 = cf[j - 1], cf[j], cf[j + 1]
    den = y0 - 2 * y1 + y2
    dd = 0.5 * (y0 - y2) / den if abs(den) > 1e-30 else 0.0
    return lags[j] + float(np.clip(dd, -1, 1))


def measure_event(S, ic, cw=16, maxlag=12, pad=64):
    """Amplitude ratios of one event at sample ic, the way the parent session's sl_events_lib.measure does it on the
    cube, here on single traces: align near and far to mid, then
        far/near = sum(far*mid)/sum(near*mid),  mid/near = sum(mid*far)/sum(near*far)   (third stack as instrument)
    and the alignment-free peak ratios (extremum near the centre, polarity of mid)."""
    lo, hi = ic - pad, ic + pad + 1
    seg = {k: np.array(S[k][lo:hi], float) for k in ("nea", "mid", "far")}
    c = pad
    ln = lag_to(seg["nea"], seg["mid"], c, cw, maxlag)
    lf = lag_to(seg["far"], seg["mid"], c, cw, maxlag)
    al = dict(mid=seg["mid"], nea=frac_shift(seg["nea"], -ln), far=frac_shift(seg["far"], -lf))
    w = slice(c - cw, c + cw + 1)
    P = lambda a, b: float(np.sum(al[a][w] * al[b][w]))
    with np.errstate(divide="ignore", invalid="ignore"):
        fn = P("far", "mid") / P("nea", "mid")
        mn = P("far", "mid") / P("nea", "far")
    seg5 = seg["mid"][c - 4:c + 5]
    pol = np.sign(seg5[np.argmax(np.abs(seg5))])
    pk = {k: float((pol * seg[k][c - hw:c + hw + 1]).max()) for k, hw in (("nea", 10), ("mid", 10), ("far", 16))}
    return dict(t_ms=1000 * ic * DT, fn=float(fn), mn=float(mn), lag_near_ms=1000 * DT * ln, lag_far_ms=1000 * DT * lf,
                pk_fn=pk["far"] / pk["nea"] if pk["nea"] > 0 else np.nan, pk_mn=pk["mid"] / pk["nea"] if pk["nea"] > 0 else np.nan,
                amp_near=float(pol * pk["nea"]), amp_mid=float(pol * pk["mid"]), amp_far=float(pol * pk["far"]))


def copy_coef(x, ic, lag_lo, lag_hi, sign=0, y=None):
    """Copy of the event at ic found between lag_lo and lag_hi (seconds), as the parent session's fm_sl_mult.py does it:
    the lag with the best wavelet match (normalised correlation with the parent segment -20..+40 ms, of the expected
    sign if sign = -1 or +1), and there the amplitude relative to the parent (projection on the parent segment).
    y: trace in which the copy is searched (default: the same trace as the parent)."""
    y = x if y is None else y
    a, b = int(0.020 / DT), int(0.040 / DT)
    par = x[ic - a:ic + b + 1]
    pn = float(np.sum(par ** 2))
    best = (-9.0, 0, 0.0)
    for L in range(int(lag_lo / DT), int(lag_hi / DT) + 1):
        if ic + L + b + 1 >= len(y):
            break
        seg = y[ic - a + L:ic + b + 1 + L]
        dot = float(np.sum(seg * par))
        cc = dot / np.sqrt(pn * float(np.sum(seg ** 2)) + 1e-300)
        score = sign * cc if sign else abs(cc)
        if score > best[0]:
            best = (score, L, dot / pn)
    return best[2], best[1] * DT


def three_point_filter(x, a1, lag1, a2, lag2):
    """y(t) = x(t) - a1 x(t - lag1) + (a1^2 - a2) x(t - lag2): removes a first copy a1 and a second copy a2."""
    y = np.array(x, float)
    l1, l2 = int(round(lag1 / DT)), int(round(lag2 / DT))
    y[l1:] += -a1 * x[:-l1]
    y[l2:] += (a1 * a1 - a2) * x[:-l2]
    return y


def summarise(events):
    """Zone against reference from a list of measured events (dicts from measure_event)."""
    t = np.array([e["t_ms"] for e in events]) / 1000.0
    zone = (t >= ZONE[0]) & (t <= ZONE[1])
    ref = ((t < REF_EXCL[0]) | (t > REF_EXCL[1]))
    out = dict(n_zone=int(zone.sum()), n_ref=int(ref.sum()))
    for key in ("fn", "mn", "pk_fn", "pk_mn"):
        v = np.array([e[key] for e in events], float)
        good = np.isfinite(v) & (v > 0) & (v < 10)
        z, r = v[zone & good], v[ref & good]
        out[key] = dict(zone_median=float(np.median(z)) if len(z) else None, ref_median=float(np.median(r)) if len(r) else None,
                        zone_over_ref=float(np.median(z) / np.median(r)) if len(z) and len(r) else None,
                        zone_values=[round(float(a), 3) for a in z], ref_min_max=[round(float(r.min()), 3), round(float(r.max()), 3)] if len(r) else None,
                        ref_sd_log=float(np.std(np.log(r))) if len(r) > 1 else None, n_bad=int((~good).sum()))
    return out
