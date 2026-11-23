#!/usr/bin/env python3
"""Shared pieces for the Sleipner near/mid/far test. Concatenated in front of each fm_sl_*.py script:
    cat sl_common.py fm_sl_model.py | ssh linux-pc '<geoscience python> -I -'
Everything is read-only on the open data; derived arrays live in AVO (outside the app folder).
"""
import json
import os
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026"
AVO = R + "/sleipner/avo"
DT = 0.002
KEYS = ("nea", "mid", "far", "ful")
# Offset ranges (m) and outer mute, worked out from the data in fm_sl_geom.py (fold weights, mute onset, NMO stretch).
X0, X1, X2, X3, MUTE_DEG = 150.0, 761.0, 1763.0, 3150.0, 48.0
WELL_IL, WELL_XL = 57, 143          # array indices of 15/9-13 (inline 1773, crossline 1099)
F0 = 35.0                           # Ricker peak frequency of the unstretched wavelet, Hz


def stacks():
    return {k: np.load(os.path.join(AVO, k + ".npy"), mmap_mode="r") for k in KEYS}


def well():
    W = json.load(open(R + "/sleipner/cubes/sleipner94/wells.json"))
    w = W["wells"][0]
    L = w["logs"]
    t = L["t0"] + L["dt"] * np.arange(L["n"])
    g = lambda k: np.array([np.nan if v is None else v for v in L[k]], float)
    return w, t, {k: g(k) for k in ("md", "gr", "vp", "rhob", "rho_used", "ai")}


class Vel:
    """Layered velocity model (20 ms layers) from the time-depth curve of 15/9-13."""

    def __init__(self, w, dtl=0.02):
        md = np.array(w["td"]["md"], float)
        twt = np.array(w["td"]["twt"], float) / 1000.0
        z = md - w["kb"]
        ok = np.isfinite(z) & np.isfinite(twt)
        self.dtl = dtl
        self.tl = np.arange(0.0, twt[ok].max() + 1e-9, dtl)
        self.zl = np.interp(self.tl, np.r_[0.0, twt[ok]], np.r_[0.0, z[ok]])
        self.h = np.diff(self.zl)
        self.v = 2 * self.h / dtl

    def _n(self, t0):
        return max(1, min(len(self.v), int(round(t0 / self.dtl))))

    def vrms(self, t0):
        n = self._n(t0)
        return float(np.sqrt(np.mean(self.v[:n] ** 2)))

    def offset(self, th_deg, t0):
        n = self._n(t0)
        v, h = self.v[:n], self.h[:n]
        p = np.sin(np.radians(th_deg)) / v[-1]
        if p * v.max() >= 1:
            return np.inf
        return float(np.sum(2 * h * p * v / np.sqrt(1 - (p * v) ** 2)))

    def sin_angle(self, x, t0):
        n = self._n(t0)
        v, h = self.v[:n], self.h[:n]
        lo, hi = 0.0, 0.999999 / v.max()
        p = 0.0
        for _ in range(50):
            p = 0.5 * (lo + hi)
            if np.sum(2 * h * p * v / np.sqrt(1 - (p * v) ** 2)) < x:
                lo = p
            else:
                hi = p
        return min(1.0, p * v[-1])


def live_ranges(vel, t0, x=(X0, X1, X2, X3), mute=MUTE_DEG):
    xm = min(x[3], vel.offset(mute, t0))
    return {"nea": (x[0], min(x[1], xm)), "mid": (x[1], min(x[2], xm)), "far": (x[2], xm)}


def stack_geometry(vel, t0, n=12, x=(X0, X1, X2, X3), mute=MUTE_DEG):
    """For each stack: sines of the incidence angles and NMO stretch factors of n offsets spread over its live range."""
    out = {}
    for k, (a, b) in live_ranges(vel, t0, x, mute).items():
        if b <= a + 1:
            out[k] = (np.array([]), np.array([]))
            continue
        xs = a + (np.arange(n) + 0.5) * (b - a) / n
        s = np.array([vel.sin_angle(xx, t0) for xx in xs])
        st = t0 / np.sqrt(t0 ** 2 + (xs / vel.vrms(t0)) ** 2)
        out[k] = (s, st)
    return out


def zoep(a1, b1, r1, a2, b2, r2, sin_th):
    """Exact P-P reflection coefficient (Zoeppritz, Aki and Richards form) for an array of sin(incidence angle)."""
    p = np.asarray(sin_th, float) / a1
    c1 = np.sqrt(1 - (p * a1) ** 2 + 0j)
    c2 = np.sqrt(1 - (p * a2) ** 2 + 0j)
    d1 = np.sqrt(1 - (p * b1) ** 2 + 0j)
    d2 = np.sqrt(1 - (p * b2) ** 2 + 0j)
    a = r2 * (1 - 2 * b2 ** 2 * p ** 2) - r1 * (1 - 2 * b1 ** 2 * p ** 2)
    b = r2 * (1 - 2 * b2 ** 2 * p ** 2) + 2 * r1 * b1 ** 2 * p ** 2
    c = r1 * (1 - 2 * b1 ** 2 * p ** 2) + 2 * r2 * b2 ** 2 * p ** 2
    d = 2 * (r2 * b2 ** 2 - r1 * b1 ** 2)
    E = b * c1 / a1 + c * c2 / a2
    F = b * d1 / b1 + c * d2 / b2
    G = a - d * c1 / a1 * d2 / b2
    H = a - d * c2 / a2 * d1 / b1
    D = E * F + G * H * p ** 2
    return (((b * c1 / a1 - c * c2 / a2) * F - (a + d * c1 / a1 * d2 / b2) * H * p ** 2) / D).real


def vs_trend(vp):
    """Local mudstone trend measured in the Lark Formation at Nini-4 (808-1700 m): Vs = -699 + 0.702 Vp."""
    return -699.0 + 0.702 * vp


def vs_castagna(vp):
    """Mudrock line of Castagna et al. (1985): Vs = -1172 + 0.862 Vp."""
    return -1172.0 + 0.862 * vp


def ricker(t, f):
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)


def synth(t_ms, vp, rho, vs, vel, t_lo=1150.0, t_hi=1988.0, f0=F0, x=(X0, X1, X2, X3), mute=MUTE_DEG, stretch=True):
    """Synthetic near, mid and far traces on the 2 ms grid. Reflectivity per interface is the mean of the exact
    reflection coefficient over the live offsets of each stack; the wavelet is a Ricker stretched by NMO."""
    n = len(t_ms)
    tr = {k: np.zeros(n) for k in ("nea", "mid", "far")}
    rc = {k: np.zeros(n) for k in ("nea", "mid", "far")}
    geo_t, geo = None, None
    tt = np.arange(n) * DT
    for i in range(n - 1):
        ti = t_ms[i + 1]
        if ti < t_lo or ti > t_hi or not np.isfinite(vp[i] * vp[i + 1] * rho[i] * rho[i + 1] * vs[i] * vs[i + 1]):
            continue
        if geo_t is None or abs(ti - geo_t) >= 10:
            geo_t, geo = ti, stack_geometry(vel, ti / 1000.0, x=x, mute=mute)
        for k in tr:
            s, st = geo[k]
            if len(s) == 0:
                continue
            r = float(np.mean(zoep(vp[i], vs[i], rho[i], vp[i + 1], vs[i + 1], rho[i + 1], s)))
            rc[k][i + 1] = r
            f = f0 * (float(np.mean(st)) if stretch else 1.0)
            j0, j1 = max(0, i + 1 - 40), min(n, i + 1 + 41)
            tr[k][j0:j1] += r * ricker(tt[j0:j1] - tt[i + 1], f)
    return tr, rc


def xc(a, b, maxlag=20):
    """Normalised cross-correlation of b against a; returns (cc at zero lag, best cc, best lag in ms; positive = b later)."""
    a = a - a.mean()
    b = b - b.mean()
    n = len(a)
    den = np.sqrt((a * a).sum() * (b * b).sum()) + 1e-30
    best = (-9, 0)
    c0 = (a * b).sum() / den
    for L in range(-maxlag, maxlag + 1):
        c = ((a[:n - L] * b[L:]).sum() if L >= 0 else (a[-L:] * b[:n + L]).sum()) / den
        if c > best[0]:
            best = (c, L)
    return c0, best[0], best[1] * DT * 1000

