"""Independent primaries-only modelling for the check of the Sleipner experiment (lens "metode").

Own Zoeppritz (4x4 matrix solve, not the closed form the experiment uses), own ray tracing in the 1D model,
stack reflectivity per interface for near / mid / far, convolution with a Ricker wavelet.
Only model_1d.npz of the experiment is read (Vp, Vs, density, two-way time on 1 m cells).
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DT = 0.001
T_OUT = np.arange(0.0, 2.2001, DT)
ZONE = (1.408, 1.492)
REF_EXCL = (1.380, 1.560)
T_EV = (1.100, 1.900)
OFF = np.arange(150.0, 3150.1, 25.0)
GROUPS = {"nea": (149.0, 761.0), "mid": (761.0, 1763.0), "far": (1763.0, 3151.0)}
W_FULL = dict(nea=0.20, mid=0.32, far=0.45)
MUTE = 48.0


def zoeppritz(vp1, vs1, r1, vp2, vs2, r2, p, full=False):
    """Zoeppritz equations for a P wave incident from medium 1, solved as a 4x4 linear system per element.
    Unknowns: [Rpp, Rps, Tpp, Tps] (displacement amplitudes). Returns complex Rpp, or all four with full=True."""
    vp1, vs1, r1, vp2, vs2, r2, p = [np.asarray(a, float) for a in np.broadcast_arrays(vp1, vs1, r1, vp2, vs2, r2, p)]
    shp = p.shape
    vp1, vs1, r1, vp2, vs2, r2, p = [a.ravel().astype(complex) for a in (vp1, vs1, r1, vp2, vs2, r2, p)]
    s1, j1, s2, j2 = p * vp1, p * vs1, p * vp2, p * vs2          # sines of the four angles
    c1, d1, c2, d2 = [np.sqrt(1 - s * s) for s in (s1, j1, s2, j2)]
    sin2t1, sin2t2 = 2 * s1 * c1, 2 * s2 * c2
    sin2f1, sin2f2 = 2 * j1 * d1, 2 * j2 * d2
    cos2f1, cos2f2 = 1 - 2 * j1 * j1, 1 - 2 * j2 * j2
    M = np.empty((len(p), 4, 4), complex)
    M[:, 0, 0], M[:, 0, 1], M[:, 0, 2], M[:, 0, 3] = -s1, -d1, s2, d2
    M[:, 1, 0], M[:, 1, 1], M[:, 1, 2], M[:, 1, 3] = c1, -j1, c2, -j2
    M[:, 2, 0], M[:, 2, 1] = sin2t1, vp1 / vs1 * cos2f1
    M[:, 2, 2], M[:, 2, 3] = r2 * vs2 ** 2 * vp1 / (r1 * vs1 ** 2 * vp2) * sin2t2, r2 * vs2 * vp1 / (r1 * vs1 ** 2) * cos2f2
    M[:, 3, 0], M[:, 3, 1] = -cos2f1, vs1 / vp1 * sin2f1
    M[:, 3, 2], M[:, 3, 3] = r2 * vp2 / (r1 * vp1) * cos2f2, -r2 * vs2 / (r1 * vp1) * sin2f2
    b = np.stack([s1, c1, sin2t1, cos2f1], axis=1)
    x = np.linalg.solve(M, b[:, :, None])[:, :, 0]
    if full:
        return x.reshape(shp + (4,))
    return x[:, 0].reshape(shp)


def aki_richards(vp1, vs1, r1, vp2, vs2, r2, theta):
    """Linearised P-P reflection (Aki & Richards three-term), for the test of the exact solver at small contrasts."""
    vp, vs, r = 0.5 * (vp1 + vp2), 0.5 * (vs1 + vs2), 0.5 * (r1 + r2)
    dvp, dvs, dr = vp2 - vp1, vs2 - vs1, r2 - r1
    k = (vs / vp) ** 2
    s2 = np.sin(theta) ** 2
    return 0.5 * (1 - 4 * k * s2) * dr / r + dvp / vp / (2 * np.cos(theta) ** 2) - 4 * k * s2 * dvs / vs


def self_test(verbose=True):
    out = {}
    rng = np.random.default_rng(1)
    # 1. normal incidence equals the impedance contrast
    a = dict(vp1=2000.0, vs1=750.0, r1=1.9, vp2=2200.0, vs2=830.0, r2=2.0)
    r0 = zoeppritz(p=0.0, **a).real
    out["normal_incidence_err"] = float(abs(r0 - (2200 * 2.0 - 2000 * 1.9) / (2200 * 2.0 + 2000 * 1.9)))
    # 2. energy flux conserved below the critical angle (random media, 2000 cases)
    n = 2000
    vp1 = rng.uniform(1700, 3000, n); vs1 = vp1 / rng.uniform(1.6, 3.5, n); r1 = rng.uniform(1.8, 2.5, n)
    vp2 = rng.uniform(1700, 3000, n); vs2 = vp2 / rng.uniform(1.6, 3.5, n); r2 = rng.uniform(1.8, 2.5, n)
    th = np.radians(rng.uniform(0, 60, n))
    p = np.sin(th) / vp1
    ok = p * np.maximum(vp1, vp2) < 0.98
    x = zoeppritz(vp1, vs1, r1, vp2, vs2, r2, p, full=True)
    c1 = np.sqrt(1 - (p * vp1) ** 2); d1 = np.sqrt(1 - (p * vs1) ** 2)
    c2 = np.sqrt(np.maximum(1 - (p * vp2) ** 2, 0)); d2 = np.sqrt(1 - (p * vs2) ** 2)
    E = (np.abs(x[:, 0]) ** 2 + np.abs(x[:, 1]) ** 2 * vs1 * d1 / (vp1 * c1) + np.abs(x[:, 2]) ** 2 * r2 * vp2 * c2 / (r1 * vp1 * c1)
         + np.abs(x[:, 3]) ** 2 * r2 * vs2 * d2 / (r1 * vp1 * c1))
    out["energy_flux_max_err"] = float(np.max(np.abs(E[ok] - 1)))
    out["energy_cases"] = int(ok.sum())
    # 3. small contrasts: exact against Aki & Richards (error must be second order)
    vp1 = rng.uniform(1900, 2400, n); vs1 = vp1 / rng.uniform(2.2, 3.0, n); r1 = rng.uniform(1.85, 2.1, n)
    e = 0.02
    vp2 = vp1 * (1 + e * rng.uniform(-1, 1, n)); vs2 = vs1 * (1 + e * rng.uniform(-1, 1, n)); r2 = r1 * (1 + e * rng.uniform(-1, 1, n))
    th = np.radians(rng.uniform(0, 45, n))
    ex = zoeppritz(vp1, vs1, r1, vp2, vs2, r2, np.sin(th) / vp1).real
    # the linearised formula uses the mean angle
    th2 = np.arcsin(np.clip(np.sin(th) * vp2 / vp1, -1, 1))
    ar = aki_richards(vp1, vs1, r1, vp2, vs2, r2, 0.5 * (th + th2))
    out["aki_richards_max_abs_diff_2pct_contrast"] = float(np.max(np.abs(ex - ar)))
    out["aki_richards_typical_R"] = float(np.median(np.abs(ex)))
    # 4. symmetry: swapping the media at normal incidence flips the sign
    out["swap_sign_err"] = float(abs(zoeppritz(2200.0, 830.0, 2.0, 2000.0, 750.0, 1.9, 0.0).real + r0))
    # 5. against the experiment's closed-form routine (independent implementation), if slp.py is beside this file
    try:
        import slp
        n = 5000
        vp1 = rng.uniform(1700, 3000, n); vs1 = vp1 / rng.uniform(1.6, 3.5, n); r1 = rng.uniform(1.8, 2.5, n)
        vp2 = rng.uniform(1700, 3000, n); vs2 = vp2 / rng.uniform(1.6, 3.5, n); r2 = rng.uniform(1.8, 2.5, n)
        p = np.sin(np.radians(rng.uniform(0, 70, n))) / vp1
        mine = zoeppritz(vp1, vs1, r1, vp2, vs2, r2, p).real
        theirs = slp.zoep_pp(vp1, vs1, r1, vp2, vs2, r2, p)
        out["max_abs_diff_vs_experiment_closed_form"] = float(np.max(np.abs(mine - theirs)))
        out["n_postcritical_in_that_test"] = int((p * vp2 > 1).sum())
    except Exception as ex_:                                   # noqa
        out["max_abs_diff_vs_experiment_closed_form"] = "slp.py not importable: %s" % ex_
    if verbose:
        print(json.dumps(out, indent=1))
    return out


# ------------------------------------------------------------------------------------------------ model and rays
def load_model(path=None):
    M = np.load(path or os.path.join(HERE, "model_1d.npz"))
    return {k: M[k].copy() for k in ("z", "ze", "vp", "vs", "rho", "twt", "soft")}


def trace_rays(vp, dz=1.0, offsets=OFF, n_p=4000):
    """Ray parameter, incidence angle and travel time of the primary P reflection from the bottom of every cell,
    for source and receiver at z = 0. Straight ray segments in 1 m cells, no smoothing.
    Returns P[k, ix], ANG[k, ix] (degrees, in the cell above the interface), T[k, ix]."""
    n = len(vp)
    pmax = 0.99999 / vp.min()
    p = pmax * np.sin(np.linspace(0, np.pi / 2, n_p))
    P = np.full((n, len(offsets)), np.nan)
    T = np.full((n, len(offsets)), np.nan)
    X = np.zeros(n_p)
    TT = np.zeros(n_p)
    alive = np.ones(n_p, bool)
    for k in range(n):
        pv = p * vp[k]
        alive &= pv < 0.9999
        c = np.sqrt(np.where(alive, 1 - pv * pv, 1.0))
        X = X + 2 * dz * pv / c
        TT = TT + 2 * dz / (vp[k] * c)
        if alive.sum() > 5:
            P[k] = np.interp(offsets, X[alive], p[alive], right=np.nan)
            T[k] = np.interp(offsets, X[alive], TT[alive], right=np.nan)
    with np.errstate(invalid="ignore"):
        ANG = np.degrees(np.arcsin(np.clip(P * vp[:, None], 0, 1)))
    return P, ANG, T


class Forward:
    """Stack reflectivity for near / mid / far and synthetic stacks; fast re-evaluation when only Vs changes."""

    def __init__(self, m, kmin=300, f0=25.0, mute=MUTE):
        self.m = m
        self.vp, self.rho, self.vs0 = m["vp"], m["rho"], m["vs"].copy()
        self.n = len(self.vp)
        self.kmin = kmin                                  # first interface used (below the sea bed and shallow section)
        self.P, self.ANG, self.T = trace_rays(self.vp)
        self.t0 = m["twt"][1:]                             # two-way time of the bottom of each cell
        self.f0 = f0
        k = np.arange(self.n - 1)
        # mute on the angle computed with a smoothed velocity (25 m running mean of the slowness): with the raw 1 m
        # velocity the angle jumps from cell to cell, and near the mute single reflectors would drop out of a trace
        ks = np.ones(25) / 25
        vsm = 1.0 / np.convolve(np.pad(1.0 / self.vp, 12, mode="edge"), ks, mode="valid")
        with np.errstate(invalid="ignore"):
            self.ANGS = np.degrees(np.arcsin(np.clip(self.P * vsm[:, None], 0, 1)))
        self.live = np.isfinite(self.P[:-1]) & (self.ANGS[:-1] <= mute) & (k[:, None] >= kmin)
        self.use = np.isfinite(self.P[:-1]) & (self.ANG[:-1] < 75.0) & (k[:, None] >= kmin)      # reflectors that are modelled
        self.gm = {g: (OFF > a) & (OFF <= b) for g, (a, b) in GROUPS.items()}
        self.R0 = self._refl(self.vs0, np.arange(self.n - 1))
        self._base = None

    def _refl(self, vs, kk):
        """R[kk, ix] of the interfaces kk (between cell k and k+1) at all offsets; 0 where not live."""
        kk = np.asarray(kk)
        P = np.nan_to_num(self.P[kk])
        out = np.zeros(P.shape)
        lv = self.use[kk]
        if lv.any():
            i, j = np.nonzero(lv)
            k = kk[i]
            out[i, j] = zoeppritz(self.vp[k], vs[k], self.rho[k], self.vp[k + 1], vs[k + 1], self.rho[k + 1], P[i, j]).real
        return out

    def refl(self, fac=None):
        """Full R[k, ix] for Vs = Vs0 * fac (fac per cell); only interfaces touched by fac != 1 are recomputed."""
        if fac is None:
            return self.R0
        ch = np.nonzero(np.abs(fac - 1.0) > 1e-12)[0]
        if len(ch) == 0:
            return self.R0
        kk = np.unique(np.clip(np.r_[ch - 1, ch], 0, self.n - 2))
        R = self.R0.copy()
        R[kk] = self._refl(self.vs0 * fac, kk)
        return R

    def wavelet(self, f=None, n=301):
        f = f or self.f0
        t = (np.arange(n) - n // 2) * DT
        a = (np.pi * f * t) ** 2
        return (1 - 2 * a) * np.exp(-a)

    def stacks(self, R, stretch=False, f0=None):
        """Near / mid / far / full synthetic stacks on T_OUT.
        stretch False: the stack reflectivity (mean over the live offsets) is put at t0 and convolved with the wavelet
                       (no NMO stretch).  stretch True: every offset trace is made in data time (spikes at the ray
                       time), convolved, and mapped back to t0 with the exact ray time (exact NMO, with its stretch)."""
        w = self.wavelet(f0)
        S = {}
        nt = len(T_OUT)
        t0 = self.t0[:self.n - 1]
        for g, gm in self.gm.items():
            lv = self.live[:, gm]
            cnt = lv.sum(1)
            if not stretch:
                rs = np.where(cnt > 0, (R[:, gm] * lv).sum(1) / np.maximum(cnt, 1), 0.0)
                x = t0 / DT
                i0 = np.floor(x).astype(int)
                fr = x - i0
                ok = (i0 >= 0) & (i0 + 1 < nt)
                spk = np.zeros(nt)
                np.add.at(spk, i0[ok], rs[ok] * (1 - fr[ok]))
                np.add.at(spk, i0[ok] + 1, rs[ok] * fr[ok])
                S[g] = np.convolve(spk, w, mode="same")
            else:
                acc = np.zeros(nt)
                num = np.zeros(nt)
                ntd = 3200
                td = np.arange(ntd) * DT
                for ix in np.nonzero(gm)[0]:
                    l = self.use[:, ix]                      # all modelled reflectors make the trace ...
                    if self.live[:, ix].sum() < 10:
                        continue
                    x = self.T[:self.n - 1, ix][l] / DT
                    i0 = np.floor(x).astype(int)
                    fr = x - i0
                    spk = np.zeros(ntd)
                    ok = i0 + 1 < ntd
                    np.add.at(spk, i0[ok], R[l, ix][ok] * (1 - fr[ok]))
                    np.add.at(spk, i0[ok] + 1, R[l, ix][ok] * fr[ok])
                    tr = np.convolve(spk, w, mode="same")
                    val = np.interp(self.T[:self.n - 1, ix][l], td, tr)
                    t0l = t0[l]
                    # ... and the mute is applied to the samples of the corrected trace (smoothed-angle criterion)
                    inside = np.interp(T_OUT, t0l, self.live[l, ix].astype(float), left=0.0, right=0.0) > 0.5
                    acc[inside] += np.interp(T_OUT[inside], t0l, val)
                    num[inside] += 1
                S[g] = np.where(num > 0, acc / np.maximum(num, 1), 0.0)
        S["ful"] = sum(W_FULL[g] * S[g] for g in W_FULL)
        return S


# ------------------------------------------------------------------------------------------------ measures
def win_masks(zone=ZONE, excl=REF_EXCL, tev=T_EV):
    t = T_OUT
    z = (t >= zone[0]) & (t <= zone[1])
    r = ((t >= tev[0]) & (t < excl[0])) | ((t > excl[1]) & (t <= tev[1]))
    return z, r


def rms(x, w):
    return float(np.sqrt(np.mean(np.asarray(x)[w] ** 2)))


def rms_measure(S, zone=ZONE, excl=REF_EXCL):
    z, r = win_masks(zone, excl)
    fz, fr = rms(S["far"], z) / rms(S["nea"], z), rms(S["far"], r) / rms(S["nea"], r)
    mz, mr = rms(S["mid"], z) / rms(S["nea"], z), rms(S["mid"], r) / rms(S["nea"], r)
    return dict(fn_zone=fz, fn_ref=fr, fn=fz / fr, mn_zone=mz, mn_ref=mr, mn=mz / mr,
                near_zone_over_ref=rms(S["nea"], z) / rms(S["nea"], r), mid_zone_over_ref=rms(S["mid"], z) / rms(S["mid"], r),
                far_zone_over_ref=rms(S["far"], z) / rms(S["far"], r))


def event_measure(S, events=None):
    """The experiment's own per-reflector measure (slp.measure_event / summarise), on a given list of event samples
    or on events picked on this full stack. Returns the summary plus the event list."""
    import slp
    ev = slp.pick_events(S["ful"]) if events is None else list(events)
    E = [slp.measure_event(S, i) for i in ev]
    s = slp.summarise(E)
    return dict(n_zone=s["n_zone"], n_ref=s["n_ref"], fn=s["fn"]["zone_over_ref"], fn_zone=s["fn"]["zone_median"], fn_ref=s["fn"]["ref_median"],
                fn_zone_values=s["fn"]["zone_values"], mn=s["mn"]["zone_over_ref"], pk_fn=s["pk_fn"]["zone_over_ref"], pk_mn=s["pk_mn"]["zone_over_ref"],
                events_ms=[int(round(1000 * i * DT)) for i in ev])


if __name__ == "__main__":
    self_test()
