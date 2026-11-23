"""Step 2: elastic finite-difference shot gathers (Devito) through the 1D Sleipner model of step 1.

The operator is the velocity-stress scheme of the Troll project (scripts/07_fd_elastic.py: staggered grid, 8th order
in space, "mask" absorbing layers of 120 cells, source injected into the normal stresses, pressure = -(txx+tzz)/2).

Grid (2 m): nodes at z = 0, 2, 4 ... m, so that the sea surface (z = 0, on a node) and the sea bed (81 m, half-way
between two nodes) are honoured. Grid (1 m, convergence check): nodes at z = 0.5, 1.5 ... m = the 1 m cells of the
model; surface and sea bed both lie half-way between nodes. Two configurations, both with the same ghosts:

  FS   reflecting sea surface by the method of images: the model is mirrored about z = 0 and a source of opposite
       sign is added at the mirror position. The surface lies in water, so pressure = 0 at z = 0 exactly.
       Data: primaries + ghosts + all sea-surface multiples.
  G    no reflecting surface (water continues upwards into an absorbing layer), but image source (opposite sign)
       and image receivers (subtracted). Data: primaries + the same ghosts, NO sea-surface multiples.

Cases: trend | m07 | m15 | m20 (Vs in the soft layer on trend, 7, 15, 20 % lower) | trend_sbl | trend_sbh (density of
the unlogged sediment 81-470 m set to 1.80 / 2.10 instead of 1.95 g/cm3: sea-bed reflection 0.37 / 0.43 instead of
0.40) | sb07 | sb15 (Vs 7 / 15 % lower in the soft beds of the zone, see softbeds_mask) | water (water only, for removing the direct wave) | val1 | val0 (small validation model with / without a reflector).

Run:  SLEIPNER_EXP_ROOT=... [SLEIPNER_DX=1] python -B s2_fd.py FS trend
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
DX = float(os.environ.get("SLEIPNER_DX", "2.0"))
HALF = abs(DX - 1.0) < 1e-9            # 1 m grid: nodes at (j + 0.5) m; 2 m grid: nodes at 2 j m
SYM = os.environ.get("SLEIPNER_SYMSLIP", "") == "1"   # exactly mirror-symmetric shear modulus at the t_xz points (see run)
SUF = ("" if abs(DX - 2.0) < 1e-9 else "_dx%g" % DX) + ("_sym" if SYM else "")
SO, NBL = 8, 120
F0 = 35.0                      # Hz, Ricker peak frequency of the modelled source (25 Hz data are derived by spectral shaping)
TN = 2530.0                    # ms
XS, ZSRC, ZREC = 300.0, 6.0, 8.0
OFF = np.arange(150.0, 3150.1, 25.0)
XMAX = 3700.0
UP_M = 60.0                    # water above z = 0 in the G configuration
CFL = 0.48                     # dt = CFL * DX / vmax (the von Neumann limit of this scheme is 0.5497)
VSCEN = {"trend": 0.0, "m07": 0.07, "m15": 0.15, "m20": 0.20, "trend_sbl": 0.0, "trend_sbh": 0.0}
SBRHO = {"trend_sbl": 1.80 / 1.95, "trend_sbh": 2.10 / 1.95}
SOFTBEDS = {"sb07": 0.07, "sb15": 0.15}   # Vs lowered in the soft beds of the zone (see softbeds_mask)
ZONE_S = (1.408, 1.492)
DT_OUT = 1.0                   # ms


def to_nodes(x1m, slow=False):
    """1 m cells [i, i+1) -> nodes. 1 m grid: the cells themselves. 2 m grid: nodes at z = 2j, each the mean of the
    two cells in [2j-1, 2j+1) (slowness mean for velocities; a node is water only if both cells are water)."""
    if HALF:
        return np.array(x1m, float)
    n = len(x1m) // 2
    j = np.arange(1, n)
    a, b = x1m[2 * j - 1], x1m[2 * j]
    if slow:
        with np.errstate(divide="ignore"):
            v = np.where((a > 0) & (b > 0), 2.0 / (1.0 / np.maximum(a, 1e-9) + 1.0 / np.maximum(b, 1e-9)), 0.0)
    else:
        v = 0.5 * (a + b)
    return np.r_[x1m[0], v]


def softbeds_mask(M, L=61):
    """The 1 m cells inside the zone (depth limits from the model's time-depth curve) whose acoustic impedance lies
    below its running mean over L m. This is the definition of the independent checker (p1_lowvs.py,
    softbeds_zone_L61): the clay package without its thin hard beds, plus one soft bed 30 m deeper."""
    ai = M["vp"] * M["rho"]
    rm = np.convolve(np.pad(ai, L // 2, mode="edge"), np.ones(L) / L, mode="valid")[:len(ai)]
    za, zb = [float(np.interp(t, M["twt"], M["ze"])) for t in ZONE_S]
    return (M["z"] >= za) & (M["z"] < zb) & (rm - ai > 0)


def node_z(n):
    return (np.arange(n) + 0.5) * DX if HALF else np.arange(n) * DX


def column(case):
    """Model column on the nodes (km/s, g/cm3) for a case, and the time step (ms)."""
    if case in ("val1", "val0"):
        z = node_z(int(round(1100.0 / DX)) + (0 if HALF else 1))
        vp = np.where(z < 81.0, 1.48, np.where((z < 800.0) | (case == "val0"), 1.80, 2.20))
        vs = np.where(z < 81.0, 0.0, np.where((z < 800.0) | (case == "val0"), 0.60, 0.80))
        rho = np.where(z < 81.0, 1.03, np.where((z < 800.0) | (case == "val0"), 1.95, 2.10))
        return vp, vs, rho, CFL * DX / 2.20
    M = np.load(os.path.join(ROOT, "model_1d.npz"))
    vp1, vs1, rho1, soft, z1 = M["vp"].copy(), M["vs"].copy(), M["rho"].copy(), M["soft"], M["z"]
    dt = CFL * DX / (vp1.max() / 1000.0)
    if case == "water":
        return None, None, None, dt
    vs1[soft] *= (1.0 - VSCEN.get(case, 0.0))
    if case in SOFTBEDS:
        vs1[softbeds_mask(M)] *= (1.0 - SOFTBEDS[case])
    if case in SBRHO:
        rho1[(z1 > 81.0) & (z1 < 470.0)] *= SBRHO[case]
    return to_nodes(vp1, True) / 1000.0, to_nodes(vs1, True) / 1000.0, to_nodes(rho1), dt


def run(cfg, case):
    from devito import Eq, Operator, TensorTimeFunction, VectorTimeFunction, configuration, diag, div, grad
    configuration["log-level"] = "WARNING"
    from examples.seismic import AcquisitionGeometry, Model, Receiver

    val = case in ("val1", "val0")
    vp, vs, rho, dt = column(case)
    xmax = 1500.0 if val else XMAX
    off = np.arange(0.0, 900.1, 25.0) if val else OFF
    tn = 1400.0 if val else TN
    if case == "water":                                   # deep enough that nothing returns from the base within the record
        nzw = int(round((1900.0 if cfg == "FS" else 2000.0) / DX)) + (0 if HALF else 1)
        vp, vs, rho = np.full(nzw, 1.48), np.zeros(nzw), np.full(nzw, 1.03)
    n = len(vp)
    if cfg == "FS":                                       # mirror about z = 0
        col = (lambda a: np.r_[a[::-1], a]) if HALF else (lambda a: np.r_[a[:0:-1], a])
        z0 = -(n - 0.5) * DX if HALF else -(n - 1) * DX
    else:                                                 # water above z = 0
        nup = int(round(UP_M / DX))
        col = lambda a: np.r_[np.full(nup, a[0]), a]
        z0 = -(nup - 0.5) * DX if HALF else -nup * DX
    vp, vs, rho = col(vp), col(vs), col(rho)
    nx = int(round(xmax / DX)) + 1
    tile = lambda a: np.repeat(a[None, :].astype(np.float32), nx, axis=0)
    model = Model(vp=tile(vp), vs=tile(vs), b=tile(1.0 / rho), origin=(0.0, z0), shape=(nx, len(vp)), spacing=(DX, DX),
                  space_order=SO, nbl=NBL, bcs="mask", dt=np.float32(dt))
    xr = XS + off
    src_xy = np.array([[XS, ZSRC], [XS, -ZSRC]])
    geo = AcquisitionGeometry(model, np.c_[xr, np.full(len(xr), ZREC)], src_xy, 0.0, tn, f0=F0 / 1000.0, src_type="Ricker")
    step = float(geo.dt)
    v = VectorTimeFunction(name="v", grid=model.grid, space_order=SO, time_order=1)
    tau = TensorTimeFunction(name="t", grid=model.grid, space_order=SO, time_order=1)
    lam, mu, b, damp = model.lam, model.mu, model.b, model.damp
    s = model.grid.stepping_dim.spacing
    u_v = Eq(v.forward, damp * (v + s * b * div(tau)))
    u_t = Eq(tau.forward, damp * (tau + s * (lam * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).transpose(inner=False)))))
    eq_t = [u_t]
    if SYM:
        # Devito averages mu to the staggered t_xz points with a one-sided guard against mu = 0, so the mirrored model
        # is not an exact mirror image at the image of the sea bed (found by the independent checker, e10-e14).
        # Here mu at the t_xz points is set explicitly: harmonic mean of the four neighbours, zero if any is fluid.
        # That is exactly symmetric and gives slip at both sea beds, so the pressure at z = 0 is zero.
        from devito import Function
        xd, yd = model.grid.dimensions
        muxz = Function(name="muxz", grid=model.grid, space_order=SO, staggered=(xd, yd), parameter=True)
        m_ = np.array(mu.data[:], np.float64)
        q = np.stack([m_[:-1, :-1], m_[1:, :-1], m_[:-1, 1:], m_[1:, 1:]])
        with np.errstate(divide="ignore"):
            hm = np.where((q > 0).all(axis=0), 4.0 / np.sum(1.0 / np.where(q > 0, q, 1.0), axis=0), 0.0)
        full = np.zeros(m_.shape)
        full[:-1, :-1] = hm
        full[-1, :], full[:, -1] = full[-2, :], full[:, -2]
        muxz.data[:] = full
        E, dv = grad(v.forward) + grad(v.forward).transpose(inner=False), div(v.forward)
        eq_t = [Eq(tau[0, 0].forward, damp * (tau[0, 0] + s * (lam * dv + mu * E[0, 0]))),
                Eq(tau[1, 1].forward, damp * (tau[1, 1] + s * (lam * dv + mu * E[1, 1]))),
                Eq(tau[0, 1].forward, damp * (tau[0, 1] + s * muxz * E[0, 1]))]
    src = geo.src
    src.data[:, 1] *= -1.0                                # image source: opposite sign
    eqs = [u_v] + eq_t + src.inject(field=tau.forward[0, 0], expr=s * src) + src.inject(field=tau.forward[1, 1], expr=s * src)
    zr = [ZREC] if cfg == "FS" else [ZREC, -ZREC]
    recs = []
    for k, zz in enumerate(zr):
        r = Receiver(name="rec%d" % k, grid=model.grid, time_range=geo.time_axis, coordinates=np.c_[xr, np.full(len(xr), zz)])
        eqs += r.interpolate(expr=-(tau[0, 0] + tau[1, 1]) / 2)
        recs.append(r)
    op = Operator(eqs, subs=model.spacing_map, name="elastic")
    t0 = time.time()
    op(dt=step, time_M=geo.nt - 2)
    sec = time.time() - t0
    d = np.array(recs[0].data, np.float64)
    if cfg == "G":
        d = d - np.array(recs[1].data, np.float64)       # receiver ghost from the image receivers
    tfd = geo.time_axis.time_values
    # the source is injected as stress RATE: integrate in time and flip the sign to get the pressure pulse of a
    # pressure source (as the Troll project does in 09_fd_figures.py); a 3 Hz low-cut removes the drift
    from scipy.signal import butter, filtfilt
    d = -np.cumsum(d, axis=0) * (step / 1000.0)
    bb, aa = butter(2, 3.0 / (0.5 * 1000.0 / step), "high")
    d = filtfilt(bb, aa, d, axis=0)
    t_out = np.arange(0.0, tn - 1000.0 / F0, DT_OUT)
    out = np.empty((len(xr), len(t_out)), np.float32)
    for i in range(len(xr)):
        out[i] = np.interp(t_out + 1000.0 / F0, tfd, d[:, i])
    info = dict(cfg=cfg, case=case, dx=DX, symslip=SYM, seconds=round(sec, 1), dt_ms=step, nt=int(geo.nt), shape=[int(nx), int(len(vp))], z0=float(z0),
                cells_with_absorbing=int((nx + 2 * NBL) * (len(vp) + 2 * NBL)), finite=bool(np.isfinite(out).all()),
                vp_max=float(vp.max() * 1000), vs_min_solid=float(vs[vs > 0].min() * 1000) if (vs > 0).any() else 0.0,
                courant=float(vp.max() * step / DX), f0=F0, threads=os.environ.get("OMP_NUM_THREADS", ""), host=os.uname().nodename[:12])
    np.savez_compressed(os.path.join(ROOT, "gather_%s_%s%s.npz" % (cfg, case, SUF)), d=out, t=t_out / 1000.0, offsets=off, info=json.dumps(info))
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
