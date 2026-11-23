"""How much does the one-sided mu averaging at the image sea bed (e12) change the surface multiples at large angles?
Small mirrored model: 81 m water, sediment 1800/600/1.95, reflector at 400 m (2200/800/2.10), offsets 0-800 m
(angle in the water up to about 35 deg).  First peg-leg / primary against own theory (-2 R_seabed x 2D spreading),
with the experiment's set-up (model.mu, 'safe_harmonic') and with a mirror-symmetric mu ('arithmetic')."""
import json

import numpy as np
from devito import Eq, Function, Operator, TensorTimeFunction, VectorTimeFunction, configuration, diag, div, grad

configuration["log-level"] = "WARNING"
from examples.seismic import AcquisitionGeometry, Model, Receiver

import e1_lib as L
import e3_zoep as Z

DX, SO, NBL, F0, XS, ZREF = 2.0, 8, 60, 35.0, 200.0, 400.0
OFF = np.arange(0.0, 800.1, 25.0)


def run(mode, refl):
    z = np.arange(0.0, 700.1, DX)
    sed, deep = z >= 81.0, (z >= ZREF) & refl
    vp = np.where(deep, 2.20, np.where(sed, 1.80, 1.48))
    vs = np.where(deep, 0.80, np.where(sed, 0.60, 0.0))
    rho = np.where(deep, 2.10, np.where(sed, 1.95, 1.03))
    n = len(vp)
    col = lambda a: np.r_[a[:0:-1], a]
    vp, vs, rho = col(vp), col(vs), col(rho)
    nx = int(round(1200.0 / DX)) + 1
    tile = lambda a: np.repeat(a[None, :].astype(np.float32), nx, axis=0)
    model = Model(vp=tile(vp), vs=tile(vs), b=tile(1.0 / rho), origin=(0.0, -(n - 1) * DX), shape=(nx, len(vp)), spacing=(DX, DX),
                  space_order=SO, nbl=NBL, bcs="mask", dt=np.float32(0.48 * DX / 2.20))
    xr = XS + OFF
    geo = AcquisitionGeometry(model, np.c_[xr, np.full(len(xr), 8.0)], np.array([[XS, 6.0], [XS, -6.0]]), 0.0, 950.0, f0=F0 / 1000.0, src_type="Ricker")
    v = VectorTimeFunction(name="v", grid=model.grid, space_order=SO, time_order=1)
    tau = TensorTimeFunction(name="t", grid=model.grid, space_order=SO, time_order=1)
    lam, mu, b, damp = model.lam, model.mu, model.b, model.damp
    if mode != "model":
        mu2 = Function(name="mu2", grid=model.grid, space_order=SO, avg_mode=mode, parameter=True)
        mu2.data_with_halo[:] = mu.data_with_halo[:]
        mu = mu2
    s = model.grid.stepping_dim.spacing
    u_v = Eq(v.forward, damp * (v + s * b * div(tau)))
    u_t = Eq(tau.forward, damp * (tau + s * (lam * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).transpose(inner=False)))))
    src = geo.src
    src.data[:, 1] *= -1.0
    eqs = [u_v, u_t] + src.inject(field=tau.forward[0, 0], expr=s * src) + src.inject(field=tau.forward[1, 1], expr=s * src)
    recs = {}
    for k, zz in enumerate((8.0, 0.0)):
        r = Receiver(name="rec%d" % k, grid=model.grid, time_range=geo.time_axis, coordinates=np.c_[xr, np.full(len(xr), zz)])
        eqs += r.interpolate(expr=-(tau[0, 0] + tau[1, 1]) / 2)
        recs[zz] = r
    op = Operator(eqs, subs=model.spacing_map, name="elastic")
    op(dt=float(geo.dt), time_M=geo.nt - 2)
    tt = geo.time_axis.time_values / 1000.0
    t1 = np.arange(0.0, tt[-1], L.DT)
    rs = lambda a: np.array([np.interp(t1, tt, a[:, i]) for i in range(a.shape[1])])
    return t1, rs(np.array(recs[8.0].data, np.float64)), rs(np.array(recs[0.0].data, np.float64))


vcol = np.where(np.arange(int(ZREF)) + 0.5 < 81.0, 1480.0, 1800.0)
offt = np.maximum(OFF, 1.0)
tx0, px0, dx0 = L.raytab(vcol, offt)
tx1, px1, dx1 = L.raytab(vcol, offt, n_water=1)
OUT = {}
for mode in ("model", "arithmetic"):
    t1, a8, a0 = run(mode, True)
    _, b8, b0 = run(mode, False)
    d = a8 - b8
    rows = []
    for x in (0, 200, 400, 500, 600, 700, 800):
        ix = int(np.argmin(np.abs(OFF - x)))
        tp = tx0[-1, ix] + 1.0 / F0 - 0.009
        w = (t1 > tp - 0.04) & (t1 < tp + 0.04)
        ic = int(np.where(w)[0][np.argmax(np.abs(d[ix][w]))])
        lag_th = float(tx1[-1, ix] - tx0[-1, ix])
        m_ = L.copy_fit(d[ix], d[ix], ic, lag_th - 0.012, lag_th + 0.012, -1)
        pm, pp = px1[-1, ix], px0[-1, ix]
        rsb = Z.fluid_solid_rpp(1480.0, 1.03, 1800.0, 600.0, 1.95, pm)
        cm, cp = np.sqrt(1 - (pm * 1480.0) ** 2), np.sqrt(1 - (pp * 1480.0) ** 2)
        th = -2 * rsb * np.sqrt((dx0[-1, ix] * cp * cp) / (dx1[-1, ix] * cm * cm))
        late = t1 > np.hypot(x, 148.0) / 1480.0 + 1.0 / F0 - 0.02
        rows.append(dict(offset=x, angle_in_water_deg=round(float(np.degrees(np.arcsin(pm * 1480.0))), 1), measured=m_["coef"], theory=round(float(th), 3), measured_over_theory=round(m_["coef"] / float(th), 3),
                         lag_ms=m_["lag_ms"], theory_lag_ms=round(1000 * lag_th, 1), ncc=float(m_["ncc"]),
                         rms_p0_over_rms_p8_after_seabed=float(np.sqrt(np.mean(a0[ix][late] ** 2) / np.mean(a8[ix][late] ** 2)))))
        print(mode, rows[-1], flush=True)
    OUT[mode] = rows
    OUT[mode + "_traces_reflector_only"] = None
    np.savez_compressed("/workspace/sleipner_verify_etterregning/e13_%s.npz" % mode, t=t1, off=OFF, d_reflector=d.astype(np.float32), p0=a0.astype(np.float32), p8=a8.astype(np.float32))
# direct comparison of the two set-ups: reflector-only data (primary + its surface multiples)
A, B = np.load("/workspace/sleipner_verify_etterregning/e13_model.npz"), np.load("/workspace/sleipner_verify_etterregning/e13_arithmetic.npz")
cmp_ = []
for x in (0, 200, 400, 600, 800):
    ix = int(np.argmin(np.abs(OFF - x)))
    a, b = A["d_reflector"][ix].astype(float), B["d_reflector"][ix].astype(float)
    tp = tx0[-1, ix] + 1.0 / F0
    wp = (A["t"] > tp - 0.05) & (A["t"] < tp + 0.04)
    wm = (A["t"] > tp + 0.06) & (A["t"] < tp + 0.40)
    cmp_.append(dict(offset=x, primary_rms_diff_rel=float(np.sqrt(np.mean((a[wp] - b[wp]) ** 2) / np.mean(b[wp] ** 2))),
                     multiples_rms_diff_rel=float(np.sqrt(np.mean((a[wm] - b[wm]) ** 2) / np.mean(b[wm] ** 2))),
                     multiples_rms_ratio=float(np.sqrt(np.mean(a[wm] ** 2) / np.mean(b[wm] ** 2)))))
    print("model vs arithmetic", cmp_[-1], flush=True)
OUT["model_vs_arithmetic"] = cmp_
del OUT["model_traces_reflector_only"], OUT["arithmetic_traces_reflector_only"]
json.dump(OUT, open("/workspace/sleipner_verify_etterregning/e13_results.json", "w"), indent=1)
print("saved")
