"""Why is the pressure at z = 0 not exactly zero in the mirror set-up (e6)?  Tiny test with the experiment's equations:
mirrored model, sources +1 at z = 6 m and -1 at z = -6 m, receivers at z = 0, +8, -8 m.
Cases: water only / water over elastic sediment (sea bed at 81 m) / water over a fluid 'sediment' (vs = 0).
Also prints which grid points of b, lam, mu the discrete equations read (are the parameters averaged symmetrically
at the staggered points?)."""
import json
import re
import sys

import numpy as np
from devito import Eq, Operator, TensorTimeFunction, VectorTimeFunction, configuration, diag, div, grad

configuration["log-level"] = "WARNING"
from examples.seismic import AcquisitionGeometry, Model, Receiver

DX, SO, NBL, F0, XS = 2.0, 8, 60, 35.0, 200.0
OFF = np.arange(0.0, 500.1, 25.0)
OUT = {}


def run(case, show=False):
    z = np.arange(0.0, 400.1, DX)
    sed = z >= 81.0
    if case == "water":
        vp, vs, rho = np.full(len(z), 1.48), np.zeros(len(z)), np.full(len(z), 1.03)
    elif case == "seabed_elastic":
        vp, vs, rho = np.where(sed, 1.80, 1.48), np.where(sed, 0.60, 0.0), np.where(sed, 1.95, 1.03)
    elif case == "seabed_fluid":
        vp, vs, rho = np.where(sed, 1.80, 1.48), np.zeros(len(z)), np.where(sed, 1.95, 1.03)
    elif case == "seabed_density_only":
        vp, vs, rho = np.full(len(z), 1.48), np.zeros(len(z)), np.where(sed, 1.95, 1.03)
    n = len(vp)
    col = lambda a: np.r_[a[:0:-1], a]
    vp, vs, rho = col(vp), col(vs), col(rho)
    assert np.array_equal(vp, vp[::-1]) and np.array_equal(vs, vs[::-1]) and np.array_equal(rho, rho[::-1])
    nx = int(round(800.0 / DX)) + 1
    tile = lambda a: np.repeat(a[None, :].astype(np.float32), nx, axis=0)
    model = Model(vp=tile(vp), vs=tile(vs), b=tile(1.0 / rho), origin=(0.0, -(n - 1) * DX), shape=(nx, len(vp)), spacing=(DX, DX),
                  space_order=SO, nbl=NBL, bcs="mask", dt=np.float32(0.48 * DX / 1.80))
    xr = XS + OFF
    geo = AcquisitionGeometry(model, np.c_[xr, np.full(len(xr), 8.0)], np.array([[XS, 6.0], [XS, -6.0]]), 0.0, 600.0, f0=F0 / 1000.0, src_type="Ricker")
    v = VectorTimeFunction(name="v", grid=model.grid, space_order=SO, time_order=1)
    tau = TensorTimeFunction(name="t", grid=model.grid, space_order=SO, time_order=1)
    lam, mu, b, damp = model.lam, model.mu, model.b, model.damp
    s = model.grid.stepping_dim.spacing
    u_v = Eq(v.forward, damp * (v + s * b * div(tau)))
    u_t = Eq(tau.forward, damp * (tau + s * (lam * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).transpose(inner=False)))))
    if show:
        ev, et = u_v.evaluate, u_t.evaluate
        for nm, e in (("v_x update", ev.rhs[0]), ("v_z update", ev.rhs[1]), ("t_xx update", et.rhs[0, 0]), ("t_zz update", et.rhs[1, 1]), ("t_xz update", et.rhs[0, 1])):
            st = str(e)
            acc = {p: sorted(set(re.findall(p + r"\([^()]*\)", st))) for p in ("b", "lam", "mu", "damp")}
            print("PARAM ACCESS", nm, {k: x for k, x in acc.items() if x}, flush=True)
            OUT.setdefault("parameter_access", {})[nm] = {k: x for k, x in acc.items() if x}
        for fn in (model.b, model.lam, model.mu):
            print("PARAM", fn.name, "staggered:", getattr(fn, "staggered", None), "avg_mode:", getattr(fn, "_avg_mode", getattr(fn, "avg_mode", None)), flush=True)
        # symmetry of the stored parameter arrays about the centre node (index NBL + n - 1 in z)
        for fn in (model.b, model.lam, model.mu, model.damp):
            a = np.array(fn.data[:])
            a = a[a.shape[0] // 2]
            print("ARRAY", fn.name, "len", len(a), "max |a - a[::-1]| =", float(np.abs(a - a[::-1]).max()), flush=True)
            OUT.setdefault("array_symmetry", {})[fn.name] = float(np.abs(a - a[::-1]).max())
    src = geo.src
    src.data[:, 1] *= -1.0
    eqs = [u_v, u_t] + src.inject(field=tau.forward[0, 0], expr=s * src) + src.inject(field=tau.forward[1, 1], expr=s * src)
    recs = {}
    for k, zz in enumerate((8.0, 0.0, -8.0)):
        r = Receiver(name="rec%d" % k, grid=model.grid, time_range=geo.time_axis, coordinates=np.c_[xr, np.full(len(xr), zz)])
        eqs += r.interpolate(expr=-(tau[0, 0] + tau[1, 1]) / 2)
        recs[zz] = r
    op = Operator(eqs, subs=model.spacing_map, name="elastic")
    op(dt=float(geo.dt), time_M=geo.nt - 2)
    p8, p0, pm8 = (np.array(recs[zz].data, np.float64) for zz in (8.0, 0.0, -8.0))
    tt = geo.time_axis.time_values
    rows = []
    for x in (0, 100, 200, 300, 400, 500):
        i = int(np.argmin(np.abs(OFF - x)))
        it = int(np.argmax(np.abs(p0[:, i])))
        late = tt > (np.hypot(x, 148.0) / 1.48 - 10.0)                 # from the sea-bed reflection onwards
        rows.append(dict(offset=x, max_p0_over_max_p8=float(np.abs(p0[:, i]).max() / np.abs(p8[:, i]).max()), t_of_max_p0_ms=float(tt[it]),
                         max_antisym_over_max_p8=float(np.abs(p8[:, i] + pm8[:, i]).max() / np.abs(p8[:, i]).max()),
                         rms_p0_over_rms_p8_after_seabed_arrival=float(np.sqrt(np.mean(p0[late, i] ** 2) / np.mean(p8[late, i] ** 2))),
                         rms_antisym_over_rms_p8_after_seabed_arrival=float(np.sqrt(np.mean((p8[late, i] + pm8[late, i]) ** 2) / np.mean(p8[late, i] ** 2)))))
        print(case, rows[-1], flush=True)
    OUT[case] = rows


for j, c in enumerate(("water", "seabed_elastic", "seabed_fluid", "seabed_density_only")):
    run(c, show=(j == 0))
json.dump(OUT, open("/workspace/sleipner_verify_etterregning/e10_results.json", "w"), indent=1)
print("saved")
