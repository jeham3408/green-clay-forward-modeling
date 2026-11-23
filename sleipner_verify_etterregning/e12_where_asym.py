"""Where does the mirror antisymmetry break when the sea bed is elastic?  Tiny mirrored model, snapshot at the end:
asymmetry of every field about z = 0 and the depth where it is largest; three ways of averaging mu at the
staggered t_xz points (the model's own 'safe_harmonic', 'arithmetic', 'harmonic')."""
import json
import re
import sys

import numpy as np
from devito import Eq, Function, Operator, TensorTimeFunction, VectorTimeFunction, configuration, diag, div, grad

configuration["log-level"] = "WARNING"
from examples.seismic import AcquisitionGeometry, Model, Receiver

DX, SO, NBL, F0, XS = 2.0, 8, 40, 35.0, 150.0
OFF = np.arange(0.0, 300.1, 25.0)
OUT = {}


def run(mode, tn):
    z = np.arange(0.0, 240.1, DX)
    sed = z >= 81.0
    vp, vs, rho = np.where(sed, 1.80, 1.48), np.where(sed, 0.60, 0.0), np.where(sed, 1.95, 1.03)
    n = len(vp)
    col = lambda a: np.r_[a[:0:-1], a]
    vp, vs, rho = col(vp), col(vs), col(rho)
    nx = int(round(500.0 / DX)) + 1
    tile = lambda a: np.repeat(a[None, :].astype(np.float32), nx, axis=0)
    model = Model(vp=tile(vp), vs=tile(vs), b=tile(1.0 / rho), origin=(0.0, -(n - 1) * DX), shape=(nx, len(vp)), spacing=(DX, DX),
                  space_order=SO, nbl=NBL, bcs="mask", dt=np.float32(0.48 * DX / 1.80))
    xr = XS + OFF
    geo = AcquisitionGeometry(model, np.c_[xr, np.full(len(xr), 8.0)], np.array([[XS, 6.0], [XS, -6.0]]), 0.0, tn, f0=F0 / 1000.0, src_type="Ricker")
    v = VectorTimeFunction(name="v", grid=model.grid, space_order=SO, time_order=1)
    tau = TensorTimeFunction(name="t", grid=model.grid, space_order=SO, time_order=1)
    lam, mu, b, damp = model.lam, model.mu, model.b, model.damp
    if mode != "model":
        mu2 = Function(name="mu2", grid=model.grid, space_order=SO, avg_mode=mode, parameter=True)
        mu2.data_with_halo[:] = mu.data_with_halo[:]
        mu = mu2
    for fn in (b, lam, mu):
        a = np.array(fn.data[:])
        OUT.setdefault(mode, {})["array_asym_" + fn.name] = float(np.abs(a - a[:, ::-1]).max())
    s = model.grid.stepping_dim.spacing
    u_v = Eq(v.forward, damp * (v + s * b * div(tau)))
    u_t = Eq(tau.forward, damp * (tau + s * (lam * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).transpose(inner=False)))))
    st = str(u_t.evaluate.rhs[0, 1])
    i = st.find("mu")
    print(mode, "t_xz update, how mu enters:", st[max(0, i - 200):i + 500].replace("\n", " "), flush=True)
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
    res = OUT.setdefault(mode, {})
    res["p0_over_p8_max_by_offset"] = {str(int(x)): float(np.abs(p0[:, i]).max() / np.abs(p8[:, i]).max()) for i, x in enumerate(OFF) if x in (0, 100, 200, 300)}
    # snapshot: the last written time buffer
    c = NBL + n - 1                                          # index of z = 0 along the second axis
    tlast = (geo.nt - 1) % 2
    def snap(f):
        return np.array(f.data[tlast], np.float64)
    fields = {"t_xx (odd, nodes)": (snap(tau[0, 0]), "node", -1), "t_zz (odd, nodes)": (snap(tau[1, 1]), "node", -1), "v_x (odd, nodes in z)": (snap(v[0]), "node", -1),
              "v_z (even, half nodes)": (snap(v[1]), "half", +1), "t_xz (even, half nodes)": (snap(tau[0, 1]), "half", +1)}
    K = 100
    for nm, (a, kind, par) in fields.items():
        if kind == "node":
            up, dn = a[:, c + 1:c + 1 + K], a[:, c - 1:c - 1 - K:-1]
            zc = (np.arange(K) + 1) * DX
        else:
            up, dn = a[:, c:c + K], a[:, c - 1:c - 1 - K:-1]
            zc = (np.arange(K) + 0.5) * DX
        d = up - par * dn                                    # zero if the field has the required parity
        prof = np.abs(d).max(axis=0)
        j = int(np.argmax(prof))
        res[nm] = dict(max_violation_over_field_max=float(np.abs(d).max() / (np.abs(a).max() + 1e-300)), depth_of_max_violation_m=float(zc[j]),
                       violation_in_water_0_79m=float(prof[zc < 79.5].max() / (np.abs(a).max() + 1e-300)), violation_at_81m_rows=float(prof[(zc > 79.5) & (zc < 83.5)].max() / (np.abs(a).max() + 1e-300)))
        print(mode, nm, res[nm], flush=True)
    print(mode, "p0/p8:", res["p0_over_p8_max_by_offset"], {k: v for k, v in res.items() if k.startswith("array")}, flush=True)


tn = float(sys.argv[1]) if len(sys.argv) > 1 else 220.0
for mode in ("model", "arithmetic", "harmonic"):
    try:
        run(mode, tn)
    except Exception as e:
        print(mode, "FAILED", repr(e)[:300], flush=True)
json.dump(OUT, open("/workspace/sleipner_verify_etterregning/e12_results.json", "w"), indent=1)
print("saved")
