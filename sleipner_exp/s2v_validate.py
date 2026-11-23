"""Validation of the free surface by the method of images (FS) and of the ghost-only set-up (G) on a small model:
81 m water over sediment (1800 / 600 m/s, 1.95) with one reflector at 800 m (2200 / 800, 2.10).

Expected (plane wave, checked earlier with an acoustic image-row surface and with ray theory): the first-order
water-layer peg-leg is a copy of the primary one water period later (109.5 ms at zero offset), amplitude
-2 R_seabed times the 2D spreading factor = -0.757 at zero offset, about -0.77 at 400 m. In G there must be no copy.
"""
import json
import os

import numpy as np

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
ld = lambda c, k: np.load(os.path.join(ROOT, "gather_%s_%s.npz" % (c, k)))
out = {}
for cfg in ("FS", "G"):
    a, b = ld(cfg, "val1"), ld(cfg, "val0")
    d = a["d"].astype(float) - b["d"].astype(float)          # everything that involves the deep reflector
    t, off = a["t"], a["offsets"]
    dt = t[1] - t[0]
    res = []
    for x in (0.0, 200.0, 400.0, 600.0):
        ix = int(np.argmin(np.abs(off - x)))
        ps = np.linspace(0, 0.999 / 2.2e3, 20000)
        q = lambda v: np.sqrt(1 - (ps * v) ** 2)
        X = (75 + 73) * ps * 1480 / q(1480) + 2 * 719 * ps * 1800 / q(1800)
        T = (75 + 73) / (1480 * q(1480)) + 2 * 719 / (1800 * q(1800))
        tp = float(np.interp(x, X, T))
        tr = d[ix]
        w1 = (t > tp - 0.03) & (t < tp + 0.05)
        i1 = np.where(w1)[0][np.argmax(np.abs(tr[w1]))]
        par = tr[i1 - 30:i1 + 31]
        best = (0.0, 0)
        for L in range(int(0.085 / dt), int(0.125 / dt)):
            cc = float(np.sum(tr[i1 - 30 + L:i1 + 31 + L] * par) / np.sum(par ** 2))
            if abs(cc) > abs(best[0]):
                best = (cc, L)
        res.append(dict(offset=x, t_primary_ray_ms=round(1000 * tp, 1), t_peak_ms=round(1000 * t[i1], 1), peak=float(tr[i1]),
                        copy_over_primary=round(best[0], 3), lag_ms=round(1000 * best[1] * dt, 1)))
    out[cfg] = dict(table=res, info=json.loads(str(a["info"])))
# primary in FS against primary in G (same ghosts): must be the same pulse
a = ld("FS", "val1")["d"].astype(float) - ld("FS", "val0")["d"].astype(float)
b = ld("G", "val1")["d"].astype(float) - ld("G", "val0")["d"].astype(float)
t = ld("FS", "val1")["t"]
w = (t > 0.86) & (t < 0.96)
out["primary_FS_vs_G_zero_offset"] = dict(corr=float(np.corrcoef(a[0][w], b[0][w])[0, 1]), amp_ratio=float(np.abs(a[0][w]).max() / np.abs(b[0][w]).max()))
out["expected"] = dict(copy_zero_offset=-0.757, lag_ms=109.5, note="-2 R_seabed (0.394) x 2D spreading 0.96; G: no copy")
json.dump(out, open(os.path.join(ROOT, "validation_images.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
