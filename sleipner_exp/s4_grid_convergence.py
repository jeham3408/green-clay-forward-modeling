"""Grid convergence: the stacks of the with-surface trend case on the 2 m grid against the 1 m grid (stacks.npz).

The source is injected per grid node, so its strength scales with the cell area: the 2 m data are 4 times the 1 m
data. That constant cancels in every ratio; it is removed here by a least-squares scalar before the difference is
measured (NRMS in per cent, 1100-1900 ms and in the zone), without and with the best common time shift.
"""
import json
import os

import numpy as np

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
Z = np.load(os.path.join(ROOT, "stacks.npz"), allow_pickle=True)
t = np.arange(0, 2.2001, 0.001)
w = (t >= 1.1) & (t <= 1.9)
zone = (t >= 1.408) & (t <= 1.492)
rms = lambda x: float(np.sqrt(np.mean(x ** 2)))
out = {}
for v in ("main", "w35"):
    for k in ("nea", "mid", "far"):
        A = Z["FS_trend_%s__%s" % (v, k)].astype(float)
        B = Z["FSdx1_trend_%s__%s" % (v, k)].astype(float)
        best = None
        for sh in np.arange(-2.0, 2.01, 0.1):
            Bs = np.interp(t, t + sh / 1000.0, B)
            a = np.dot(A[w], Bs[w]) / np.dot(Bs[w], Bs[w])
            n = 200 * rms(A[w] - a * Bs[w]) / (rms(A[w]) + rms(a * Bs[w]))
            if best is None or n < best[0]:
                best = (n, sh)
        a0 = np.dot(A[w], B[w]) / np.dot(B[w], B[w])
        r = A - a0 * B
        out["%s_%s" % (v, k)] = dict(scale_2m_over_1m=round(float(a0), 4),
                                     corr=round(float(np.corrcoef(A[w], B[w])[0, 1]), 4),
                                     nrms_pct_scaled=round(200 * rms(r[w]) / (rms(A[w]) + rms(a0 * B[w])), 2),
                                     nrms_pct_scaled_zone=round(200 * rms(r[zone]) / (rms(A[zone]) + rms(a0 * B[zone])), 2),
                                     best_shift_ms=round(float(best[1]), 2), nrms_pct_scaled_shifted=round(best[0], 2))
        print(v, k, out["%s_%s" % (v, k)])
json.dump(out, open(os.path.join(ROOT, "grid_convergence_detail.json"), "w"), indent=1)
