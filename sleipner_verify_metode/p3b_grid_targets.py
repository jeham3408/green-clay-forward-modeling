"""Point 3, continued: the saved 3D grid of multiple strengths (p3_grid3d.npy) against both mid/near observations."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
a = np.load(os.path.join(HERE, "p3_grid3d.npy"))
cols = dict(rms=(3, 4), evG=(5, 6), evOwn=(7, 8), pkG=(9, 10))
out = {}
for q, (i, j) in cols.items():
    fn, mn = a[:, i], a[:, j]
    ok = np.isfinite(fn) & np.isfinite(mn)
    row = {}
    for nm, (flo, fhi, mlo, mhi) in dict(main_method=(0.40, 0.56, 0.90, 1.10), other_method=(0.40, 0.56, 0.60, 0.75), deficit_only=(0.40, 0.56, 0.0, 9.0)).items():
        h = ok & (fn >= flo) & (fn <= fhi) & (mn >= mlo) & (mn <= mhi)
        ordered = h & (a[:, 0] >= a[:, 1]) & (a[:, 1] >= a[:, 2])
        row[nm] = dict(hits=int(h.sum()), hits_near_ge_mid_ge_far=int(ordered.sum()),
                       examples_ordered=[[float(v) for v in r[:3]] + [round(float(r[i]), 2), round(float(r[j]), 2)] for r in a[ordered][:6]])
    # highest mid/near among the combinations that have the far/near deficit
    h = ok & (fn >= 0.40) & (fn <= 0.56)
    row["max_mid_near_given_far_near_0p40_0p56"] = float(mn[h].max()) if h.any() else None
    k = int(np.argmax(np.where(h, mn, -1)))
    row["at_alpha"] = [float(v) for v in a[k, :3]]
    out[q] = row
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, "p3b_results.json"), "w"), indent=1)
