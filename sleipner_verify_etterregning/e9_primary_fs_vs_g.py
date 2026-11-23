"""Is the primary (with its ghosts) the same with and without the sea surface at ALL offsets of the validation model,
and in the full model before the first surface multiple of anything can arrive? (the report checks zero offset only)"""
import numpy as np

import e1_lib as L

R = {}
f = {c: (np.load(L.OUT + "/theirs/gather_%s_val1.npz" % c), np.load(L.OUT + "/theirs/gather_%s_val0.npz" % c)) for c in ("FS", "G")}
t, off = f["FS"][0]["t"], f["FS"][0]["offsets"]
a = f["FS"][0]["d"].astype(float) - f["FS"][1]["d"].astype(float)
b = f["G"][0]["d"].astype(float) - f["G"][1]["d"].astype(float)
vcol = np.where(np.arange(800) + 0.5 < 81.0, 1480.0, 1800.0)
tx = L.raytab(vcol, np.maximum(off, 1.0))[0][-1]
rows = []
for x in (0, 200, 400, 600, 800, 900):
    ix = int(np.argmin(np.abs(off - x)))
    w = (t > tx[ix] - 0.04) & (t < tx[ix] + 0.05)
    rows.append(dict(offset=x, corr=float(np.corrcoef(a[ix][w], b[ix][w])[0, 1]), amp_ratio_FS_over_G=float(np.abs(a[ix][w]).max() / np.abs(b[ix][w]).max()),
                     rms_diff_rel=float(np.sqrt(np.mean((a[ix][w] - b[ix][w]) ** 2) / np.mean(b[ix][w] ** 2)))))
    print(rows[-1])
R["validation_primary_FS_vs_G"] = rows
# full model: sea-bed reflection window on the raw gathers (direct wave removed), before its first multiple (+109 ms)
wat = {c: L.load_gather(c, "water") for c in ("FS", "G")}
dF, tt, o = L.load_gather("FS", "trend")
dG = L.load_gather("G", "trend")[0]
dF, dG = dF - wat["FS"][0], dG - wat["G"][0]
rows = []
for x in (150, 300, 450):
    ix = int(np.argmin(np.abs(o - x)))
    te = np.hypot(x, 162.0) / 1480.0
    w = (tt > te - 0.03) & (tt < te + 0.04)
    rows.append(dict(offset=x, corr=float(np.corrcoef(dF[ix][w], dG[ix][w])[0, 1]), amp_ratio_FS_over_G=float(np.abs(dF[ix][w]).max() / np.abs(dG[ix][w]).max())))
    print(rows[-1])
R["full_model_seabed_reflection_FS_vs_G"] = rows
L.jdump(R, "e9_results.json")
