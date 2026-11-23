"""Grid convergence with own processing: the with-surface trend case on the 1 m grid (gather_FS_trend_dx1.npz,
finished on the rented machine about 21:35 UTC) against the 2 m grid."""
import json

import numpy as np

import e1_lib as L

m = L.load_model()
z1, zw = np.load(L.SRC + "/gather_FS_trend_dx1.npz"), np.load(L.SRC + "/gather_FS_water_dx1.npz")
d, t, off = z1["d"].astype(float) - zw["d"].astype(float), z1["t"].astype(float), z1["offsets"].astype(float)
print("info:", str(z1["info"])[:400])
kin = L.Kin(m, off, smooth=11)
c, live = L.nmo(L.to_wavelet(d, 35.0, 25.0, point=True), t, kin, ("angle", 48.0))
S1 = L.stack3(c, live, off)
E = np.load(L.OUT + "/e2_stacks.npz")
S2 = {k: E["FS_trend_w25_a48__" + k] for k in ("near", "mid", "far")}
G2 = {k: E["G_trend_w25_a48__" + k] for k in ("near", "mid", "far")}
w = L.win(1.1, 1.9)
R = {"stacks_1m_vs_2m": {}}
for k in ("near", "mid", "far"):
    R["stacks_1m_vs_2m"][k] = dict(corr=float(np.corrcoef(S1[k][w], S2[k][w])[0, 1]), rms_ratio_1m_over_2m=float(L.rms(S1[k], w) / L.rms(S2[k], w)),
                                   rms_diff_rel=float(L.rms(S1[k] - S2[k], w) / L.rms(S2[k], w)))
    print(k, R["stacks_1m_vs_2m"][k])
pk = L.picks(G2["near"])
for nm, S in (("1m", S1), ("2m", S2)):
    rows = L.per_reflector(S, pk, "far", "near")
    R["zone_ref_" + nm] = dict(rms_fn=L.est_rms(S, "far", "near"), rms_mn=L.est_rms(S, "mid", "near"), peak_fn=L.zone_ref(rows, "pk")["ratio"], env_fn=L.zone_ref(rows, "env")["ratio"])
    print(nm, R["zone_ref_" + nm])
L.jdump(R, "e15_results.json")
print("saved")
