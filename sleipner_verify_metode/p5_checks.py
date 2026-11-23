"""Numbers of RESULTAT.md that are not in results.json / tables.md, recomputed from the experiment's stacks.npz."""
import json, os
import numpy as np
import slp, vz
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, "stacks.npz"))
t = Z["t"]
g = lambda cfg, case, k, var="main": Z["%s_%s_%s__%s" % (cfg, case, var, k)].astype(float)
rms = lambda x, w: float(np.sqrt(np.mean(x[w] ** 2)))
out = {}
w = (t >= 1.380) & (t <= 1.500)
for case in ("m07", "m15", "m20"):
    row = {}
    for k in ("nea", "mid", "far"):
        dfd = g("G", case, k) - g("G", "trend", k)
        dz = g("C", case, k) - g("C", "trend", k)
        row[k] = dict(corr_fd_vs_zoeppritz=float(np.corrcoef(dfd[w], dz[w])[0, 1]), strength_fd_over_zoeppritz=rms(dfd, w) / rms(dz, w),
                      change_over_level_fd=rms(dfd, w) / rms(g("G", "trend", k), w), level_ratio_fd=rms(g("G", case, k), w) / rms(g("G", "trend", k), w))
    out["vs_effect_1380_1500_" + case] = row
# zone and reference level of far and near, without sea surface, by Vs case (RMS)
zm, rm = vz.win_masks()
out["far_rms_zone_over_ref_by_case"] = {c: rms(g("G", c, "far"), zm) / rms(g("G", c, "far"), rm) for c in ("trend", "m07", "m15", "m20")}
# per-reflector far/near, without sea surface, by Vs case, at the true reflectors
G = {k: g("G", "trend", k) for k in ("nea", "mid", "far")}
ev = slp.pick_events(sum(slp.W_FULL[k] * G[k] for k in G))
tab = {}
for c in ("trend", "m07", "m15", "m20"):
    S = {k: g("G", c, k) for k in ("nea", "mid", "far")}
    tab[c] = {int(round(1000 * i * slp.DT)): round(slp.measure_event(S, i)["fn"], 3) for i in ev if 1.38 <= i * slp.DT <= 1.56}
out["far_near_per_reflector_G_by_case"] = tab
# FD without sea surface against the Zoeppritz synthetic: the same measure gives different zone/reference
out["zone_over_ref_G_vs_C_trend"] = dict(G=vz.event_measure({**G, "ful": sum(slp.W_FULL[k] * G[k] for k in G)}, ev)["fn"],
                                         C=vz.event_measure({**{k: g("C", "trend", k) for k in G}, "ful": sum(slp.W_FULL[k] * g("C", "trend", k) for k in G)}, ev)["fn"])
print(json.dumps(out, indent=1, default=float))
json.dump(out, open(os.path.join(HERE, "p5_checks.json"), "w"), indent=1, default=float)
