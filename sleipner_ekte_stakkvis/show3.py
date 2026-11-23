import json
R = json.load(open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis/results.json"))
S = R["oppsummering"]
allf = ["p07_dmE", "p10_dmE", "p07_dm", "p07_dmD", "p07_dmA", "p07_dmB"]
E = ("rms", "abs", "env", "coh", "env_topp_reflektor")
for v in ("p07", "p10"):
    for k in ("nea", "mid", "far"):
        n = S[v]["nivaa"][k]
        print(v, k, {e: n[e] for e in E}, "min", min(n[e] for e in E), "max", max(n[e] for e in E), "null", S[v]["nullvinduer"][k]["min"], S[v]["nullvinduer"][k]["max"])
for k in ("nea", "mid", "far"):
    vals = [S[v]["nivaa"][k][e] for v in allf for e in E]
    print("filt", k, "all est min/max", min(vals), max(vals), "rms", min(S[v]["nivaa"][k]["rms"] for v in allf), max(S[v]["nivaa"][k]["rms"] for v in allf), "null", min(S[v]["nullvinduer"][k]["min"] for v in allf), max(S[v]["nullvinduer"][k]["max"] for v in allf))
for name in ("mid_near", "far_near"):
    for v in ("p07", "p10"):
        d = S[v]["forhold"][name]
        print(v, name, {q: d[q] for q in d if q.startswith("vindu") or q.startswith("refl") or q.startswith("alle") or q.startswith("ampl")})
    print("filt", name, "alle", min(S[v]["forhold"][name]["alle_min"] for v in allf), max(S[v]["forhold"][name]["alle_max"] for v in allf))
    print("null unf", [(S[v]["nullvinduer"][name]["min"], S[v]["nullvinduer"][name]["max"], S[v]["nullvinduer"][name]["z_score"]) for v in ("p07", "p10")], "filt", min(S[v]["nullvinduer"][name]["min"] for v in allf), max(S[v]["nullvinduer"][name]["max"] for v in allf))
print("far_mid null filt", min(S[v]["nullvinduer"]["far_mid"]["min"] for v in allf), max(S[v]["nullvinduer"]["far_mid"]["max"] for v in allf))
print("endring", {k: (S["etter_filter_spenn"][k]["endring_rms_min"], S["etter_filter_spenn"][k]["endring_rms_max"]) for k in ("nea", "mid", "far")})
print("mid null", [S[v]["nullvinduer"]["mid"] for v in ("p07", "p10")])
print("trace", {v: {k: (S[v]["nivaa"][k]["trase_median"], S[v]["nivaa"][k]["trase_p16"], S[v]["nivaa"][k]["trase_p84"]) for k in ("nea", "mid", "far")} for v in ("p07", "p10")})
print("trace ratios", {v: {q: R[v]["vindu_forhold"]["kube"][q]["trace"]["rms"] for q in ("mid_near", "far_near")} for v in ("p07", "p10")})
print("bronn forhold", {v: {q: (S[v]["forhold"][q]["bronn_rms"], S[v]["forhold"][q]["bronn_coh"], S[v]["forhold"][q]["bronn_refl_rms"], S[v]["forhold"][q]["bronn_refl_ls"]) for q in ("mid_near", "far_near")} for v in ("p07", "p10", "p07_dmE", "p10_dmE")})
print("filt bronn rms", {q: (S["etter_filter_spenn"][q]["bronn_rms_min"], S["etter_filter_spenn"][q]["bronn_rms_max"]) for q in ("mid_near", "far_near")})
