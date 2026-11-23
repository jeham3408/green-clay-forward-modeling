import json
R = json.load(open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig_kontroll/results.json"))
def c(d):
    s = f"r={d['pearson']} rs={d['spearman']} n={d['n']}"
    if "pearson_ci95" in d: s += f" ci={d['pearson_ci95']} lat={d['lattice3']['mean_r']} [{d['lattice3']['min_r']},{d['lattice3']['max_r']}]"
    return s
print(R["time_origin_check"])
for p in ("p07", "p10", "p07_blocks300m", "p10_blocks300m"):
    r = R[p]
    print("=====", p, {k: r[k] for k in ("block_m", "n_valid", "area_km2", "flatten_shift_ms")})
    for k, v in r["global"].items(): print("  G", k, v)
    print("  copy", r["copy_global"]); print("  copyblk", r["copy_blocks"])
    for k, v in r["test2"].items():
        if "groups" in k: print("  T2", k, json.dumps(v))
        else: print("  T2", k, c(v))
    for w, v in r["placebo"].items():
        print("  PL", w, {k: v[k] for k in ("far_near_all", "n_negative", "coef_global_min", "coef_global_max")})
        for k, vv in v.items():
            if isinstance(vv, dict) and "pearson" in vv: print("     ", k, c(vv))
            elif "groups" in k: print("     ", k, json.dumps(vv))
    for w, v in r["independent_predictor"].items():
        for k, vv in v.items():
            if isinstance(vv, dict): print("  IND", w, k, c(vv))
    for w, v in r["test3"].items(): print("  T3", w, json.dumps(v))
    print("  T1", r["test1"])
    if "time_profile" in r:
        print("  TPS", r["time_profile_summary"])
        print("  TP", [(q["t"], q["far_nea"], q["mid_nea"]) for q in r["time_profile"]])
        print("  LV", [(q["t"], q["level_nea"], q["level_mid"], q["level_far"]) for q in r["time_profile"][::2]])
        for q in r["lag_profile"]: print("  LAG", q)
        print("  TR", r["area_mean_near_trace"])
