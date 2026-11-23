#!/usr/bin/env python3
"""Adds the few derived numbers quoted in NOTAT.md to results.json (ratios of numbers already there)."""
import json
P = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig/results.json"
R = json.load(open(P))
for p in ("p07", "p10"):
    g = R[p]["test2"]["copy_strength_groups_with_ci"]
    d = {f"level_{k}_strong_over_weak_tercile": round(g["strong_tercile"][f"level_{k}_refbelow"]["value"] / g["weak_tercile"][f"level_{k}_refbelow"]["value"], 4)
         for k in ("near", "mid", "far")}
    d["area_used_km2"] = round(R[p]["n_blocks_used"] * 0.2 * 0.225, 2)
    lp = {round(r["approx_time_ms"]): r["lag_far_after_near_ms"] for r in R[p]["lag_profile"]}
    d["lag_far_after_near_ms_by_time"] = lp
    ctl = R[p]["level_covariation_between_blocks"]
    d["sd_ln_far_level_zone"] = ctl["1400-1500 zone"]["sd_ln_level"]["far"]
    d["sd_ln_far_level_controls_min_max"] = [min(v["sd_ln_level"]["far"] for k, v in ctl.items() if "zone" not in k), max(v["sd_ln_level"]["far"] for k, v in ctl.items() if "zone" not in k)]
    R[p]["derived_for_notat"] = d
    print(p, json.dumps(d))
json.dump(R, open(P, "w"), indent=1)
