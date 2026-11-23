"""Point 4: numbers of RESULTAT.md against results.json / validation_images.json (the experiment's own files)."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "results.json")))
Vv = json.load(open(os.path.join(HERE, "validation_images.json")))
out = {}
def five(C, key):
    """the five measures of the report: two reflector lists x two estimators, plus RMS"""
    v = []
    for ev in ("events_G", "events_FS"):
        if ev in C[key]:
            s = C[key][ev]["summary"]
            v += [s["fn"]["zone_over_ref"], s["pk_fn"]["zone_over_ref"]]
    v.append(C[key]["rms"]["far_near_zone_over_ref"])
    return [round(x, 2) for x in v if x is not None]
def three(C, key):
    s = C[key]["events_G"]["summary"]
    return [round(s["fn"]["zone_over_ref"], 2), round(s["pk_fn"]["zone_over_ref"], 2), round(C[key]["rms"]["far_near_zone_over_ref"], 2)]
sens = {}
for var, C in R["cases"].items():
    row = {}
    for lab, g, f in (("trend", "G_trend", "FS_trend"), ("sbl", "G_trend_sbl", "FS_trend_sbl"), ("sbh", "G_trend_sbh", "FS_trend_sbh")):
        if g in C:
            a, b = three(C, g), five(C, f)
            c = five(C, f + "_filt") if f + "_filt" in C else []
            row[lab] = dict(G=[min(a), max(a)], FS=[min(b), max(b)], FS_filt=[min(c), max(c)] if c else None, FS_all=b, filt_all=c)
    # all G cases: range over Vs 0..20 %
    allg = []
    for k in C:
        if k.startswith("G_"):
            allg += three(C, k)
    row["G_all_cases_minmax"] = [min(allg), max(allg)]
    allfs = []
    for k in C:
        if k.startswith("FS_") and not k.endswith("_filt") and k.split("_")[1] in ("trend",):
            allfs += five(C, k)
    row["FS_trend_and_seabed_minmax"] = [min(allfs), max(allfs)]
    mn = []
    for k in ("FS_trend", "FS_trend_sbl", "FS_trend_sbh"):
        if k in C:
            for ev in ("events_G", "events_FS"):
                s = C[k][ev]["summary"]
                mn += [s["mn"]["zone_over_ref"], s["pk_mn"]["zone_over_ref"]]
            mn.append(C[k]["rms"]["mid_near_zone_over_ref"])
    row["FS_mid_near_minmax"] = [round(min(mn), 2), round(max(mn), 2)]
    sens[var] = row
out["sensitivity"] = sens
C = R["cases"]["main"]
out["G_main_by_vs"] = {k: three(C, k) for k in ("G_trend", "G_m07", "G_m15", "G_m20")}
out["G_zone_values_by_vs"] = {k: C[k]["events_G"]["summary"]["fn"]["zone_values"] for k in ("G_trend", "G_m07", "G_m15", "G_m20")}
out["C_main_by_vs"] = {k: three(C, k) for k in C if k.startswith("C_")}
out["FS_events_FS_zone_values"] = dict(values=C["FS_trend"]["events_FS"]["summary"]["fn"]["zone_values"], n_zone=C["FS_trend"]["events_FS"]["summary"]["n_zone"], n_bad=C["FS_trend"]["events_FS"]["summary"]["fn"]["n_bad"])
out["ref_spread_G"] = dict(min_max=C["G_trend"]["events_G"]["summary"]["fn"]["ref_min_max"], sd_log=C["G_trend"]["events_G"]["summary"]["fn"]["ref_sd_log"])
out["validation_images"] = Vv
txt = json.dumps(R) + json.dumps(Vv)
out["searched_strings_found"] = {s: (s in txt) for s in ("0.995", "0.93", "courant", "Courant", "cfl", "0.68", "1.10", "0.37", "dx1", "grid_convergence")}
out["run_info_example"] = R["runs"].get("FS_trend")
out["has_dx1"] = "FSdx1_trend" in C
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, "p4_audit.json"), "w"), indent=1)
