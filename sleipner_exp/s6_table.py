"""Markdown tables from results.json (and validation_images.json)."""
import json
import os

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, "results.json")))
f2 = lambda v: "–" if v is None else ("%.2f" % v).replace(".", ",")
NAMES = [("G_trend", "Uten havflate, Vs på trend"), ("G_m07", "Uten havflate, Vs −7 %"), ("G_m15", "Uten havflate, Vs −15 %"),
         ("G_m20", "Uten havflate, Vs −20 %"), ("G_trend_sbl", "Uten havflate, havbunn-R 0,37"), ("G_trend_sbh", "Uten havflate, havbunn-R 0,43"),
         ("FS_trend", "Med havflate, Vs på trend, ufiltrert"), ("FS_m07", "Med havflate, Vs −7 %, ufiltrert"),
         ("FS_m15", "Med havflate, Vs −15 %, ufiltrert"), ("FS_m20", "Med havflate, Vs −20 %, ufiltrert"),
         ("FS_trend_sbl", "Med havflate, havbunn-R 0,37, ufiltrert"), ("FS_trend_sbh", "Med havflate, havbunn-R 0,43, ufiltrert"),
         ("FSdx1_trend", "Med havflate, Vs på trend, 1 m rutenett, ufiltrert"),
         ("FSsym_trend", "Med havflate, Vs på trend, eksakt speil, ufiltrert"),
         ("G_sb07", "Uten havflate, Vs −7 % i de myke lagene"), ("G_sb15", "Uten havflate, Vs −15 % i de myke lagene"),
         ("FS_sb07", "Med havflate, Vs −7 % i de myke lagene, ufiltrert"), ("FS_sb15", "Med havflate, Vs −15 % i de myke lagene, ufiltrert"),
         ("FSs_trend", "Med havflate, multipler skalert til ekte kopistyrke, Vs på trend"),
         ("FSs_sb07", "Skalerte multipler, Vs −7 % i de myke lagene"), ("FSs_sb15", "Skalerte multipler, Vs −15 % i de myke lagene"),
         ("C_sb07", "Bare primærer (Zoeppritz), Vs −7 % i de myke lagene"), ("C_sb15", "Bare primærer (Zoeppritz), Vs −15 % i de myke lagene"),
         ("C_sb20", "Bare primærer (Zoeppritz), Vs −20 % i de myke lagene"),
         ("FS_trend_filt", "Med havflate, Vs på trend, etter multippelfilter"), ("FS_m07_filt", "Med havflate, Vs −7 %, etter filter"),
         ("FS_m15_filt", "Med havflate, Vs −15 %, etter multippelfilter"), ("FS_m20_filt", "Med havflate, Vs −20 %, etter filter"),
         ("FS_trend_sbl_filt", "Med havflate, havbunn-R 0,37, etter filter"), ("FS_trend_sbh_filt", "Med havflate, havbunn-R 0,43, etter filter"),
         ("FSdx1_trend_filt", "Med havflate, 1 m rutenett, etter filter"), ("C_trend", "Bare primærer (Zoeppritz), Vs på trend"),
         ("C_m15", "Bare primærer (Zoeppritz), Vs −15 %"), ("C_m20", "Bare primærer (Zoeppritz), Vs −20 %"),
         ("C_m30", "Bare primærer (Zoeppritz), Vs −30 %"), ("C_m40", "Bare primærer (Zoeppritz), Vs −40 %"),
         ("CM_trend", "Stråleteori med vannlagsmultipler, Vs på trend"), ("CM_m15", "Stråleteori med vannlagsmultipler, Vs −15 %")]
for vname, C in R["cases"].items():
    print("\n### Variant %s: %s\n" % (vname, R["variants"][vname]))
    print("Hendelser (ms):", C["events_ms"])
    print("\n| Tilfelle | far/near sone | far/near referanse | **sone/referanse** | n sone / ref | mid/near sone/referanse | toppmål far/near sone/ref | RMS far/near sone/ref | RMS mid/near sone/ref |")
    print("|---|---:|---:|---:|---|---:|---:|---:|---:|")
    for key, lab in NAMES:
        if key not in C:
            continue
        ev = "events_FS" if (key.startswith("FS") or key.startswith("CM")) and not key.startswith("FSs") and "events_FS" in C[key] else "events_G"
        s = C[key][ev]["summary"]
        r = C[key]["rms"]
        print("| %s | %s | %s | **%s** | %d / %d | %s | %s | %s | %s |" % (lab, f2(s["fn"]["zone_median"]), f2(s["fn"]["ref_median"]), f2(s["fn"]["zone_over_ref"]),
              s["n_zone"], s["n_ref"], f2(s["mn"]["zone_over_ref"]), f2(s["pk_fn"]["zone_over_ref"]), f2(r["far_near_zone_over_ref"]), f2(r["mid_near_zone_over_ref"])))
    if "FS_trend" in C and "events_G" in C["FS_trend"]:
        print("\nMed havflate, målt på de sanne reflektorene (hendelsene fra kjøringen uten havflate):")
        for key in ("FS_trend", "FS_trend_filt", "FS_m15", "FS_m15_filt"):
            if key in C:
                s = C[key]["events_G"]["summary"]
                print("- %s: far/near sone/referanse %s, mid/near %s, toppmål %s (sone %s, n = %d)" % (key, f2(s["fn"]["zone_over_ref"]), f2(s["mn"]["zone_over_ref"]), f2(s["pk_fn"]["zone_over_ref"]), s["fn"]["zone_values"], s["n_zone"]))
    for k in ("filter_from_data", "true_copies_of_strongest_primaries", "multiples_over_primaries_rms", "multiples_over_primaries_rms_by_case", "grid_convergence_2m_vs_1m", "exact_mirror_vs_standard", "rms_far_over_near", "rms_mid_over_near"):
        if k in C:
            print("\n%s:\n%s" % (k, json.dumps(C[k], indent=1, default=float)))
print("\nchecks:", json.dumps(R["checks"], indent=1))
print("\nruns:", json.dumps(R["runs"], indent=1))
