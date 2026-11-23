#!/usr/bin/env python3
"""post3.py (kontroll): compact summary from results.json (no cube access). Adds the key "oppsummering"."""
import json, os
import numpy as np
OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
R = json.load(open(os.path.join(OUT, "results.json")))
DM = ("p07_dmE", "p10_dmE", "p07_dm", "p07_dmA", "p07_dmB", "p07_dmD")
Z = "sone_1408_1492"
g = lambda v, key, e="rms", w=Z: R[v]["nivaa_og_forhold"][w][key]["hele"][e]
rng = lambda xs: [round(float(min(xs)), 3), round(float(max(xs)), 3)]
S = {}
EST = ("rms", "abs", "nabo_xl", "nabo_il", "blokktrase_rms", "blokktrase_env_median")
for key in ("nea", "mid", "far", "mid/nea", "far/nea", "far/mid"):
    o = {}
    for v in ("p07", "p10"):
        o[v] = dict(rms=g(v, key), spenn_seks_maal=rng([g(v, key, e) for e in EST]), ren_referanse=g(v, key, "rms_ren_ref"),
                    delomr_median=R[v]["nivaa_og_forhold"][Z][key]["delomraader_rms"]["median"],
                    delomr_se=R[v]["nivaa_og_forhold"][Z][key]["delomraader_rms"]["se"],
                    delomr_min=R[v]["nivaa_og_forhold"][Z][key]["delomraader_rms"]["min"],
                    delomr_max=R[v]["nivaa_og_forhold"][Z][key]["delomraader_rms"]["max"])
    o["etter_filter_rms"] = rng([g(v, key) for v in DM])
    o["etter_filter_alle_maal"] = rng([g(v, key, e) for v in DM for e in EST])
    S[key] = o
# neutral-pick envelope peaks
S["envelopetopper"] = {v: {k: R[v]["reflektor_envelopetopper_noytrale_plukk"][k] for k in
                           ("nea", "mid", "far", "mid/nea_median_av_forhold", "far/nea_median_av_forhold", "far/mid_median_av_forhold")}
                       for v in ("p07", "p10")}
S["envelopetopper"]["etter_filter"] = {k: rng([R[v]["reflektor_envelopetopper_noytrale_plukk"][k] for v in DM]) for k in
                                       ("nea", "mid", "far", "mid/nea_median_av_forhold", "far/nea_median_av_forhold")}
# decomposition of ln(far/near zone/ref): part removed by each filter, far/mid after filter, mid/near after filter
dec = {}
for v in DM:
    base = "p10" if v.startswith("p10") else "p07"
    tot = np.log(g(base, "far/nea"))
    aft = np.log(g(v, "far/nea"))
    dec[v] = dict(ln_far_near_raa=tot, andel_fjernet_av_filter=(tot - aft) / tot, andel_far_mid_etter_filter=np.log(g(v, "far/mid")) / tot,
                  andel_mid_near_etter_filter=np.log(g(v, "mid/nea")) / tot)
S["dekomponering_ln_far_near"] = dec
S["dekomponering_spenn"] = {k: rng([dec[v][k] for v in DM]) for k in ("andel_fjernet_av_filter", "andel_far_mid_etter_filter", "andel_mid_near_etter_filter")}
# copy
cp = {}
for v in ("p07", "p10"):
    cp[v] = {k: dict(lag_ms=R[v]["kopi"][k]["lag_ms"], koef_traser=R[v]["kopi"][k]["koef"], andel_traser=R[v]["kopi"][k]["andel_energi"],
                     koef_blokk=R[v]["kopi"][k]["blokktrase_koef"], andel_blokk=R[v]["kopi"][k]["blokktrase_andel"],
                     andel_sone_traser=R[v]["kopi"][k]["andel_sone_1408_1492"],
                     andel_leirvindu_traser=R["ekstra"][v]["kopiandel_maalvindu"][k]["leirlaget_1412_1446"]["andel"],
                     referanse_andel_spenn=rng([-x for x in R["ekstra"][v]["kopiandel_maalvindu"][k]["referanse_signert_andel"]]),
                     referanse_andel_median=float(np.median([-x for x in R["ekstra"][v]["kopiandel_maalvindu"][k]["referanse_signert_andel"]])))
             for k in ("nea", "mid", "far")}
S["kopi"] = cp
# coefficient after filter at the lag the copy had before filtering, and the most negative one
af = {}
for v in DM:
    base = "p10" if v.startswith("p10") else "p07"
    o = {}
    for k in ("nea", "mid", "far"):
        lags = R[v]["kopi"][k]["lag_ms_alle"]
        j = lags.index(R[base]["kopi"][k]["lag_ms"])
        o[k] = dict(ved_opprinnelig_lag=R[v]["kopi"][k]["koef_alle_lag"][j], mest_negative=R[v]["kopi"][k]["koef"], lag_mest_negative=R[v]["kopi"][k]["lag_ms"])
    af[v] = o
S["kopi_etter_filter"] = af
S["kopi_etter_filter_near_ved_opprinnelig_lag_spenn"] = rng([af[v]["nea"]["ved_opprinnelig_lag"] for v in DM])
# alignment
S["innretting"] = {v: {p: {w: (R[v]["innretting_traser"][p][w]["lag_ms"], R[v]["innretting_traser"][p][w]["korr_best"], R[v]["innretting_traser"][p][w]["korr_null"])
                           for w in ("sone", "ref")} for p in ("nea-mid", "nea-far", "mid-far")} for v in ("p07", "p10", "p07_dmE", "p10_dmE")}
S["blokk_korr"] = {v: {w: R[v]["blokktraser"][w]["korr_nf"] for w in ("sone", "ref")} for v in ("p07", "p10") + DM}
S["instrument_mid_near_sone_over_ref"] = {v: dict(traser=R[v]["minste_kvadrater_traser"]["instrument_null_lag"]["sone_over_ref"]["mid_near"],
                                                  blokk=R[v]["blokktraser"]["sone_over_ref"]["iv_mid_near"]) for v in ("p07", "p10", "p07_dmE", "p10_dmE")}
S["minste_kvadrater_far_near_sone_over_ref"] = {v: dict(null_lag=R[v]["minste_kvadrater_traser"]["far/nea_null"]["sone_over_ref"]["ls"],
                                                        beste_lag=R[v]["minste_kvadrater_traser"]["far/nea_best"]["sone_over_ref"]["ls"]) for v in ("p07", "p10")}
# location relative to the soft layers
E = R["ekstra"]
S["myke_lag_mot_naboene"] = {v: E[v]["myke_lag_mot_naboene_over_og_under"] for v in ("p07", "p10") + DM}
S["myke_lag_mot_naboene_delomr"] = {v: E[v]["myke_lag_mot_naboene_delomraader_rms"] for v in ("p07", "p10") + DM}
S["far_mid_vinduer"] = {v: {w: (E[v]["vinduer_mot_myke_lag"][w]["far/mid"]["rms"], E[v]["vinduer_mot_myke_lag"][w]["far/mid"]["nabo"]) for w in E[v]["vinduer_mot_myke_lag"]}
                        for v in ("p07", "p10", "p07_dmE", "p10_dmE")}
S["filter_near_sone"] = rng([E["filter"][v]["nea"]["sone"] for v in DM])
S["filter_near_ref"] = rng([E["filter"][v]["nea"]["ref"] for v in DM])
S["filter_near_andel_fra_ref"] = rng([E["filter"][v]["nea"]["andel_av_endring_fra_referansen"] for v in DM])
S["filter_mid_sone"] = rng([E["filter"][v]["mid"]["sone"] for v in DM])
S["filter_far_sone"] = rng([E["filter"][v]["far"]["sone"] for v in DM])
S["filter_ln_avvik_igjen_mid_near"] = rng([E["filter"][v]["mid/nea_andel_av_ln_avvik_igjen"] for v in DM])
S["filter_ln_avvik_igjen_far_near"] = rng([E["filter"][v]["far/nea_andel_av_ln_avvik_igjen"] for v in DM])
S["near_uten_kopi"] = {v: dict(bare_sone=R[v]["sone_uten_kopiandel"]["nea"], sone_og_referanse=E[v]["near_sone_uten_kopi_i_sone_og_referanse"]["near_korrigert"],
                               mid=R[v]["sone_uten_kopiandel"]["mid"], far=R[v]["sone_uten_kopiandel"]["far"],
                               mid_near=R[v]["sone_uten_kopiandel"]["mid/nea"], far_near=R[v]["sone_uten_kopiandel"]["far/nea"]) for v in ("p07", "p10")}
S["delomr_kopi"] = {v: E[v]["delomraader_kopi"] for v in ("p07", "p10")}
S["glidende"] = {v: {q: {k: R[v]["glidende_vindu_84ms"][q][k] for k in ("min_verdi", "min_vindu_start_ms", "sone_verdi", "z_sone", "uavhengige_vinduer_verdier")}
                     for q in ("mid/nea", "far/nea", "far/mid")} for v in ("p07", "p10", "p07_dmE", "p10_dmE")}
S["kopilast_spearman"] = {v: {k: R[v]["kopilast_mot_forhold_20ms"][k] for k in R[v]["kopilast_mot_forhold_20ms"] if k.startswith("spearman")} for v in ("p07", "p10", "p07_dmE", "p10_dmE")}
S["per_trase"] = {v: E[v]["per_trase_1400_1500"] for v in ("p07", "p10")}
S["bronn"] = {v: E[v]["rundt_bronnen_1400_1500"] for v in ("p07", "p10")}
S["gain"] = {v: {k: R[v]["gain_100ms"][k] for k in ("nea", "mid", "far", "mid/nea", "far/nea", "far/mid")} for v in ("p07", "p10")}
S["sjekk0"] = json.load(open(os.path.join(OUT, "check0.json")))
R["oppsummering"] = S
json.dump(R, open(os.path.join(OUT, "results.json"), "w"), indent=1, ensure_ascii=False)
def show(o, ind=0):
    for k, v in o.items():
        if isinstance(v, dict) and any(isinstance(x, dict) for x in v.values()):
            print(" " * ind + str(k) + ":")
            show(v, ind + 2)
        else:
            print(" " * ind + f"{k}: {json.dumps(v)}")
show({k: v for k, v in S.items() if k not in ("sjekk0",)})
