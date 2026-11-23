"""Predictions table: what each explanation predicts for near alone, mid alone, far alone, far/near and mid/near
(zone against reference), for two sets of offset ranges. Reads results.json, writes prediksjoner.json and prints the
Markdown tables used in RESULTAT.md.

Run:  SLEIPNER_EXP_ROOT=... python -B s7_predictions.py
"""
import json
import os

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(ROOT, "results.json")))
SETS = [("modellens_offsetomraader", "main", "near 150-760 m, mid 760-1760 m, far 1760-3150 m (områdene den første målingen på de ekte dataene antok)"),
        ("offsetomraader_tilpasset_forsinkelsene", "reallag", "near 150-550 m, mid 575-1000 m, far 1200-1800 m (valgt slik at forsinkelsen på multippelkopiene blir som målt i de ekte stakkene: 110 / 112 / 119 ms)")]
EXPL = [
    ("ingen_forklaring", "G_trend", "Ingen av forklaringene: ingen multipler fra havflaten, Vs på trenden"),
    ("b_multipler_modellstyrke", "FS_trend", "(b) Multipler fra havflaten, så sterke som modellen gir dem. Vs på trenden"),
    ("b_multipler_skalert", "FSs_trend", "(b) Multipler skalert til kopistyrken i de ekte dataene (near 0,73, mid 0,62, far 0,55 av modellens). Vs på trenden"),
    ("a_myke_lag_vs07", "G_sb07", "(a) Vs 7 % under trenden i de myke lagene i sonen. Ingen multipler"),
    ("a_myke_lag_vs15", "G_sb15", "(a) Vs 15 % under trenden i de myke lagene i sonen. Ingen multipler"),
    ("ab_vs07_multipler_modellstyrke", "FS_sb07", "(a)+(b) Vs 7 % lavere i de myke lagene og multipler med modellens styrke"),
    ("ab_vs07_multipler_skalert", "FSs_sb07", "(a)+(b) Vs 7 % lavere i de myke lagene og skalerte multipler"),
    ("ab_vs15_multipler_modellstyrke", "FS_sb15", "(a)+(b) Vs 15 % lavere i de myke lagene og multipler med modellens styrke"),
    ("ab_vs15_multipler_skalert", "FSs_sb15", "(a)+(b) Vs 15 % lavere i de myke lagene og skalerte multipler"),
    ("a_ett_lag_vs15", "G_m15", "(a) Vs 15 % under trenden i ett homogent lag 1358-1390 m (det første forsøket). Ingen multipler"),
    ("kontroll_zoeppritz_myke_lag_vs07", "C_sb07", "Kontroll uten FD: bare primærer (Zoeppritz), Vs 7 % lavere i de myke lagene"),
    ("kontroll_zoeppritz_myke_lag_vs15", "C_sb15", "Kontroll uten FD: bare primærer (Zoeppritz), Vs 15 % lavere i de myke lagene"),
    ("kontroll_zoeppritz_myke_lag_vs20", "C_sb20", "Kontroll uten FD: bare primærer (Zoeppritz), Vs 20 % lavere i de myke lagene"),
]
r2 = lambda v: None if v is None else round(float(v), 3)
OUT = {
    "hva_dette_er": "Syntetisk 1D-modell av 15/9-13 (Sleipner), elastisk FD. Tallene er sone/referanse: verdien i sonen 1408-1492 ms delt på verdien i referansen (1100-1380 ms og 1560-1900 ms). 1,00 betyr at sonen er som resten.",
    "maal": {
        "near_alene, mid_alene, far_alene": "RMS-amplituden til den ene stakken i sonen delt på RMS-amplituden til samme stakk i referansen",
        "far_over_near.rms og mid_over_near.rms": "forholdet mellom RMS-amplitudene til to stakker i sonen, delt på det samme forholdet i referansen (= far_alene / near_alene)",
        "far_over_near.per_reflektor": "median over reflektorene i sonen delt på median over reflektorene i referansen, der far/near per reflektor er sum(far*mid)/sum(near*mid) i +-16 ms etter innretting mot mid. Reflektorene er de som plukkes på fullstakken uten multipler og med Vs på trenden (samme liste for alle forklaringer)",
        "far_over_near.per_reflektor_topp": "som over, men forholdet mellom toppamplitudene",
        "n_reflektorer": "antall reflektorer i sonen / i referansen. Få reflektorer i sonen gjør medianene grove; RMS-målet er derfor hovedmålet",
    },
    "multiplenes_styrke": "modellstyrke = som FD-modellen gir dem (havbunn-R 0,40). skalert = primærer + faktor x (multipler), faktor near 0,73, mid 0,62, far 0,55, som er kopistyrken målt i de ekte stakkene delt på modellens. De samme faktorene er brukt for begge offsetsettene",
    "myke_lag": "cellene i sonen der akustisk impedans ligger under glidende middel over 61 m: 1353-1362, 1366-1373, 1375-1399 og 1425-1435 m under havnivå (den uavhengige kontrollens definisjon)",
    "observert_i_de_ekte_stakkene": {"far_over_near": "0,48 (0,62-0,65 etter multippelfilter på stakkene)", "mid_over_near": "0,94-1,00 med den første målemetoden, ca. 0,6-0,75 med den andre", "near_alene, mid_alene, far_alene": "ikke målt da dette ble skrevet"},
    "slik_skilles_forklaringene": "Multipler løfter near i sonen (near_alene ca. 1,4-1,5 mot 1,1 uten) og senker far moderat. Lav Vs i de myke lagene lar near stå (ca. 1,0-1,1) og senker far kraftig. near_alene er derfor det tydeligste skillet. far_alene faller med begge (0,6-0,9) og skiller dårlig. mid_over_near skiller i de tilpassede offsetområdene (multipler ca. 0,8, lav Vs ca. 1,0).",
    "prediksjoner": {},
}
for key, vname, text in SETS:
    C = R["cases"].get(vname)
    if C is None:
        continue
    blk = {"offsetomraader": text, "forklaringer": {}}
    for ek, name, etext in EXPL:
        if name not in C:
            continue
        r, s = C[name]["rms"], C[name]["events_G"]["summary"]
        blk["forklaringer"][ek] = {
            "beskrivelse": etext, "datasett_i_results_json": "cases.%s.%s" % (vname, name),
            "near_alene": r2(r["near_zone_over_ref"]), "mid_alene": r2(r["mid_zone_over_ref"]), "far_alene": r2(r["far_zone_over_ref"]),
            "far_over_near": {"rms": r2(r["far_near_zone_over_ref"]), "per_reflektor": r2(s["fn"]["zone_over_ref"]), "per_reflektor_topp": r2(s["pk_fn"]["zone_over_ref"])},
            "mid_over_near": {"rms": r2(r["mid_near_zone_over_ref"]), "per_reflektor": r2(s["mn"]["zone_over_ref"]), "per_reflektor_topp": r2(s["pk_mn"]["zone_over_ref"])},
            "n_reflektorer": [s["n_zone"], s["n_ref"]],
        }
    tc = C.get("true_copies_of_strongest_primaries")
    if tc:
        blk["forsinkelse_paa_forste_kopi_i_modellen_ms"] = {k: tc[0][k]["first"][1] for k in ("nea", "mid", "far")}
        blk["forsinkelse_maalt_i_de_ekte_stakkene_ms"] = {"nea": 110, "mid": 112, "far": 119}
    OUT["prediksjoner"][key] = blk
json.dump(OUT, open(os.path.join(ROOT, "prediksjoner.json"), "w"), indent=1, ensure_ascii=False)
f = lambda v: "–" if v is None else ("%.2f" % v).replace(".", ",")
for key, blk in OUT["prediksjoner"].items():
    print("\n**%s**: %s" % (key, blk["offsetomraader"]))
    if "forsinkelse_paa_forste_kopi_i_modellen_ms" in blk:
        print("Forsinkelse på første kopi i modellen (near / mid / far):", blk["forsinkelse_paa_forste_kopi_i_modellen_ms"])
    print("\n| Forklaring | near alene | mid alene | far alene | far/near (RMS) | far/near (per reflektor) | mid/near (RMS) | mid/near (per reflektor) |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    for ek, e in blk["forklaringer"].items():
        print("| %s | %s | %s | %s | %s | %s | %s | %s |" % (e["beskrivelse"], f(e["near_alene"]), f(e["mid_alene"]), f(e["far_alene"]), f(e["far_over_near"]["rms"]),
                                                          f(e["far_over_near"]["per_reflektor"]), f(e["mid_over_near"]["rms"]), f(e["mid_over_near"]["per_reflektor"])))
