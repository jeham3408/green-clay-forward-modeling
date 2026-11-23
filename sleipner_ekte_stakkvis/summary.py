#!/usr/bin/env python3
"""summary.py: derived numbers from results.json and profiles_*.npz (no cube access). Adds the key "oppsummering"."""
import json
import numpy as np

D = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis"
R = json.load(open(D + "/results.json"))
VERS = ["p07", "p10", "p07_dmE", "p10_dmE", "p07_dm", "p07_dmD", "p07_dmA", "p07_dmB"]
KS = ("nea", "mid", "far")
tms = 800.0 + 2.0 * np.arange(600)
ZM = (tms >= 1408) & (tms <= 1492)
RM = ((tms >= 1100) & (tms < 1380)) | ((tms > 1560) & (tms <= 1900))
WINS = [(1100, 1184), (1186, 1270), (1272, 1356), (1562, 1646), (1648, 1732), (1734, 1818), (1816, 1900)]
PAIRS = (("mid_near", "mid", "nea", "mn"), ("far_near", "far", "nea", "fn"))
r4 = lambda x: None if x is None or not np.isfinite(x) else round(float(x), 4)
S = {}
for v in VERS:
    r = R[v]
    P = np.load(f"{D}/profiles_{v}.npz")
    e = {k: P[f"e2t_{k}_kube"] for k in KS}
    ref = {k: e[k][RM].mean() for k in KS}
    lev = lambda k, m: float(np.sqrt(e[k][m].mean() / ref[k]))
    zone = {k: lev(k, ZM) for k in KS}
    null = [{k: lev(k, (tms >= a) & (tms <= b)) for k in KS} for a, b in WINS]
    s = dict(nivaa={}, forhold={}, nullvinduer=dict(vinduer_ms=WINS))
    for k in KS:
        n = r["nivaa"]
        nl = [w[k] for w in null]
        s["nivaa"][k] = dict(rms=n["kube"][k]["pooled"]["rms"], abs=n["kube"][k]["pooled"]["abs"], env=n["kube"][k]["pooled"]["env"],
                             coh=n["kube"][k]["pooled"]["coh"], env_topp_reflektor=r["reflektor"]["ful"]["kube"]["zr"][{"nea": "near_pk", "mid": "mid_pk", "far": "far_pk"}[k]],
                             trase_median=n["kube"][k]["trace"]["rms"]["median"], trase_p16=n["kube"][k]["trace"]["rms"]["p16"], trase_p84=n["kube"][k]["trace"]["rms"]["p84"],
                             delomr_median=n["delomraader"][k]["rms"]["median"], delomr_se=n["delomraader"][k]["rms"]["se"],
                             delomr_min=n["delomraader"][k]["rms"]["min"], delomr_max=n["delomraader"][k]["rms"]["max"],
                             bronn_rms=n["bronn"][k]["pooled"]["rms"], bronn_coh=n["bronn"][k]["pooled"]["coh"],
                             kontroll_rms_fra_profil=r4(zone[k]))
        s["nullvinduer"][k] = dict(sone=r4(zone[k]), vinduer=[r4(x) for x in nl], min=r4(min(nl)), max=r4(max(nl)))
    for name, k1, k2, pre in PAIRS:
        w = r["vindu_forhold"]
        rf = r["reflektor"]
        d = {"vindu_rms": w["kube"][name]["pooled"]["rms"], "vindu_abs": w["kube"][name]["pooled"]["abs"],
             "vindu_env": w["kube"][name]["pooled"]["env"], "vindu_coh": w["kube"][name]["pooled"]["coh"],
             "trase_median_rms": w["kube"][name]["trace"]["rms"]["median"]}
        for q in ("ls", "rms", "inv", "iv", "ls_med", "rms_med", "env_med", "env_pool"):
            d["refl_" + q] = rf["ful"]["kube"]["zr"][pre + "_" + q]
        amp = [d[q] for q in ("vindu_rms", "vindu_abs", "vindu_env", "vindu_coh", "refl_rms", "refl_rms_med", "refl_env_med", "refl_env_pool")]
        allv = [d[q] for q in d]
        picks = [rf[src]["kube"]["zr"][pre + "_" + q] for src in ("nea", "mid", "far", "ful") for q in ("ls", "rms", "inv", "iv", "ls_med", "rms_med", "env_med", "env_pool")]
        d.update(amplitudemaal_min=min(amp), amplitudemaal_max=max(amp), alle_min=min(allv + picks), alle_max=max(allv + picks),
                 delomr_rms_median=w["delomraader"][name]["rms"]["median"], delomr_rms_se=w["delomraader"][name]["rms"]["se"],
                 delomr_rms_min=w["delomraader"][name]["rms"]["min"], delomr_rms_max=w["delomraader"][name]["rms"]["max"],
                 delomr_coh_median=w["delomraader"][name]["coh"]["median"], delomr_coh_se=w["delomraader"][name]["coh"]["se"],
                 bronn_rms=w["bronn"][name]["pooled"]["rms"], bronn_coh=w["bronn"][name]["pooled"]["coh"],
                 bronn_refl_rms=rf["ful"]["bronn"]["zr"][pre + "_rms"], bronn_refl_ls=rf["ful"]["bronn"]["zr"][pre + "_ls"],
                 iv_nullforsinkelse=r["blokktraser_vindu_nullforsinkelse"]["zr"][pre + "_iv"])
        s["forhold"][name] = d
        z = zone[k1] / zone[k2]
        nl = [x[k1] / x[k2] for x in null]
        ln = np.log(nl)
        s["nullvinduer"][name] = dict(sone=r4(z), vinduer=[r4(x) for x in nl], min=r4(min(nl)), max=r4(max(nl)),
                                      sd_ln=r4(ln.std(ddof=1)), z_score=r4((np.log(z) - ln.mean()) / ln.std(ddof=1)))
    z = zone["far"] / zone["mid"]
    nl = [x["far"] / x["mid"] for x in null]
    ln = np.log(nl)
    s["nullvinduer"]["far_mid"] = dict(sone=r4(z), vinduer=[r4(x) for x in nl], min=r4(min(nl)), max=r4(max(nl)),
                                       sd_ln=r4(ln.std(ddof=1)), z_score=r4((np.log(z) - ln.mean()) / ln.std(ddof=1)))
    for q in ("vindu_rms", "vindu_coh", "refl_rms", "refl_ls", "refl_env_med"):
        s.setdefault("andel_av_ln_underskudd_paa_mid", {})[q] = r4(np.log(s["forhold"]["mid_near"][q]) / np.log(s["forhold"]["far_near"][q]))
    s["far_mid"] = {q: r4(s["forhold"]["far_near"][q] / s["forhold"]["mid_near"][q]) for q in ("vindu_rms", "vindu_coh", "refl_rms", "refl_ls", "refl_env_med")}
    c = r["korrelasjon_blokktraser_nullforsinkelse"]
    s["korrelasjon_nullforsinkelse"] = c
    s["kopi"] = {k: {q: r["kopi_1366_1466_av_1256_1356"][k][q] for q in ("lag_ms", "coef", "rho2")} for k in KS}
    z, rr = r["reflektor"]["ful"]["kube"]["zone"], r["reflektor"]["ful"]["kube"]["ref"]
    s["innretting"] = dict(mid_mot_near_ms=dict(sone_middel=z["lag_mid_mean_ms"], ref_middel=rr["lag_mid_mean_ms"], sone_median=z["lag_mid_ms"], ref_median=rr["lag_mid_ms"], sone_andel_paa_kanten=z["lag_mid_edge"], ref_andel_paa_kanten=rr["lag_mid_edge"]),
                           far_mot_near_ms=dict(sone_middel=z["lag_far_mean_ms"], ref_middel=rr["lag_far_mean_ms"], sone_median=z["lag_far_ms"], ref_median=rr["lag_far_ms"], sone_andel_paa_kanten=z["lag_far_edge"], ref_andel_paa_kanten=rr["lag_far_edge"]),
                           korr_innrettet_sone=dict(nm=z["c_nm"], nf=z["c_nf"], mf=z["c_mf"]), korr_innrettet_ref=dict(nm=rr["c_nm"], nf=rr["c_nf"], mf=rr["c_mf"]))
    S[v] = s
# 100 ms windows relative to the median of the clean reference windows (1100-1300 and 1600-1900 ms)
for v in VERS:
    g = R[v]["gain_100ms"]["kube"]
    refw = [3, 4, 8, 9, 10]
    S[v]["vindu100_rel"] = {name: {str(g["t_start_ms"][j]): r4(g[name][j] / np.median([g[name][i] for i in refw])) for j in range(3, 11)}
                            for name in ("nea", "mid", "far", "mid_near", "far_near")}
    S[v]["vindu100_rel"]["ref_spenn"] = {name: [r4(min(g[name][i] for i in refw) / np.median([g[name][i] for i in refw])), r4(max(g[name][i] for i in refw) / np.median([g[name][i] for i in refw]))] for name in ("mid_near", "far_near")}
# change produced by the multiple filters, stack by stack (filtered level / unfiltered level), and ranges over the filters
filt07 = ["p07_dmE", "p07_dm", "p07_dmD", "p07_dmA", "p07_dmB"]
allf = filt07 + ["p10_dmE"]
F = {}
for k in KS:
    F[k] = dict(rms_filtrert_min=min(S[v]["nivaa"][k]["rms"] for v in allf), rms_filtrert_max=max(S[v]["nivaa"][k]["rms"] for v in allf),
                coh_filtrert_min=min(S[v]["nivaa"][k]["coh"] for v in allf), coh_filtrert_max=max(S[v]["nivaa"][k]["coh"] for v in allf),
                endring_rms_min=r4(min(S[v]["nivaa"][k]["rms"] / S[v[:3]]["nivaa"][k]["rms"] for v in allf)),
                endring_rms_max=r4(max(S[v]["nivaa"][k]["rms"] / S[v[:3]]["nivaa"][k]["rms"] for v in allf)))
for name in ("mid_near", "far_near"):
    F[name] = {q + "_" + m: (min if m == "min" else max)(S[v]["forhold"][name][q] for v in allf)
               for q in ("vindu_rms", "vindu_abs", "vindu_env", "vindu_coh", "refl_ls", "refl_rms", "refl_inv", "refl_iv", "refl_env_med", "amplitudemaal_min", "amplitudemaal_max", "alle_min", "alle_max", "bronn_rms")
               for m in ("min", "max")}
    F[name]["nullvindu_z_min"] = min(S[v]["nullvinduer"][name]["z_score"] for v in allf)
    F[name]["nullvindu_z_max"] = max(S[v]["nullvinduer"][name]["z_score"] for v in allf)
F["andel_paa_mid_vindu_rms"] = [min(S[v]["andel_av_ln_underskudd_paa_mid"]["vindu_rms"] for v in allf), max(S[v]["andel_av_ln_underskudd_paa_mid"]["vindu_rms"] for v in allf)]
F["korr_nf_sone"] = [min(S[v]["korrelasjon_nullforsinkelse"]["nf_zone"] for v in allf), max(S[v]["korrelasjon_nullforsinkelse"]["nf_zone"] for v in allf)]
F["korr_nf_ref"] = [min(S[v]["korrelasjon_nullforsinkelse"]["nf_ref"] for v in allf), max(S[v]["korrelasjon_nullforsinkelse"]["nf_ref"] for v in allf)]
for name in ("mid_near", "far_near"):
    rest = [np.log(S[v]["forhold"][name]["vindu_rms"]) / np.log(S[v[:3]]["forhold"][name]["vindu_rms"]) for v in allf]
    F[name]["andel_av_ln_underskudd_som_staar_igjen_min"] = r4(min(rest))
    F[name]["andel_av_ln_underskudd_som_staar_igjen_max"] = r4(max(rest))
F["far_mid_vindu_rms"] = dict(ufiltrert=[S[v]["far_mid"]["vindu_rms"] for v in ("p07", "p10")], filtrert_min=min(S[v]["far_mid"]["vindu_rms"] for v in allf), filtrert_max=max(S[v]["far_mid"]["vindu_rms"] for v in allf),
                              null_z_ufiltrert=[S[v]["nullvinduer"]["far_mid"]["z_score"] for v in ("p07", "p10")],
                              null_z_filtrert_min=min(S[v]["nullvinduer"]["far_mid"]["z_score"] for v in allf), null_z_filtrert_max=max(S[v]["nullvinduer"]["far_mid"]["z_score"] for v in allf))
F["nullvinduer_spenn"] = {name: dict(ufiltrert_min=min(S[v]["nullvinduer"][name]["min"] for v in ("p07", "p10")), ufiltrert_max=max(S[v]["nullvinduer"][name]["max"] for v in ("p07", "p10")),
                                     filtrert_min=min(S[v]["nullvinduer"][name]["min"] for v in allf), filtrert_max=max(S[v]["nullvinduer"][name]["max"] for v in allf)) for name in ("mid_near", "far_near", "far_mid")}
F["vindu100_rel"] = {name: {w: [min(S[v]["vindu100_rel"][name][w] for v in allf), max(S[v]["vindu100_rel"][name][w] for v in allf)] for w in ("1300", "1400", "1500", "1600")} for name in ("mid_near", "far_near")}
S["etter_filter_spenn"] = F
# predictions of tonight's synthetic test, read from its file (RMS measures, both offset sets)
PF = json.load(open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_exp/prediksjoner.json"))
PR = {"kilde": "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_exp/prediksjoner.json (RMS-maalene)", "sett": {}}
meas = {v: [S[v]["nivaa"][k]["rms"] for k in KS] for v in ("p07", "p10")}
measf = {k: float(np.median([S[v]["nivaa"][k]["rms"] for v in allf])) for k in KS}
for sett, d in PF["prediksjoner"].items():
    o = {}
    for q, e in d["forklaringer"].items():
        p = [e["near_alene"], e["mid_alene"], e["far_alene"]]
        o[q] = dict(near=p[0], mid=p[1], far=p[2], mid_near=e["mid_over_near"]["rms"], far_near=e["far_over_near"]["rms"],
                    far_mid=r4(p[2] / p[1]),
                    avstand_ln_ufiltrert_p07=r4(np.sqrt(np.mean((np.log(meas["p07"]) - np.log(p)) ** 2))),
                    avstand_ln_ufiltrert_p10=r4(np.sqrt(np.mean((np.log(meas["p10"]) - np.log(p)) ** 2))),
                    avstand_ln_filtrert_median=r4(np.sqrt(np.mean((np.log([measf[k] for k in KS]) - np.log(p)) ** 2))))
    PR["sett"][sett] = dict(offsetomraader=d.get("offsetomraader"), forklaringer=o)
PR["maalt_filtrert_median_rms"] = {k: r4(measf[k]) for k in KS}
S["prediksjoner"] = PR
R["oppsummering"] = S
json.dump(R, open(D + "/results.json", "w"), indent=1)

f = lambda x: "  nan" if x is None else f"{x:5.2f}"
print("version   | near mid far (rms; coh; env-peak) | trace median near (p16-p84) | sub-area near med±se mid far")
for v in VERS:
    n = S[v]["nivaa"]
    print(f"{v:9s} | {f(n['nea']['rms'])} {f(n['mid']['rms'])} {f(n['far']['rms'])} ; {f(n['nea']['coh'])} {f(n['mid']['coh'])} {f(n['far']['coh'])} ; {f(n['nea']['env_topp_reflektor'])} {f(n['mid']['env_topp_reflektor'])} {f(n['far']['env_topp_reflektor'])} | " + " ".join(f"{f(n[k]['trase_median'])} ({f(n[k]['trase_p16'])}-{f(n[k]['trase_p84'])})" for k in KS) + " | " + " ".join(f"{f(n[k]['delomr_median'])}±{n[k]['delomr_se']:.3f} [{f(n[k]['delomr_min'])},{f(n[k]['delomr_max'])}]" for k in KS) + " | check " + " ".join(f(n[k]["kontroll_rms_fra_profil"]) for k in KS))
for name in ("mid_near", "far_near"):
    print(name)
    for v in VERS:
        d = S[v]["forhold"][name]
        nz = S[v]["nullvinduer"][name]
        print(f"{v:9s} | amp {f(d['amplitudemaal_min'])}-{f(d['amplitudemaal_max'])} all {f(d['alle_min'])}-{f(d['alle_max'])} | sub rms {f(d['delomr_rms_median'])}±{d['delomr_rms_se']:.3f} [{f(d['delomr_rms_min'])},{f(d['delomr_rms_max'])}] coh {f(d['delomr_coh_median'])}±{d['delomr_coh_se']:.3f} | bronn rms {f(d['bronn_rms'])} coh {f(d['bronn_coh'])} refl rms {f(d['bronn_refl_rms'])} ls {f(d['bronn_refl_ls'])} | iv0 {f(d['iv_nullforsinkelse'])} | null windows {nz['vinduer']} zone {nz['sone']} z {nz['z_score']}")
for v in VERS:
    print(v, "share on mid", S[v]["andel_av_ln_underskudd_paa_mid"], "far/mid", S[v]["far_mid"])
    print("   null near", S[v]["nullvinduer"]["nea"], "\n   null far", S[v]["nullvinduer"]["far"])
for q in ("mid_near", "far_near", "far_mid_vindu_rms", "nullvinduer_spenn", "andel_paa_mid_vindu_rms", "korr_nf_sone", "korr_nf_ref"):
    print(q, {k: v for k, v in F[q].items() if "staar" in k or "null" in k or "filtr" in k} if isinstance(F[q], dict) else F[q])
for v in VERS:
    print(v, "far/mid null", S[v]["nullvinduer"]["far_mid"])
for sett, d in PR["sett"].items():
    print("==", sett)
    for q, e in d["forklaringer"].items():
        print(f"  {q:36s} near {e['near']:.2f} mid {e['mid']:.2f} far {e['far']:.2f} mn {e['mid_near']:.2f} fn {e['far_near']:.2f} fm {e['far_mid']:.2f} | dist p07 {e['avstand_ln_ufiltrert_p07']:.3f} p10 {e['avstand_ln_ufiltrert_p10']:.3f} filt {e['avstand_ln_filtrert_median']:.3f}")
print(PR["maalt_filtrert_median_rms"])

for v in VERS:
    w = S[v]["vindu100_rel"]
    print(v, "100ms rel: mid/near", w["mid_near"], "| far/near", w["far_near"], "| ref span", w["ref_spenn"])
print(F["vindu100_rel"])
print("innretting", {v: S[v]["innretting"] for v in ("p07", "p10", "p07_dmE")})
print("gain p07/p10 raw", {v: {k: R[v]["gain_100ms"]["kube"][k] for k in ("nea", "mid", "far", "mid_near", "far_near")} for v in ("p07", "p10")})
print("sub 1400-1500", {v: R[v]["gain_100ms"]["vindu_1400_1500_mot_rene_referansevinduer_delomraader"] for v in ("p07", "p10")})
print("bronn", {v: {k: (S[v]["nivaa"][k]["bronn_rms"], S[v]["nivaa"][k]["bronn_coh"]) for k in KS} for v in VERS})
print("n traces", R["p07"]["nivaa"]["kube"]["n_traces"], R["p07"]["nivaa"]["bronn"]["n_traces"], "picks", R["p07"]["reflektor"]["ful"]["kube"]["zone"]["n"], R["p07"]["reflektor"]["ful"]["kube"]["ref"]["n"])
print("p07 vs p10", R["p07_mot_p10_korrelasjon_blokktraser_1100_1900_etter_forsinkelse_ms"])
print("kopi", {v: S[v]["kopi"] for v in VERS})
print("korr", {v: S[v]["korrelasjon_nullforsinkelse"] for v in VERS})
