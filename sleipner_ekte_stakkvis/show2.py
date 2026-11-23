import json, sys
R = json.load(open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis/results.json"))
vers = [k for k in R if not k.startswith("_") and "korrelasjon" not in k and k != "oppsummering"]
f = lambda x: "  nan" if x is None else f"{x:5.2f}"
print("version   | level z/r pooled rms: near mid far | coh: near mid far | sub-area rms med±se near, mid, far")
for v in vers:
    n = R[v]["nivaa"]; p = n["kube"]; d = n["delomraader"]
    print(f"{v:9s} | {f(p['nea']['pooled']['rms'])} {f(p['mid']['pooled']['rms'])} {f(p['far']['pooled']['rms'])} | {f(p['nea']['pooled']['coh'])} {f(p['mid']['pooled']['coh'])} {f(p['far']['pooled']['coh'])} | " + " ".join(f"{f(d[k]['rms']['median'])}±{d[k]['rms']['se']:.3f}" for k in ("nea","mid","far")) + " | bronn rms " + " ".join(f(n['bronn'][k]['pooled']['rms']) for k in ("nea","mid","far")))
print("\nversion   | window mid/near: rms abs env coh | far/near: rms abs env coh | trace-median rms mn fn | sub rms mn med±se, fn med±se | bronn rms mn fn")
for v in vers:
    w = R[v]["vindu_forhold"]; k = w["kube"]; d = w["delomraader"]; b = w["bronn"]
    print(f"{v:9s} | " + " ".join(f(k['mid_near']['pooled'][e]) for e in ("rms","abs","env","coh")) + " | " + " ".join(f(k['far_near']['pooled'][e]) for e in ("rms","abs","env","coh")) + f" | {f(k['mid_near']['trace']['rms']['median'])} {f(k['far_near']['trace']['rms']['median'])} | {f(d['mid_near']['rms']['median'])}±{d['mid_near']['rms']['se']:.3f} {f(d['far_near']['rms']['median'])}±{d['far_near']['rms']['se']:.3f} | {f(b['mid_near']['pooled']['rms'])} {f(b['far_near']['pooled']['rms'])}")
print("\nreflector estimators zone/ref, kube, pick=ful (and range over pick sources)")
keys = ("ls","rms","inv","iv","ls_med","rms_med","env_med","env_pool")
for v in vers:
    x = R[v]["reflektor"]
    for pre in ("mn","fn"):
        s = []
        for e in keys:
            vals = [x[src]["kube"]["zr"][pre+"_"+e] for src in ("nea","mid","far","ful")]
            s.append(f"{e} {f(x['ful']['kube']['zr'][pre+'_'+e])}[{f(min(vals))},{f(max(vals))}]")
        print(f"{v:9s} {pre}: " + " ".join(s))
    z = x["ful"]["kube"]
    print(f"          env-peak z/r near {f(z['zr']['near_pk'])} mid {f(z['zr']['mid_pk'])} far {f(z['zr']['far_pk'])} | lags mid,far zone {z['zone']['lag_mid_mean_ms']:+.1f},{z['zone']['lag_far_mean_ms']:+.1f} ref {z['ref']['lag_mid_mean_ms']:+.1f},{z['ref']['lag_far_mean_ms']:+.1f} med zone {z['zone']['lag_mid_ms']:+.0f},{z['zone']['lag_far_ms']:+.0f} ref {z['ref']['lag_mid_ms']:+.0f},{z['ref']['lag_far_ms']:+.0f} edge mid z/r {z['zone']['lag_mid_edge']:.2f}/{z['ref']['lag_mid_edge']:.2f} far {z['zone']['lag_far_edge']:.2f}/{z['ref']['lag_far_edge']:.2f} | corr zone nm nf mf {f(z['zone']['c_nm'])} {f(z['zone']['c_nf'])} {f(z['zone']['c_mf'])} ref {f(z['ref']['c_nm'])} {f(z['ref']['c_nf'])} {f(z['ref']['c_mf'])}")
    s = x["ful"]["delomraader"]
    print("          sub: " + " ".join(f"{k} {f(s[k]['median'])}±{s[k]['se']:.2f}" for k in ("mn_ls","mn_rms","mn_iv","mn_env_med","fn_ls","fn_rms","fn_iv","fn_env_med")))
print("\nzero-lag block traces: corr and IV")
for v in vers:
    c = R[v]["korrelasjon_blokktraser_nullforsinkelse"]; iv = R[v]["blokktraser_vindu_nullforsinkelse"]
    print(f"{v:9s} corr zone nm {c['nm_zone']} nf {c['nf_zone']} mf {c['mf_zone']} | ref nm {c['nm_ref']} nf {c['nf_ref']} mf {c['mf_ref']} | zr: " + " ".join(f"{q} {f(iv['zr'][q])}" for q in iv["zr"]) + f" | zone mn_iv {iv['zone']['mn_iv']} ref {iv['ref']['mn_iv']}")
print("\ncopy")
for v in vers:
    c = R[v]["kopi_1366_1466_av_1256_1356"]
    print(f"{v:9s} " + " | ".join(f"{k}: lag {c[k]['lag_ms']:.0f} coef {c[k]['coef']:+.2f} rho2 {c[k]['rho2']:.2f} sub med {c[k]['coef_delomraader']['median']:+.2f} sd {c[k]['coef_delomraader']['sd']:.2f}" for k in ("nea","mid","far")))
print("\ngain 100 ms windows (kube)")
for v in vers[:2]:
    g = R[v]["gain_100ms"]["kube"]
    print(v, "t   :", " ".join(f"{t:5d}" for t in g["t_start_ms"]))
    for k in ("nea", "mid", "far", "mid_near", "far_near"):
        print(f"  {k:9s}:", " ".join(f"{x:5.3f}" for x in g[k]))
    print("  1400-1500 vs clean ref:", {k: (x["median"], x["se"]) for k, x in R[v]["gain_100ms"]["vindu_1400_1500_mot_rene_referansevinduer_delomraader"].items()})
for k in R:
    if "korrelasjon" in k and k.startswith("p07_mot"): print(k, R[k])
