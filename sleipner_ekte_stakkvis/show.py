import json, sys
R = json.load(open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis/results.json"))
vers = sys.argv[1:] or [k for k in R if not k.startswith("_") and "korrelasjon" not in k]
f = lambda x: "  nan" if x is None else f"{x:5.2f}"
for v in vers:
    r = R[v]
    print(f"\n===== {v}  live {r['live_fraction']}")
    for an in ("kube", "bronn"):
        n = r["nivaa"][an]
        print(f" [{an}] n={n['n_traces']}  level zone/ref  pooled rms/abs/env/coh | trace median rms (p16-p84), abs, env")
        for k in ("nea", "mid", "far"):
            p, t = n[k]["pooled"], n[k]["trace"]
            print(f"   {k}: {f(p['rms'])} {f(p['abs'])} {f(p['env'])} {f(p['coh'])} | {f(t['rms']['median'])} ({f(t['rms']['p16'])}-{f(t['rms']['p84'])}) {f(t['abs']['median'])} {f(t['env']['median'])}")
        for q in ("mid_near", "far_near", "far_mid"):
            p, t = r["vindu_forhold"][an][q]["pooled"], r["vindu_forhold"][an][q]["trace"]
            print(f"   {q}: pooled {f(p['rms'])} {f(p['abs'])} {f(p['env'])} {f(p['coh'])} | trace {f(t['rms']['median'])} ({f(t['rms']['p16'])}-{f(t['rms']['p84'])}) {f(t['abs']['median'])} {f(t['env']['median'])}")
    d = r["nivaa"]["delomraader"]
    for k in ("nea", "mid", "far"):
        print(f"   sub {k}: " + "  ".join(f"{e} med {f(d[k][e]['median'])} sd {f(d[k][e]['sd'])} se {d[k][e]['se']:.3f} [{f(d[k][e]['min'])},{f(d[k][e]['max'])}]" for e in ("rms", "env", "coh")))
    d = r["vindu_forhold"]["delomraader"]
    for q in ("mid_near", "far_near"):
        print(f"   sub {q}: " + "  ".join(f"{e} med {f(d[q][e]['median'])} sd {f(d[q][e]['sd'])} se {d[q][e]['se']:.3f} [{f(d[q][e]['min'])},{f(d[q][e]['max'])}]" for e in ("rms", "abs", "env", "coh")))
    g = r["gain_100ms"]["kube"]
    print("  gain t   :", " ".join(f"{t:5d}" for t in g["t_start_ms"]))
    for k in ("nea", "mid", "far", "mid_near", "far_near"):
        print(f"  {k:9s}:", " ".join(f"{x:5.3f}" for x in g[k]))
    print("  1400-1500 vs clean ref windows, sub-areas:", {k: (v["median"], v["sd"], v["se"]) for k, v in r["gain_100ms"]["vindu_1400_1500_mot_rene_referansevinduer_delomraader"].items()})
    for src in ("nea", "mid", "far", "ful"):
        x = r["reflektor"][src]
        for an in ("kube", "bronn"):
            z = x[an]
            print(f"  refl pick={src} [{an}] n zone/ref {z['zone']['n']}/{z['ref']['n']}  lag mid,far ms zone {z['zone']['lag_mid_mean_ms']:+.1f},{z['zone']['lag_far_mean_ms']:+.1f} ref {z['ref']['lag_mid_mean_ms']:+.1f},{z['ref']['lag_far_mean_ms']:+.1f} edge far z/r {z['zone']['lag_far_edge']:.2f}/{z['ref']['lag_far_edge']:.2f} corr z nm,nf,mf {f(z['zone']['c_nm'])} {f(z['zone']['c_nf'])} {f(z['zone']['c_mf'])} ref {f(z['ref']['c_nm'])} {f(z['ref']['c_nf'])} {f(z['ref']['c_mf'])}")
            print("      mid/near z/r: " + " ".join(f"{k[3:]} {f(z['zr'][k])}" for k in z["zr"] if k.startswith("mn_")))
            print("      far/near z/r: " + " ".join(f"{k[3:]} {f(z['zr'][k])}" for k in z["zr"] if k.startswith("fn_")))
            print(f"      env peak z/r near {f(z['zr']['near_pk'])} mid {f(z['zr']['mid_pk'])} far {f(z['zr']['far_pk'])}   raw zone mid/near: ls {f(z['zone']['mn_ls'])} rms {f(z['zone']['mn_rms'])} iv {f(z['zone']['mn_iv'])} | ref: ls {f(z['ref']['mn_ls'])} rms {f(z['ref']['mn_rms'])} iv {f(z['ref']['mn_iv'])}")
        s = x["delomraader"]
        print("      sub mid/near: " + " ".join(f"{k[3:]} {f(s[k]['median'])}±{s[k]['se']:.2f}" for k in s if k.startswith("mn_")))
        print("      sub far/near: " + " ".join(f"{k[3:]} {f(s[k]['median'])}±{s[k]['se']:.2f}" for k in s if k.startswith("fn_")))
    print("  corr block traces zero lag:", r["korrelasjon_blokktraser_nullforsinkelse"])
    for k in ("nea", "mid", "far"):
        c = r["kopi_1366_1466_av_1256_1356"][k]
        print(f"  copy {k}: lag {c['lag_ms']} coef {c['coef']} rho {c['rho']} rho2 {c['rho2']} sub coef med {c['coef_delomraader']['median']} sd {c['coef_delomraader']['sd']}")
for k in R:
    if "korrelasjon" in k: print(k, R[k])
