"""Point 2, continued: the multiples of the model with the lags and strengths of the copies in the real data.

Model copies (true copies under the reflector at 1156 ms): lag 111 / 118 / 131 ms (near / mid / far), strength -0.93 / -0.60 / -0.51.
Real data (parent session): 110 / 112 / 119 ms, -0.68 / -0.37 / -0.28.
Crude emulation: the multiples-only trace (FS - G) of each stack is shifted by the lag difference (-1 / -6 / -12 ms) and
scaled by the strength ratio. Also: the far multiples shifted from +6 to -24 ms in 2 ms steps (others fixed), and all
stacks shifted together by -10 .. +10 ms in 2 ms steps (tie / water-depth uncertainty). Writes p8_results.json.
"""
import json, os
import numpy as np
import slp, vz
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, "stacks.npz"))
KS = ("nea", "mid", "far")
G = {k: Z["G_trend_main__" + k].astype(float) for k in KS}
FS = {k: Z["FS_trend_main__" + k].astype(float) for k in KS}
M = {k: FS[k] - G[k] for k in KS}
ful = lambda S: sum(slp.W_FULL[k] * S[k] for k in KS)
evG = slp.pick_events(ful(G))
A = dict(nea=0.68 / 0.928, mid=0.37 / 0.596, far=0.28 / 0.509)
ONE = dict(nea=1.0, mid=1.0, far=1.0)
def meas(al, sh):
    S = {k: G[k] + al[k] * np.roll(M[k], int(sh[k])) for k in KS}
    S["ful"] = ful(S)
    r = vz.rms_measure(S); a = vz.event_measure(S, evG); b = vz.event_measure(S, None)
    return dict(rms_fn=r["fn"], evG_fn=a["fn"], evOwn_fn=b["fn"], pkG_fn=a["pk_fn"], rms_mn=r["mn"], evG_mn=a["mn"], evOwn_mn=b["mn"], near_zr=r["near_zone_over_ref"], far_zr=r["far_zone_over_ref"], zone_values=a["fn_zone_values"], n_own=[b["n_zone"], b["n_ref"]])
OUT = {}
real = dict(nea=-1, mid=-6, far=-12)
zero = dict(nea=0, mid=0, far=0)
OUT["cases"] = {"alpha1_model_lags": meas(ONE, zero), "alpha1_real_lags": meas(ONE, real), "setA_model_lags": meas(A, zero), "setA_real_lags": meas(A, real)}
OUT["far_shift_scan"] = {nm: [dict(far_shift_ms=s, **{q: v for q, v in meas(al, dict(nea=0, mid=0, far=s)).items() if q in ("rms_fn", "evG_fn", "evOwn_fn", "rms_mn", "evG_mn")}) for s in range(6, -25, -2)] for nm, al in (("alpha1", ONE), ("setA", A))}
OUT["common_shift_scan"] = {nm: [dict(shift_ms=s, **{q: v for q, v in meas(al, dict(nea=s, mid=s, far=s)).items() if q in ("rms_fn", "evG_fn", "evOwn_fn", "rms_mn", "evG_mn")}) for s in range(-10, 11, 2)] for nm, al in (("alpha1", ONE), ("setA", A))}
json.dump(OUT, open(os.path.join(HERE, "p8_results.json"), "w"), indent=1, default=float)
f = lambda v: "  nan" if v is None else "%5.2f" % v
print("case: far/near RMS evG evOwn pkG | mid/near RMS evG evOwn | near z/r far z/r")
for nm, r in OUT["cases"].items():
    print("  %-20s %s %s %s %s | %s %s %s | %s %s  %s n %s" % (nm, f(r["rms_fn"]), f(r["evG_fn"]), f(r["evOwn_fn"]), f(r["pkG_fn"]), f(r["rms_mn"]), f(r["evG_mn"]), f(r["evOwn_mn"]), f(r["near_zr"]), f(r["far_zr"]), r["zone_values"], r["n_own"]))
for key in ("far_shift_scan", "common_shift_scan"):
    for nm, rows in OUT[key].items():
        print(key, nm, "(shift: far/near RMS evG evOwn | mid/near RMS evG)")
        print("   " + "  ".join("%d: %s %s %s | %s %s" % (list(r.values())[0], f(r["rms_fn"]).strip(), f(r["evG_fn"]).strip(), f(r["evOwn_fn"]).strip(), f(r["rms_mn"]).strip(), f(r["evG_mn"]).strip()) for r in rows))
