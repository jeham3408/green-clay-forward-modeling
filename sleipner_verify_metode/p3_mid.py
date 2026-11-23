"""Point 3 and extras for point 2.

(1) The copy measure (the experiment's copy_coef, as used for the filter and for the comparison with the real data)
    applied to the stacks WITHOUT sea surface: what "copies" does it report when there are no multiples at all?
(2) Any per-stack strength of the multiples (alpha_near, alpha_mid, alpha_far in 0..1.2) that gives BOTH the observed
    far/near deficit (zone/reference 0.40-0.56) and a normal mid/near (>= 0.90)?
(3) The multiples shifted in time against the primaries (other reflectors above the zone / tie error): range of the deficit.
Writes p3_results.json.
"""
import json
import os

import numpy as np

import slp
import vz

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, "stacks.npz"))
KS = ("nea", "mid", "far")
get = lambda cfg, case, var: {k: Z["%s_%s_%s__%s" % (cfg, case, var, k)].astype(float) for k in KS}
G, FS = get("G", "trend", "main"), get("FS", "trend", "main")
M = {k: FS[k] - G[k] for k in KS}
ful = lambda S: sum(slp.W_FULL[k] * S[k] for k in KS)
evG = slp.pick_events(ful(G))
evFS = slp.pick_events(ful(FS))
OUT = {}

# ---- (1) copy measure without any multiples
rows = {}
for nm, S, ev in (("G_no_multiples", G, evG), ("FS_with_multiples", FS, evFS)):
    par_all = [i for i in ev if 1.15 <= i * slp.DT <= 1.38]
    r = []
    for i in par_all:
        row = dict(parent_ms=int(round(1000 * i * slp.DT)))
        for k in KS:
            c1 = slp.copy_coef(S[k], i, 0.080, 0.150, sign=-1)
            c2 = slp.copy_coef(S[k], i, 0.190, 0.300, sign=+1)
            row[k] = dict(first=[round(c1[0], 2), round(1000 * c1[1])], second=[round(c2[0], 2), round(1000 * c2[1])])
        r.append(row)
    rows[nm] = r
OUT["copy_measure_by_parent"] = rows
summ = {}
for nm, r in rows.items():
    summ[nm] = {k: dict(first_median=float(np.median([q[k]["first"][0] for q in r])), first_min_max=[float(min(q[k]["first"][0] for q in r)), float(max(q[k]["first"][0] for q in r))],
                        first_lag_ms=[q[k]["first"][1] for q in r], second_median=float(np.median([q[k]["second"][0] for q in r]))) for k in KS}
OUT["copy_measure_summary"] = summ


def meas(al, shift_ms=0):
    S = {}
    for k in KS:
        mk = M[k] if shift_ms == 0 else np.roll(M[k], int(shift_ms))
        S[k] = G[k] + al[k] * mk
    S["ful"] = ful(S)
    r = vz.rms_measure(S)
    a = vz.event_measure(S, evG)
    b = vz.event_measure(S, None)
    return dict(rms_fn=r["fn"], rms_mn=r["mn"], evG_fn=a["fn"], evG_mn=a["mn"], evOwn_fn=b["fn"], evOwn_mn=b["mn"], evG_pk_fn=a["pk_fn"], evG_pk_mn=a["pk_mn"], evOwn_n=[b["n_zone"], b["n_ref"]])


# ---- (2) 3D grid
gr = np.round(np.arange(0.0, 1.2001, 0.1), 1)
hits = {q: [] for q in ("rms", "evG", "evOwn", "pkG")}
best = {q: (9.0, None) for q in hits}
n = 0
allrows = []
for an in gr:
    for am in gr:
        for af in gr:
            r = meas(dict(nea=an, mid=am, far=af))
            n += 1
            allrows.append([float(an), float(am), float(af)] + [None if r[q] is None else round(float(r[q]), 3) for q in ("rms_fn", "rms_mn", "evG_fn", "evG_mn", "evOwn_fn", "evOwn_mn", "evG_pk_fn", "evG_pk_mn")])
            for q, (fk, mk) in dict(rms=("rms_fn", "rms_mn"), evG=("evG_fn", "evG_mn"), evOwn=("evOwn_fn", "evOwn_mn"), pkG=("evG_pk_fn", "evG_pk_mn")).items():
                fn, mn = r[fk], r[mk]
                if fn is None or mn is None:
                    continue
                if 0.40 <= fn <= 0.56 and mn >= 0.90:
                    hits[q].append([float(an), float(am), float(af), round(fn, 2), round(mn, 2)])
                d = abs(np.log(fn / 0.48)) + abs(np.log(min(mn, 0.97) / 0.97))
                if d < best[q][0]:
                    best[q] = (float(d), [float(an), float(am), float(af), round(fn, 2), round(mn, 2)])
OUT["grid3d"] = dict(n=n, alpha_values=[float(v) for v in gr], target="far/near zone/ref 0.40-0.56 and mid/near zone/ref >= 0.90",
                     n_hits={q: len(v) for q, v in hits.items()}, hits={q: v[:40] for q, v in hits.items()}, closest={q: v[1] for q, v in best.items()})
np.save(os.path.join(HERE, "p3_grid3d.npy"), np.array(allrows, dtype=float))
# for the hits: are they physically sensible (multiples on near weaker than on mid)?
OUT["grid3d"]["hits_with_near_ge_mid_ge_far"] = {q: [h for h in v if h[0] >= h[1] >= h[2]] for q, v in hits.items()}

# ---- (3) shift of the multiples against the primaries
A = dict(nea=0.68 / 0.928, mid=0.37 / 0.596, far=0.28 / 0.509)
sh = {}
for nm, al in (("alpha1", dict(nea=1.0, mid=1.0, far=1.0)), ("setA", A)):
    rr = []
    for s in range(-60, 61, 10):
        r = meas(al, s)
        rr.append(dict(shift_ms=s, rms_fn=r["rms_fn"], evG_fn=r["evG_fn"], evOwn_fn=r["evOwn_fn"], rms_mn=r["rms_mn"], evG_mn=r["evG_mn"], evOwn_n=r["evOwn_n"]))
    sh[nm] = rr
OUT["shifted_multiples"] = sh
json.dump(OUT, open(os.path.join(HERE, "p3_results.json"), "w"), indent=1, default=float)

print("copy measure per parent (coef, lag ms):")
for nm, r in rows.items():
    print(" ", nm)
    for q in r:
        print("    %d  near %s / %s   mid %s / %s   far %s / %s" % (q["parent_ms"], q["nea"]["first"], q["nea"]["second"], q["mid"]["first"], q["mid"]["second"], q["far"]["first"], q["far"]["second"]))
print(json.dumps(summ, default=lambda v: round(float(v), 2)))
print("\n3D grid: %d combinations; hits (target %s): %s" % (n, OUT["grid3d"]["target"], OUT["grid3d"]["n_hits"]))
for q, v in hits.items():
    print("  ", q, v[:12])
print("  closest:", OUT["grid3d"]["closest"])
print("  hits with near >= mid >= far:", {q: len(v) for q, v in OUT["grid3d"]["hits_with_near_ge_mid_ge_far"].items()})
f = lambda v: "  nan" if v is None else "%5.2f" % v
print("\nshifted multiples (shift ms: far/near RMS evG evOwn | mid/near RMS evG):")
for nm, rr in sh.items():
    print(" ", nm)
    for r in rr:
        print("    %4d : %s %s %s | %s %s  n %s" % (r["shift_ms"], f(r["rms_fn"]), f(r["evG_fn"]), f(r["evOwn_fn"]), f(r["rms_mn"]), f(r["evG_mn"]), r["evOwn_n"]))
