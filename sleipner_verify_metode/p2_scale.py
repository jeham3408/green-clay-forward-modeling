"""Points 2 and 3: scale the surface multiples of the experiment down to the strength seen in the real data.

Input: the experiment's stacks.npz (FD stacks with sea surface = FS, without = G; the multiples alone = FS - G).
Scaled stack per offset class k:  S_k = G_k + alpha_k (FS_k - G_k).
Measures: the experiment's per-reflector measure (slp.measure_event / summarise, unchanged) on the true reflectors
(picked on G) and on reflectors picked on the scaled full stack; and an RMS measure in the windows (own).
Also: how special is the zone? The same RMS measure in an 84 ms window slid through 1100-1900 ms.
Writes p2_results.json and prints tables.
"""
import json
import os

import numpy as np

import slp
import vz

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, "stacks.npz"))
T = Z["t"]
assert len(T) == len(vz.T_OUT)
KS = ("nea", "mid", "far")
REAL = dict(nea=(-0.68, 0.35), mid=(-0.37, 0.28), far=(-0.28, 0.12))        # copies measured on the real stacks (parent session)
OBS = dict(fn=0.48, fn_filtered=(0.62, 0.65), mn_main=(0.94, 1.00), mn_other=(0.60, 0.75))


def get(cfg, case, var):
    return {k: Z["%s_%s_%s__%s" % (cfg, case, var, k)].astype(float) for k in KS}


def scaled(G, FS, al):
    S = {k: G[k] + al[k] * (FS[k] - G[k]) for k in KS}
    S["ful"] = sum(slp.W_FULL[k] * S[k] for k in KS)
    return S


def measures(S, evG):
    r = vz.rms_measure(S)
    a = vz.event_measure(S, evG)
    b = vz.event_measure(S, None)
    return dict(rms_fn=r["fn"], rms_mn=r["mn"], rms_near_zone_over_ref=r["near_zone_over_ref"], rms_mid_zone_over_ref=r["mid_zone_over_ref"], rms_far_zone_over_ref=r["far_zone_over_ref"],
                evG_fn=a["fn"], evG_mn=a["mn"], evG_pk_fn=a["pk_fn"], evG_pk_mn=a["pk_mn"], evG_zone_values=a["fn_zone_values"], evG_n=[a["n_zone"], a["n_ref"]],
                evOwn_fn=b["fn"], evOwn_mn=b["mn"], evOwn_pk_fn=b["pk_fn"], evOwn_pk_mn=b["pk_mn"], evOwn_zone_values=b["fn_zone_values"], evOwn_n=[b["n_zone"], b["n_ref"]], evOwn_events_ms=b["events_ms"])


def copies_like_filter(S, parents):
    """First and second copy the way the experiment's filter measures them: on the stack itself, below the parents."""
    out = {}
    for k in KS:
        c1 = [slp.copy_coef(S[k], i, 0.080, 0.150, sign=-1) for i in parents]
        c2 = [slp.copy_coef(S[k], i, 0.190, 0.300, sign=+1) for i in parents]
        out[k] = dict(a1=float(np.mean([c[0] for c in c1])), a2=float(np.mean([c[0] for c in c2])), lag1_ms=float(1000 * np.mean([c[1] for c in c1])),
                      per_parent_a1=[round(c[0], 3) for c in c1])
    return out


def sliding(S, step=0.010):
    """RMS far/near and mid/near, window / reference, for an 84 ms window slid through the section. The reference
    excludes the window with the same margins as the real zone (28 ms before, 68 ms after)."""
    rows = []
    L = vz.ZONE[1] - vz.ZONE[0]
    for a in np.arange(1.130, 1.900 - L - 0.03 + 1e-9, step):
        r = vz.rms_measure(S, zone=(a, a + L), excl=(a - 0.028, a + L + 0.068))
        rows.append((round(1000 * a), r["fn"], r["mn"]))
    return rows


OUT = dict(real_copies=REAL, observed=OBS)
G, FS = get("G", "trend", "main"), get("FS", "trend", "main")
G["ful"] = sum(slp.W_FULL[k] * G[k] for k in KS)
FS["ful"] = sum(slp.W_FULL[k] * FS[k] for k in KS)
evG = slp.pick_events(G["ful"])
evFS = slp.pick_events(FS["ful"])
OUT["events_ms"] = dict(G=[int(round(1000 * i * slp.DT)) for i in evG], FS=[int(round(1000 * i * slp.DT)) for i in evFS])

# ---- model copies: true copies (multiples-only trace against the primary wavelet) for the strong primaries
M = {k: FS[k] - G[k] for k in KS}
tc = {}
for i in evG:
    if 1.15 <= i * slp.DT <= 1.36:
        tc[int(round(1000 * i * slp.DT))] = {k: [round(v, 3) for v in slp.copy_coef(G[k], i, 0.080, 0.150, sign=-1, y=M[k])] for k in KS}
OUT["model_true_first_copies"] = tc
parents = sorted([i for i in evFS if 1.24 <= i * slp.DT <= 1.38], key=lambda i: -abs(FS["ful"][i]))[:2]
OUT["filter_parents_ms"] = [int(round(1000 * i * slp.DT)) for i in parents]
OUT["model_copies_like_filter_alpha1"] = copies_like_filter(FS, parents)

# ---- alpha sets
p1156 = tc.get(1156) or list(tc.values())[0]
A = {}
A["alpha1_as_experiment"] = dict(nea=1.0, mid=1.0, far=1.0)
A["alpha0_no_multiples"] = dict(nea=0.0, mid=0.0, far=0.0)
A["A_ratio_to_true_copy_1156"] = {k: REAL[k][0] / p1156[k][0] for k in KS}
two = [tc[t] for t in tc if t in (1156, 1176)]
A["B_ratio_to_mean_true_copy_1156_1176"] = {k: REAL[k][0] / float(np.mean([q[k][0] for q in two])) for k in KS}
A["D_near_from_single_reflector_check_0p757"] = dict(nea=0.68 / 0.757, mid=A["A_ratio_to_true_copy_1156"]["mid"], far=A["A_ratio_to_true_copy_1156"]["far"])
# C: self-consistent. alpha per stack so that the copy measured the filter's way on the scaled stack equals the real one
grid = np.round(np.arange(0.0, 1.5001, 0.01), 2)
selfc, curve = {}, {}
for k in KS:
    vals = []
    for a in grid:
        S = {k: G[k] + a * M[k]}
        c1 = [slp.copy_coef(S[k], i, 0.080, 0.150, sign=-1)[0] for i in parents]
        vals.append(float(np.mean(c1)))
    vals = np.array(vals)
    j = int(np.argmin(np.abs(vals - REAL[k][0])))
    selfc[k] = float(grid[j])
    curve[k] = {str(a): round(float(v), 3) for a, v in zip(grid[::10], vals[::10])}
A["C_selfconsistent_filter_measure"] = selfc
OUT["copy_measured_filter_way_vs_alpha"] = curve
A["E_uniform_0p25"] = dict(nea=0.25, mid=0.25, far=0.25)
A["E_uniform_0p50"] = dict(nea=0.50, mid=0.50, far=0.50)
A["E_uniform_0p75"] = dict(nea=0.75, mid=0.75, far=0.75)
A["F_near_only_A"] = dict(nea=A["A_ratio_to_true_copy_1156"]["nea"], mid=0.0, far=0.0)
A["F_far_only_A"] = dict(nea=0.0, mid=0.0, far=A["A_ratio_to_true_copy_1156"]["far"])
OUT["alpha_sets"] = A

RES = {}
for nm, al in A.items():
    S = scaled(G, FS, al)
    r = measures(S, evG)
    r["copies_like_filter"] = {k: v["a1"] for k, v in copies_like_filter(S, parents).items()}
    r["copies2_like_filter"] = {k: v["a2"] for k, v in copies_like_filter(S, parents).items()}
    RES[nm] = r
OUT["main"] = RES

# ---- the same scaling in the other variants and sea-bed cases (sets A and C)
VAR = {}
for var in ("main", "stretch50", "hyper", "w35"):
    for case in ("trend", "trend_sbl", "trend_sbh", "m15"):
        try:
            g, f = get("G", case, var), get("FS", case, var)
        except KeyError:
            continue
        g["ful"] = sum(slp.W_FULL[k] * g[k] for k in KS)
        ev = slp.pick_events(get_ful := sum(slp.W_FULL[k] * get("G", "trend", var)[k] for k in KS))
        row = {}
        for nm in ("alpha1_as_experiment", "A_ratio_to_true_copy_1156", "C_selfconsistent_filter_measure", "D_near_from_single_reflector_check_0p757"):
            r = measures(scaled(g, f, A[nm]), ev)
            row[nm] = {q: r[q] for q in ("rms_fn", "evG_fn", "evOwn_fn", "evG_pk_fn", "rms_mn", "evG_mn", "evOwn_mn", "evG_pk_mn", "evOwn_n")}
        VAR["%s|%s" % (var, case)] = row
OUT["variants"] = VAR

# ---- grid over alpha_near and alpha_far (alpha_mid = geometric mean): which strengths give the observed 0.48?
GR = []
for an in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
    for af in (0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0):
        al = dict(nea=an, mid=float(np.sqrt(an * af)), far=af)
        r = measures(scaled(G, FS, al), evG)
        GR.append(dict(a_near=an, a_far=af, rms_fn=r["rms_fn"], evG_fn=r["evG_fn"], evOwn_fn=r["evOwn_fn"], rms_mn=r["rms_mn"], evG_mn=r["evG_mn"]))
OUT["grid_alpha"] = GR

# ---- how special is the zone: sliding window
SL = {}
for nm in ("alpha0_no_multiples", "alpha1_as_experiment", "A_ratio_to_true_copy_1156", "C_selfconsistent_filter_measure"):
    rows = sliding(scaled(G, FS, A[nm]))
    a = np.array(rows)
    jz = int(np.argmin(np.abs(a[:, 0] - 1408)))
    SL[nm] = dict(windows=len(rows), fn_min=float(a[:, 1].min()), fn_at_ms_of_min=float(a[int(np.argmin(a[:, 1])), 0]), fn_p10=float(np.percentile(a[:, 1], 10)), fn_median=float(np.median(a[:, 1])),
                  fn_p90=float(np.percentile(a[:, 1], 90)), fn_max=float(a[:, 1].max()), fn_zone_window=float(a[jz, 1]), rank_of_zone_from_lowest=int((a[:, 1] < a[jz, 1]).sum()) + 1,
                  share_windows_fn_below_0p65=float(np.mean(a[:, 1] < 0.65)), mn_min=float(a[:, 2].min()), mn_median=float(np.median(a[:, 2])), mn_zone_window=float(a[jz, 2]),
                  table=[[int(r[0]), round(r[1], 2), round(r[2], 2)] for r in rows])
OUT["sliding_window"] = SL
json.dump(OUT, open(os.path.join(HERE, "p2_results.json"), "w"), indent=1, default=float)

# ---- print
f = lambda v: "  nan" if v is None else "%5.2f" % v
print("events:", OUT["events_ms"], "filter parents:", OUT["filter_parents_ms"])
print("model true first copies (coef, lag s):", json.dumps(tc))
print("model copies the filter's way at alpha 1:", json.dumps(OUT["model_copies_like_filter_alpha1"], default=lambda v: round(float(v), 3)))
print("copy (filter's way) against alpha:", json.dumps(curve))
print("\n%-44s %-17s | far/near zone/ref: RMS evG evOwn pkG | mid/near: RMS evG evOwn pkG | near z/r, far z/r | copies measured" % ("set", "alpha n/m/f"))
for nm, r in RES.items():
    al = A[nm]
    print("%-44s %.2f/%.2f/%.2f    | %s %s %s %s | %s %s %s %s | %s %s | %s  zoneG %s nOwn %s" % (
        nm, al["nea"], al["mid"], al["far"], f(r["rms_fn"]), f(r["evG_fn"]), f(r["evOwn_fn"]), f(r["evG_pk_fn"]), f(r["rms_mn"]), f(r["evG_mn"]), f(r["evOwn_mn"]), f(r["evG_pk_mn"]),
        f(r["rms_near_zone_over_ref"]), f(r["rms_far_zone_over_ref"]), {k: round(v, 2) for k, v in r["copies_like_filter"].items()}, r["evG_zone_values"], r["evOwn_n"]))
print("\nvariants (far/near RMS evG evOwn pkG | mid/near RMS evG evOwn pkG):")
for kk, row in VAR.items():
    for nm, r in row.items():
        print("%-22s %-42s %s %s %s %s | %s %s %s %s  n %s" % (kk, nm, f(r["rms_fn"]), f(r["evG_fn"]), f(r["evOwn_fn"]), f(r["evG_pk_fn"]), f(r["rms_mn"]), f(r["evG_mn"]), f(r["evOwn_mn"]), f(r["evG_pk_mn"]), r["evOwn_n"]))
print("\ngrid (a_near a_far: far/near RMS evG evOwn | mid/near RMS evG):")
for g in GR:
    print("  %.1f %.1f : %s %s %s | %s %s" % (g["a_near"], g["a_far"], f(g["rms_fn"]), f(g["evG_fn"]), f(g["evOwn_fn"]), f(g["rms_mn"]), f(g["evG_mn"])))
print("\nsliding window (RMS far/near window/reference):")
for nm, s in SL.items():
    print("%-40s min %.2f (at %d ms) p10 %.2f median %.2f p90 %.2f max %.2f | zone window %.2f rank %d of %d | share < 0.65: %.2f | mid/near min %.2f median %.2f zone %.2f" % (
        nm, s["fn_min"], s["fn_at_ms_of_min"], s["fn_p10"], s["fn_median"], s["fn_p90"], s["fn_max"], s["fn_zone_window"], s["rank_of_zone_from_lowest"], s["windows"], s["share_windows_fn_below_0p65"],
        s["mn_min"], s["mn_median"], s["mn_zone_window"]))
    print("   ", " ".join("%d:%.2f" % (r[0], r[1]) for r in s["table"][::2]))
