"""Point 1: is there a low-Vs geometry (at most 20 % below trend) that lowers far/near in the zone to 0.5-0.65 of the reference?

Primaries only: own Zoeppritz at every 1 m interface, own rays, near / mid / far stack reflectivity, 25 Hz Ricker.
Two versions of every scenario: without NMO stretch ("nmofree") and with exact NMO of every offset trace ("nmo").
Measures: RMS in the windows (own), and the experiment's per-reflector measure on a fixed list of reflectors (picked on
the unperturbed full stack) and on reflectors picked on the scenario's own full stack.
Writes p1_results.json and prints tables.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import lsq_linear

import vz

HERE = os.path.dirname(os.path.abspath(__file__))
t_start = time.time()
m = vz.load_model()
z = m["z"]
F = vz.Forward(m)
n = len(z)
SOFT = (1358.0, 1390.0)
zone_z = [float(np.interp(t, m["twt"], m["ze"])) for t in vz.ZONE]          # depth of the zone limits
OUT = dict(zone_depth_m=zone_z, soft_layer_m=SOFT, soft_layer_twt_ms=[float(1000 * np.interp(a, m["ze"], m["twt"])) for a in SOFT])

# ---- model in and around the zone (8 m blocks), and the angles
prof = []
for a in np.arange(1280.0, 1500.0, 8.0):
    c = (z >= a) & (z < a + 8)
    prof.append(dict(z=float(a), twt_ms=round(float(1000 * np.interp(a, m["ze"], m["twt"])), 1), vp=round(float(m["vp"][c].mean())), vs=round(float(m["vs"][c].mean())),
                     rho=round(float(m["rho"][c].mean()), 3), ai=round(float((m["vp"][c] * m["rho"][c]).mean()))))
OUT["profile_8m"] = prof
k_top = int(np.searchsorted(m["ze"], SOFT[0])) - 1
ang = {}
for x in (150, 750, 775, 1750, 1775, 2500, 3150):
    ix = int(np.argmin(np.abs(vz.OFF - x)))
    ang[str(x)] = None if not np.isfinite(F.ANG[k_top, ix]) else round(float(F.ANG[k_top, ix]), 1)
OUT["angle_deg_at_soft_top_by_offset"] = ang
OUT["live_offsets_at_soft_top"] = {g: int(F.live[k_top, gm].sum()) for g, gm in F.gm.items()}
OUT["mean_sin2_at_soft_top"] = {g: float(np.mean(np.sin(np.radians(F.ANG[k_top, gm][F.live[k_top, gm]])) ** 2)) for g, gm in F.gm.items()}

# ---- single interface (blocked media around the base of the soft layer): stack-averaged R, far over near
def blk(a, b):
    c = (z >= a) & (z < b)
    return float(1 / np.mean(1 / m["vp"][c])), float(m["vs"][c].mean()), float(m["rho"][c].mean())
k_base = int(np.searchsorted(m["ze"], SOFT[1])) - 1
up, lo, ab = blk(*SOFT), blk(SOFT[1], SOFT[1] + 30), blk(SOFT[0] - 30, SOFT[0])
single = {}
for nm, (m1, m2, kk, which) in {"base": (up, lo, k_base, 1), "top": (ab, up, k_top, 2)}.items():
    row = {}
    for X in (0.0, 0.07, 0.15, 0.20):
        a1 = list(m1); a2 = list(m2)
        if which == 1:
            a1[1] *= (1 - X)
        else:
            a2[1] *= (1 - X)
        r = {}
        for g, gm in F.gm.items():
            lv = F.live[kk, gm]
            r[g] = float(np.mean(vz.zoeppritz(a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], F.P[kk, gm][lv]).real))
        row["%d" % round(100 * X)] = dict(R_near=r["nea"], R_mid=r["mid"], R_far=r["far"], far_over_near=r["far"] / r["nea"])
    single[nm] = dict(media=dict(upper=m1, lower=m2), by_vs_reduction_pct=row)
OUT["single_interface_blocked"] = single

# ---- scenarios
def box(z1, z2, X):
    f = np.ones(n)
    f[(z >= z1) & (z < z2)] = 1 - X
    return f

def ramp_top(z0, z1, X):          # gradual top (from z0), sharp base at z1
    f = np.ones(n)
    c = (z >= z0) & (z < z1)
    f[c] = 1 - X * (z[c] - z0) / (z1 - z0)
    return f

def ramp_base(z1, z2, X):         # sharp top at z1, gradual return to the trend at z2
    f = np.ones(n)
    c = (z >= z1) & (z < z2)
    f[c] = 1 - X * (z2 - z[c]) / (z2 - z1)
    return f

def tri(zc, half, X):
    f = np.ones(n)
    c = np.abs(z - zc) < half
    f[c] = 1 - X * (1 - np.abs(z[c] - zc) / half)
    return f

ai = m["vp"] * m["rho"]
def runmean(x, L):
    k = np.ones(L) / L
    return np.convolve(np.pad(x, L // 2, mode="edge"), k, mode="valid")[:len(x)]

def softbeds(za, zb, X, L=61, prop=False, key="ai"):
    """Vs lowered in every bed that is softer than the running mean over L m, between za and zb.
    prop: the reduction is proportional to how soft the bed is (full X at the 90th percentile of the softness)."""
    q = ai if key == "ai" else m["vp"]
    d = runmean(q, L) - q
    c = (z >= za) & (z < zb)
    f = np.ones(n)
    if prop:
        s = np.percentile(d[c & (d > 0)], 90)
        f[c] = 1 - X * np.clip(d[c] / s, 0, 1)
    else:
        f[c & (d > 0)] = 1 - X
    return f

zt, zb_ = zone_z
SC = {}
for X in (0.07, 0.15, 0.20, 0.30, 0.40):
    SC["box32_softlayer_%02d" % round(100 * X)] = box(SOFT[0], SOFT[1], X)
for X in (0.07, 0.15, 0.20):
    p = "%02d" % round(100 * X)
    SC["box_wholezone_" + p] = box(zt, zb_, X)
    for T in (60, 90, 150):
        SC["box%d_topfixed_%s" % (T, p)] = box(SOFT[0], SOFT[0] + T, X)
        SC["box%d_basefixed_%s" % (T, p)] = box(SOFT[1] - T, SOFT[1], X)
        SC["box%d_centred_%s" % (T, p)] = box(0.5 * sum(SOFT) - T / 2, 0.5 * sum(SOFT) + T / 2, X)
    for T in (32, 60, 90):
        SC["above%d_%s" % (T, p)] = box(SOFT[0] - T, SOFT[0], X)
        SC["below%d_%s" % (T, p)] = box(SOFT[1], SOFT[1] + T, X)
    for Lr in (32, 90, 150, 300):
        SC["gradtop%d_sharpbase_%s" % (Lr, p)] = ramp_top(SOFT[1] - Lr, SOFT[1], X)
        SC["sharptop_gradbase%d_%s" % (Lr, p)] = ramp_base(SOFT[0], SOFT[0] + Lr, X)
    for h in (45, 90, 150):
        SC["triangle_half%d_%s" % (h, p)] = tri(0.5 * sum(SOFT), h, X)
    SC["all_above_softbase_from1000m_%s" % p] = box(1000.0, SOFT[1], X)
    SC["all_below_softtop_to2000m_%s" % p] = box(SOFT[0], 2000.0, X)
    for L in (31, 61, 101):
        SC["softbeds_zone_L%d_%s" % (L, p)] = softbeds(zt, zb_, X, L)
        SC["softbeds_zone_pm30_L%d_%s" % (L, p)] = softbeds(zt - 30, zb_ + 30, X, L)
    SC["softbeds_zone_prop_L61_%s" % p] = softbeds(zt, zb_, X, 61, prop=True)
    SC["softbeds_zone_pm30_prop_L61_%s" % p] = softbeds(zt - 30, zb_ + 30, X, 61, prop=True)
    SC["softbeds_zone_pm30_lowVp_L61_%s" % p] = softbeds(zt - 30, zb_ + 30, X, 61, key="vp")

BASE = {False: F.stacks(F.R0, stretch=False), True: F.stacks(F.R0, stretch=True)}
import slp
EV0 = {s: slp.pick_events(BASE[s]["ful"]) for s in (False, True)}

def evaluate(fac, stretch):
    S = F.stacks(F.refl(fac), stretch=stretch)
    r = vz.rms_measure(S)
    e1 = vz.event_measure(S, EV0[stretch])
    e2 = vz.event_measure(S, None)
    return dict(rms_fn=r["fn"], rms_mn=r["mn"], rms_far_zone_over_ref=r["far_zone_over_ref"], rms_near_zone_over_ref=r["near_zone_over_ref"],
                ev_fixed_fn=e1["fn"], ev_fixed_mn=e1["mn"], ev_fixed_pk_fn=e1["pk_fn"], ev_fixed_zone_values=e1["fn_zone_values"], ev_fixed_n=[e1["n_zone"], e1["n_ref"]],
                ev_own_fn=e2["fn"], ev_own_mn=e2["mn"], ev_own_n=[e2["n_zone"], e2["n_ref"]], ev_own_events_ms=e2["events_ms"])

RES = {}
for st in (False, True):
    tag = "nmo" if st else "nmofree"
    RES[tag] = {"baseline": evaluate(None, st)}
    RES[tag]["baseline"]["events_ms"] = [int(round(1000 * i * vz.DT)) for i in EV0[st]]
    for nm, fac in SC.items():
        RES[tag][nm] = evaluate(fac, st)
OUT["scenarios"] = RES

# ---- scan of boxes: every top depth and thickness, 20 % and 7 %
scan = {}
for X in (0.20, 0.07):
    rows = []
    for T in (8, 12, 16, 24, 32, 46, 60, 76, 90, 120, 150, 200, 300, 500):
        for top in np.arange(1200.0, 1500.1, 4.0):
            fac = box(top, top + T, X)
            S = F.stacks(F.refl(fac), stretch=False)
            r = vz.rms_measure(S)
            e1 = vz.event_measure(S, EV0[False])
            rows.append((T, float(top), r["fn"], e1["fn"]))
    a = np.array(rows)
    b0 = RES["nmofree"]["baseline"]
    i1, i2 = int(np.argmin(a[:, 2])), int(np.argmin(a[:, 3]))
    scan["%02d" % round(100 * X)] = dict(n=len(rows), baseline_rms=b0["rms_fn"], baseline_ev=b0["ev_fixed_fn"],
                                         min_rms=dict(thickness=a[i1, 0], top=a[i1, 1], rms_fn=a[i1, 2], ev_fixed_fn=a[i1, 3]),
                                         min_ev=dict(thickness=a[i2, 0], top=a[i2, 1], rms_fn=a[i2, 2], ev_fixed_fn=a[i2, 3]),
                                         max_rms=float(a[:, 2].max()), max_ev=float(a[:, 3].max()),
                                         share_rms_below_0p8x_baseline=float(np.mean(a[:, 2] < 0.8 * b0["rms_fn"])),
                                         share_ev_below_0p8x_baseline=float(np.mean(a[:, 3] < 0.8 * b0["ev_fixed_fn"])))
    # the best boxes also with NMO stretch
    for key, i in (("min_rms", i1), ("min_ev", i2)):
        e = evaluate(box(a[i, 1], a[i, 1] + a[i, 0], X), True)
        scan["%02d" % round(100 * X)][key]["with_nmo"] = dict(rms_fn=e["rms_fn"], ev_fixed_fn=e["ev_fixed_fn"], ev_own_fn=e["ev_own_fn"])
OUT["box_scan"] = scan

# ---- bound: the Vs profile (0 to Xmax below trend, free per bed of B m) that minimises the far stack in the zone
zmask, rmask = vz.win_masks()
def far_zone(fac, stretch=False):
    return F.stacks(F.refl(fac), stretch=stretch)["far"][zmask]

bound = {}
for Xmax in (0.20, 0.07):
    for B in (4, 8, 16, 32):
        for rng_name, (za, zb2) in {"zone_only": (zt, zb_), "zone_pm60": (zt - 60, zb_ + 60)}.items():
            edges = np.arange(za, zb2 + 1e-6, B)
            nb = len(edges) - 1
            cells = [(z >= edges[j]) & (z < edges[j + 1]) for j in range(nb)]
            def fac_of(x):
                f = np.ones(n)
                for j in range(nb):
                    f[cells[j]] = 1 - x[j]
                return f
            x = np.zeros(nb)
            for it in range(3):
                b = far_zone(fac_of(x))
                A = np.zeros((len(b), nb))
                h = 0.02
                for j in range(nb):
                    xp = x.copy()
                    xp[j] = x[j] + h if x[j] + h <= Xmax else x[j] - h
                    A[:, j] = (far_zone(fac_of(xp)) - b) / (xp[j] - x[j])
                sol = lsq_linear(A, -(b - A @ x), bounds=(0.0, Xmax))
                x = np.clip(sol.x, 0, Xmax)
            fac = fac_of(x)
            e0, e1 = evaluate(fac, False), evaluate(fac, True)
            ai_b = np.array([ai[c].mean() for c in cells])
            cc = float(np.corrcoef(x, ai_b - runmean(ai, 61)[[int(np.nonzero(c)[0].mean()) for c in cells]])[0, 1]) if nb > 2 and x.std() > 0 else None
            bound["X%02d_B%d_%s" % (round(100 * Xmax), B, rng_name)] = dict(
                n_beds=nb, x_pct=[round(float(100 * v), 1) for v in x], bed_tops_m=[float(e) for e in edges[:-1]],
                share_at_bound=float(np.mean((x < 1e-4) | (x > Xmax - 1e-4))), mean_reduction_pct=float(100 * x.mean()),
                corr_reduction_with_ai_anomaly=cc,
                nmofree=dict(rms_fn=e0["rms_fn"], ev_fixed_fn=e0["ev_fixed_fn"], ev_own_fn=e0["ev_own_fn"], ev_fixed_pk_fn=e0["ev_fixed_pk_fn"], zone_values=e0["ev_fixed_zone_values"]),
                nmo=dict(rms_fn=e1["rms_fn"], ev_fixed_fn=e1["ev_fixed_fn"], ev_own_fn=e1["ev_own_fn"], ev_fixed_pk_fn=e1["ev_fixed_pk_fn"], zone_values=e1["ev_fixed_zone_values"]))
OUT["bound_free_profile"] = bound
OUT["runtime_s"] = time.time() - t_start
json.dump(OUT, open(os.path.join(HERE, "p1_results.json"), "w"), indent=1, default=float)

# ---- print
print("zone depth %s m; soft layer %s m = %s ms" % (np.round(zone_z, 1), SOFT, np.round(OUT["soft_layer_twt_ms"], 1)))
print("angles at soft top:", ang, "live offsets:", OUT["live_offsets_at_soft_top"], "mean sin2:", {k: round(v, 3) for k, v in OUT["mean_sin2_at_soft_top"].items()})
print("profile (8 m blocks): z twt vp vs rho ai")
for p in prof:
    print("  %6.0f %7.1f %5d %4d %.3f %5d" % (p["z"], p["twt_ms"], p["vp"], p["vs"], p["rho"], p["ai"]))
print("single interface:", json.dumps(single, indent=None, default=lambda v: round(float(v), 4)))
f = lambda v: "  nan" if v is None else "%5.2f" % v
for tag in ("nmofree", "nmo"):
    print("\n=== %s  (zone/reference of far/near; mid/near in brackets)   baseline events %s" % (tag, RES[tag]["baseline"]["events_ms"]))
    print("%-42s %6s %6s %6s %6s | %6s %6s | zone values (fixed reflectors)" % ("scenario", "RMS", "evFix", "evOwn", "pkFix", "RMSmn", "evFmn"))
    for nm, r in RES[tag].items():
        print("%-42s %s %s %s %s | %s %s | %s n=%s own n=%s" % (nm, f(r["rms_fn"]), f(r["ev_fixed_fn"]), f(r["ev_own_fn"]), f(r["ev_fixed_pk_fn"]), f(r["rms_mn"]), f(r["ev_fixed_mn"]),
                                                             r["ev_fixed_zone_values"], r["ev_fixed_n"], r["ev_own_n"]))
print("\nbox scan:", json.dumps(scan, indent=1, default=lambda v: round(float(v), 3)))
print("\nbound (free profile):")
for k, b in bound.items():
    print("%-22s beds %3d mean red %4.1f %% at-bound %.2f corrAI %s | nmofree RMS %s evFix %s evOwn %s | nmo RMS %s evFix %s evOwn %s" % (
        k, b["n_beds"], b["mean_reduction_pct"], b["share_at_bound"], f(b["corr_reduction_with_ai_anomaly"]), f(b["nmofree"]["rms_fn"]), f(b["nmofree"]["ev_fixed_fn"]), f(b["nmofree"]["ev_own_fn"]),
        f(b["nmo"]["rms_fn"]), f(b["nmo"]["ev_fixed_fn"]), f(b["nmo"]["ev_own_fn"])))
print("runtime %.0f s" % OUT["runtime_s"])
