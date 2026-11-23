"""Primaries only, own Zoeppritz (e3_zoep.py, energy-conserving 4x4 solve): what does lowering Vs by 7, 15, 20 % in the
soft layer as built (cells 1358-1390 m) do to far/near, reflector by reflector?
(1) single interfaces from 30 m block averages: top (above -> soft) and base (soft -> below);
(2) the 1 m model as built, 1250-1560 m: reflectivity per interface and offset (own ray parameters), mean over the
    live offsets of each stack (48 deg mute), convolved in vertical time with a 25 Hz Ricker with vertical ghosts."""
import numpy as np

import e1_lib as L
import e3_zoep as Z

m = L.load_model()
off = np.arange(150.0, 3150.1, 25.0)
kin = L.Kin(m, off, smooth=11)
vp, vs, rho, soft = m["vp"], m["vs"], m["rho"], m["soft"]
i_s = np.where(soft)[0]
k0, k1 = int(i_s.min()), int(i_s.max()) + 1                 # soft cells k0..k1-1 ; top interface above cell k0, base below cell k1-1
R = dict(soft_cells=[k0, k1 - 1], t0_top_ms=float(1000 * kin.t0[k0 - 1]), t0_base_ms=float(1000 * kin.t0[k1 - 1]))
blk = lambda a, i, j: float(np.mean(a[i:j]))
B = {"above": (k0 - 30, k0), "soft": (k0, k1), "below": (k1, k1 + 30)}
P = {n: dict(vp=blk(vp, *ij), vs=blk(vs, *ij), rho=blk(rho, *ij)) for n, ij in B.items()}
R["blocks"] = P
print("blocks:", {n: {k: round(v, 1) if k != "rho" else round(v, 3) for k, v in x.items()} for n, x in P.items()})
sel = {}
for sk, (a, b) in L.RANGES.items():
    s = (off >= a - 0.5) & (off <= b + 0.5 if sk == "far" else off < b - 0.5) & np.isfinite(kin.px[k0 - 1]) & (kin.ang[k0 - 1] <= 48.0)
    sel[sk] = s
R["angles_at_soft_top_deg"] = {sk: [float(kin.ang[k0 - 1][s].min()), float(kin.ang[k0 - 1][s].mean()), float(kin.ang[k0 - 1][s].max()), int(s.sum())] for sk, s in sel.items()}
print("incidence angle at the soft-layer top (min, mean, max, n live offsets):", {k: [round(x, 1) for x in v] for k, v in R["angles_at_soft_top_deg"].items()})

CASES = {"trend": 0.0, "m07": 0.07, "m15": 0.15, "m20": 0.20}
single = {}
for cn, f in CASES.items():
    row = {}
    for iface, (u, l) in (("top", ("above", "soft")), ("base", ("soft", "below"))):
        a, b = dict(P[u]), dict(P[l])
        if u == "soft":
            a["vs"] *= 1 - f
        if l == "soft":
            b["vs"] *= 1 - f
        pk = kin.px[(k0 - 1) if iface == "top" else (k1 - 1)]
        r = np.array([Z.rpp(a["vp"], a["vs"], a["rho"], b["vp"], b["vs"], b["rho"], p).real if np.isfinite(p) else np.nan for p in pk])
        row[iface] = {sk: float(np.nanmean(r[s])) for sk, s in sel.items()}
        row[iface]["R0"] = float(Z.rpp(a["vp"], a["vs"], a["rho"], b["vp"], b["vs"], b["rho"], 0.0).real)
        row[iface]["far_over_near"] = row[iface]["far"] / row[iface]["near"]
        row[iface]["mid_over_near"] = row[iface]["mid"] / row[iface]["near"]
    single[cn] = row
R["single_interfaces"] = single
print("\nsingle interfaces (30 m block averages), mean Rpp per stack:")
for cn, row in single.items():
    for iface in ("top", "base"):
        x = row[iface]
        print("  %-5s %-4s R0=%+.4f near=%+.4f mid=%+.4f far=%+.4f  far/near=%+.2f  mid/near=%+.2f" % (cn, iface, x["R0"], x["near"], x["mid"], x["far"], x["far_over_near"], x["mid_over_near"]))
for cn in ("m07", "m15", "m20"):
    dt_, db_ = single[cn]["top"]["far"] - single["trend"]["top"]["far"], single[cn]["base"]["far"] - single["trend"]["base"]["far"]
    print("  %s: change of far Rpp at top %+.4f, at base %+.4f (sum %+.4f)" % (cn, dt_, db_, dt_ + db_))
    R.setdefault("far_change", {})[cn] = dict(top=dt_, base=db_, sum=dt_ + db_)

# (2) as built
ka, kb = 1200, 1620
tgrid = L.T0
def ghost_ricker(f0):
    tw = np.arange(-150, 151) * L.DT
    rk = lambda tt: (1 - 2 * (np.pi * f0 * tt) ** 2) * np.exp(-(np.pi * f0 * tt) ** 2)
    t14, t2 = 14.0 / 1480.0, 2.0 / 1480.0
    return rk(tw + t14) - rk(tw + t2) - rk(tw - t2) + rk(tw - t14)
wv = ghost_ricker(25.0)
syn = {}
for cn, f in CASES.items():
    v_s = vs.copy()
    v_s[soft] *= 1 - f
    tr = {sk: np.zeros(len(tgrid)) for sk in sel}
    for k in range(ka, kb):                                   # interface below cell k
        if abs(vp[k + 1] - vp[k]) + abs(v_s[k + 1] - v_s[k]) + abs(rho[k + 1] - rho[k]) < 1e-9:
            continue
        pk = kin.px[k]
        live = np.isfinite(pk) & (kin.ang[k] <= 48.0)
        r = np.full(len(off), np.nan)
        for ix in np.where(live)[0]:
            r[ix] = Z.rpp(vp[k], v_s[k], rho[k], vp[k + 1], v_s[k + 1], rho[k + 1], pk[ix]).real
        tt = kin.t0[k] / L.DT
        i0 = int(np.floor(tt)); fr = tt - i0
        for sk, (a, b) in L.RANGES.items():
            s = (off >= a - 0.5) & (off <= b + 0.5 if sk == "far" else off < b - 0.5) & live
            if s.sum() == 0:
                continue
            val = float(np.mean(r[s]))
            tr[sk][i0] += val * (1 - fr); tr[sk][i0 + 1] += val * fr
    syn[cn] = {sk: np.convolve(x, wv, mode="same") for sk, x in tr.items()}
np.savez_compressed(L.OUT + "/e5_zoep_stacks.npz", t=tgrid, **{"%s__%s" % (c, k): v for c, S in syn.items() for k, v in S.items()})
# reflectors of the zone: pick on the trend near trace inside 1380-1500 ms (model edge effects are far away)
w = L.win(1.380, 1.500)
pk = L.picks(syn["trend"]["near"], rel=0.25, tlim=(1.300, 1.540))
rows = {}
for cn in CASES:
    rr = L.per_reflector(syn[cn], pk, "far", "near")
    rows[cn] = rr
print("\nas built (1 m cells, primaries only, no stretch): far/near per reflector (signed peak / envelope)")
print("  t_ms   " + "   ".join("%-13s" % c for c in CASES))
for j, i in enumerate(pk):
    print("  %6.0f " % (1000 * tgrid[i]) + "   ".join("%+5.2f / %4.2f " % (rows[c][j]["pk"], rows[c][j]["env"]) for c in CASES))
R["as_built_far_over_near_per_reflector"] = rows
lev = {}
for cn in ("m07", "m15", "m20"):
    lev[cn] = {sk: dict(level_ratio=L.rms(syn[cn][sk], w) / L.rms(syn["trend"][sk], w), change_rel=L.rms(syn[cn][sk] - syn["trend"][sk], w) / L.rms(syn["trend"][sk], w)) for sk in sel}
    d = syn[cn]["far"] - syn["trend"]["far"]
    ww = L.win(1.380, 1.480)
    scale = np.abs(syn["trend"]["far"][w]).max()
    ipos, ineg = np.where(ww)[0][np.argmax(d[ww])], np.where(ww)[0][np.argmin(d[ww])]
    lev[cn]["far_difference_extrema"] = dict(max_rel=float(d[ipos] / scale), t_max_ms=float(1000 * tgrid[ipos]), min_rel=float(d[ineg] / scale), t_min_ms=float(1000 * tgrid[ineg]))
R["as_built_level_1380_1500"] = lev
print("\nas built, 1380-1500 ms: RMS level ratio to trend and relative RMS change of the waveform")
for cn, x in lev.items():
    print("  %s " % cn + "  ".join("%s: level %.3f change %.3f" % (sk, x[sk]["level_ratio"], x[sk]["change_rel"]) for sk in sel), " far difference extrema:", {k: round(v, 3) for k, v in x["far_difference_extrema"].items()})
# far/near RMS in the zone window against the same in 1200-1380 (reference inside the modelled interval)
zr = {}
for cn in CASES:
    S = syn[cn]
    zw, rw = L.win(*L.ZONE), L.win(1.250, 1.380) | L.win(1.560, 1.600)
    zr[cn] = (L.rms(S["far"], zw) / L.rms(S["near"], zw)) / (L.rms(S["far"], rw) / L.rms(S["near"], rw))
R["as_built_zone_over_ref_rms_far_near"] = zr
print("\nas built: zone/reference (RMS far/near; reference 1250-1380 and 1560-1600 ms):", {k: round(v, 3) for k, v in zr.items()})
L.jdump(R, "e5_results.json")
print("saved")
