
# ---- fm_sl_iface.py: what would the top and base of the softest clay layer in 15/9-13 do on near, mid and far? ----
# Run:  cat sl_common.py fm_sl_iface.py | ssh linux-pc '<geoscience python> -I -'
# Single interfaces, exact Zoeppritz averaged over the offsets of each stack, plus a VTI term (Rueger 1997):
#     dR = 0.5*d_delta*sin2 + 0.5*d_epsilon*sin2*tan2
w, t, L = well()
vel = Vel(w)
vp, rho, gr, md = L["vp"], L["rho_used"], L["gr"], L["md"]


def med(a, b):
    m = (t >= a) & (t < b)
    return float(np.nanmedian(vp[m])), float(np.nanmedian(rho[m])), float(np.nanmin(md[m])), float(np.nanmax(md[m]))


ABOVE = med(1350, 1400)
SOFT = med(1412, 1446)
BELOW = med(1454, 1510)
print("== layers from the log (checkshot-corrected sonic, measured density)")
for nm, x in (("above (1350-1400 ms)", ABOVE), ("softest layer (1412-1446 ms)", SOFT), ("below (1454-1510 ms)", BELOW)):
    print(f"  {nm:30s} MD {x[2]:.0f}-{x[3]:.0f} m: Vp {x[0]:.0f}, rho {x[1]:.3f}, Vs trend {vs_trend(x[0]):.0f} (Vp/Vs {x[0]/vs_trend(x[0]):.2f}), Vs Castagna {vs_castagna(x[0]):.0f} (Vp/Vs {x[0]/vs_castagna(x[0]):.2f})")


def stack_r(a1, b1, r1, a2, b2, r2, t0, dd=0.0, de=0.0):
    g = stack_geometry(vel, t0)
    out = {}
    for k in g:
        s = g[k][0]
        s2 = s ** 2
        t2 = s2 / (1 - s2)
        out[k] = float(np.mean(zoep(a1, b1, r1, a2, b2, r2, s) + 0.5 * dd * s2 + 0.5 * de * s2 * t2))
    return out


# reference: what do ordinary interfaces in this well give with trend Vs? 20 ms blocks, 1270-1960 ms
ref = []
tb = np.arange(1270, 1960, 20)
bv = [(np.nanmedian(vp[(t >= a) & (t < a + 20)]), np.nanmedian(rho[(t >= a) & (t < a + 20)])) for a in tb]
for i in range(len(tb) - 1):
    (a1, r1), (a2, r2) = bv[i], bv[i + 1]
    r0 = 0.5 * np.log(a2 * r2 / (a1 * r1))
    if abs(r0) < 0.015:
        continue
    rr = stack_r(a1, vs_trend(a1), r1, a2, vs_trend(a2), r2, (tb[i] + 20) / 1000.0)
    ref.append((tb[i] + 20, r0, rr["mid"] / rr["nea"], rr["far"] / rr["nea"]))
ref = np.array(ref)
print(f"\n== reference: {len(ref)} log interfaces with |R0| > 0.015 between 1270 and 1960 ms, trend Vs on both sides")
print(f"   mid/near median {np.median(ref[:,2]):.3f} (p10 {np.percentile(ref[:,2],10):.3f}, p90 {np.percentile(ref[:,2],90):.3f});  far/near median {np.median(ref[:,3]):.3f} (p10 {np.percentile(ref[:,3],10):.3f}, p90 {np.percentile(ref[:,3],90):.3f})")
RM, RF = float(np.median(ref[:, 2])), float(np.median(ref[:, 3]))

MEAS = {"top": (1.06, 0.73), "base": (0.81, 0.45)}   # measured mid/near and far/near relative to the median event (fm_sl_events.py)
for nm, (U, D, t0) in (("top of the softest layer (event near 1386-1410 ms)", (ABOVE, SOFT, 1.41)), ("base of the softest layer (event at 1456 ms)", (SOFT, BELOW, 1.45))):
    side = "top" if nm.startswith("top") else "base"
    a1, r1, a2, r2 = U[0], U[1], D[0], D[1]
    print(f"\n== {nm}: R0 = {0.5*np.log(a2*r2/(a1*r1)):+.4f}; measured relative to the median event: mid/near {MEAS[side][0]:.2f}, far/near {MEAS[side][1]:.2f}")
    print("   case                                   |  R near   R mid    R far  | mid/near far/near | relative to reference: mid/near far/near")
    cases = []
    for red in (0, 5, 10, 15, 20, 25, 30):
        f = 1 - red / 100
        b1 = vs_trend(a1) * (f if side == "base" else 1)
        b2 = vs_trend(a2) * (f if side == "top" else 1)
        cases.append((f"isotropic, soft-layer Vs {red:2d} % below trend", b1, b2, 0.0, 0.0))
    b1c, b2c = vs_castagna(a1), vs_castagna(a2)
    cases.append(("isotropic, Castagna line both sides", b1c, b2c, 0.0, 0.0))
    if side == "base":
        cases.append(("soft layer Castagna, layer below trend", vs_castagna(a1), vs_trend(a2), 0.0, 0.0))
    else:
        cases.append(("layer above trend, soft layer Castagna", vs_trend(a1), vs_castagna(a2), 0.0, 0.0))
    sgn = 1 if side == "top" else -1
    for de, dd in ((0.05, 0.01), (0.10, 0.02), (0.15, 0.03), (0.20, 0.05)):
        cases.append((f"trend Vs, soft layer more anisotropic: d_eps {de:.2f}, d_delta {dd:.2f}", vs_trend(a1), vs_trend(a2), sgn * dd, sgn * de))
    for lab, b1, b2, dd, de in cases:
        rr = stack_r(a1, b1, r1, a2, b2, r2, t0, dd, de)
        mn, fn_ = rr["mid"] / rr["nea"], rr["far"] / rr["nea"]
        print(f"   {lab:60s} | {rr['nea']:+.4f} {rr['mid']:+.4f} {rr['far']:+.4f} | {mn:7.3f} {fn_:7.3f} | {mn/RM:14.2f} {fn_/RF:8.2f}")
