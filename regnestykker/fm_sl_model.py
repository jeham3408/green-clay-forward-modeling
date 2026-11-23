
# ---- fm_sl_model.py: forward model at 15/9-13 for near, mid and far, and comparison with the real traces ----
# Run:  cat sl_common.py fm_sl_model.py | ssh linux-pc '<geoscience python> -I -'
S = stacks()
w, t, L = well()
vel = Vel(w)
vp, rho, gr, md = L["vp"], L["rho_used"], L["gr"], L["md"]

SOFT = (1275.0, 1448.0)      # soft, low-density mudstone unit in the upper Hordaland (MD about 1255-1440 m)
GC = (1855.0, 1990.0)        # the part of the predicted Green Clay interval that the log and the cube reach


def scen(name):
    vs = vs_trend(vp)
    if name == "trend":
        return vs
    if name == "soft-5":
        return np.where((t >= SOFT[0]) & (t < SOFT[1]), 0.95 * vs, vs)
    if name == "soft-10":
        return np.where((t >= SOFT[0]) & (t < SOFT[1]), 0.90 * vs, vs)
    if name == "gc-5":
        return np.where(t >= GC[0], 0.95 * vs, vs)
    if name == "gc-10":
        return np.where(t >= GC[0], 0.90 * vs, vs)
    if name == "castagna":
        return vs_castagna(vp)
    raise ValueError(name)


print("== unit medians used in the model (trend Vs)")
for nm, (a, b) in (("above soft unit 1180-1270", (1180, 1270)), ("soft unit", SOFT), ("below soft unit 1450-1520", (1450, 1520)), ("1520-1850", (1520, 1850)), ("predicted Green Clay part", GC)):
    m = (t >= a) & (t < b)
    v, r_ = np.nanmedian(vp[m]), np.nanmedian(rho[m])
    print(f"  {nm:28s} TWT {a:.0f}-{b:.0f} ms, MD {np.nanmin(md[m]):.0f}-{np.nanmax(md[m]):.0f} m: Vp {v:.0f}, rho {r_:.3f}, GR {np.nanmedian(gr[m]):.0f}, "
          f"Vs trend {vs_trend(v):.0f} (Vp/Vs {v/vs_trend(v):.2f}), Vs Castagna {vs_castagna(v):.0f} (Vp/Vs {v/vs_castagna(v):.2f})")

print("\n== angle ranges of the stacks at the times that matter (model offsets near 150-761, mid 761-1763, far 1763-3150 m, mute 48 deg)")
for t0 in (1.30, 1.45, 1.60, 1.85):
    g = stack_geometry(vel, t0)
    lr = live_ranges(vel, t0)
    print(f"  t0 {t0:.2f} s: " + "; ".join(f"{k} {np.degrees(np.arcsin(vel.sin_angle(lr[k][0], t0))):.0f}-{np.degrees(np.arcsin(vel.sin_angle(lr[k][1], t0))):.0f} deg, mean sin2 {np.mean(g[k][0]**2):.3f}, stretch {np.mean(g[k][1]):.3f}" for k in ("nea", "mid", "far")))

names = ("trend", "soft-5", "soft-10", "gc-5", "gc-10", "castagna")
SY, RC = {}, {}
for nm in names:
    SY[nm], RC[nm] = synth(t, vp, rho, scen(nm), vel)

print("\n== single interface: base of the soft unit (medians 60 ms above and 60 ms below 1448 ms)")
ma, mb = (t >= 1388) & (t < 1448), (t >= 1452) & (t < 1512)
a1, r1, a2, r2 = np.nanmedian(vp[ma]), np.nanmedian(rho[ma]), np.nanmedian(vp[mb]), np.nanmedian(rho[mb])
g = stack_geometry(vel, 1.45)
print(f"  above: Vp {a1:.0f}, rho {r1:.3f}; below: Vp {a2:.0f}, rho {r2:.3f}; R0 = {0.5*np.log(a2*r2/(a1*r1)):+.4f}")
print("  Vs of the soft unit        |  R near   R mid    R far  | mid/near  far/near")
for lab, f1, fn in (("trend", 1.0, vs_trend), ("5 % below trend", 0.95, vs_trend), ("10 % below trend", 0.90, vs_trend), ("Castagna line", 1.0, vs_castagna)):
    b1, b2 = f1 * fn(a1), fn(a2)
    rr = {k: float(np.mean(zoep(a1, b1, r1, a2, b2, r2, g[k][0]))) for k in g}
    print(f"  {lab:26s} | {rr['nea']:+.4f} {rr['mid']:+.4f} {rr['far']:+.4f} | {rr['mid']/rr['nea']:7.3f} {rr['far']/rr['nea']:8.3f}   (Vs above {b1:.0f}, below {b2:.0f})")

print("\n== Nini-4 top Horda with the Sleipner stack geometry at 1.45 s (for scale): Lark 2438/1012/2.295 over Horda 2215/x/2.171")
for lab, b2 in (("measured Vs 811", 811.0), ("local trend Vs 855", 855.0), ("Vp/Vs of layer above 912", 912.0)):
    rr = {k: float(np.mean(zoep(2438.0, 1012.0, 2.295, 2215.0, b2, 2.171, g[k][0]))) for k in g}
    print(f"  {lab:26s} | {rr['nea']:+.4f} {rr['mid']:+.4f} {rr['far']:+.4f} | mid/near {rr['mid']/rr['nea']:.3f}  far/near {rr['far']/rr['nea']:.3f}")

# real traces at the well: mean of 5x5 traces
real = {k: np.asarray(S[k][WELL_IL - 2:WELL_IL + 3, WELL_XL - 2:WELL_XL + 3, :], np.float64).mean((0, 1)) for k in KEYS}
n = min(len(t), real["nea"].shape[0])
print("\n== well tie per stack (trend scenario): correlation real vs synthetic, lags +-40 ms")
print("  window ms   | stack |  cc0   best cc  lag ms | scale real/synthetic at best lag")
for win in ((1270, 1560), (1560, 1960), (1270, 1960)):
    i0, i1 = int(win[0] / 2), int(win[1] / 2)
    for k in ("nea", "mid", "far"):
        c0, cb, lag = xc(SY["trend"][k][i0:i1], real[k][i0:i1])
        Lg = int(round(lag / 2))
        a, b = SY["trend"][k][i0:i1], real[k][i0 + Lg:i1 + Lg]
        sc = (a * b).sum() / (a * a).sum()
        print(f"  {win[0]}-{win[1]} | {k:5s} | {c0:+.2f}  {cb:+.2f}  {lag:+6.0f} | {sc:8.2f}")
    # full stack against the fold-weighted sum of the three synthetics
    fw = np.array([0.2, 0.33, 0.45])
    sf = fw[0] * SY["trend"]["nea"] + fw[1] * SY["trend"]["mid"] + fw[2] * SY["trend"]["far"]
    c0, cb, lag = xc(sf[i0:i1], real["ful"][i0:i1])
    print(f"  {win[0]}-{win[1]} | full  | {c0:+.2f}  {cb:+.2f}  {lag:+6.0f} |")

print("\n== predicted amplitude ratios from the synthetics (RMS in window), per scenario")
print("  scenario  | window ms   | rms mid/near  rms far/near | change against trend: mid/near  far/near")
base = {}
for nm in names:
    for win in ((1280, 1440), (1400, 1500), (1500, 1840), (1840, 1988)):
        i0, i1 = int(win[0] / 2), int(win[1] / 2)
        r = {k: np.sqrt(np.mean(SY[nm][k][i0:i1] ** 2)) for k in ("nea", "mid", "far")}
        mn, fn_ = r["mid"] / r["nea"], r["far"] / r["nea"]
        if nm == "trend":
            base[win] = (mn, fn_)
        print(f"  {nm:9s} | {win[0]}-{win[1]} | {mn:10.3f} {fn_:13.3f} | {100*(mn/base[win][0]-1):+18.1f} % {100*(fn_/base[win][1]-1):+8.1f} %")

print("\n== strongest events on the real stacks near the well (mean of 5x5 traces), 1150-1990 ms")
ful = real["ful"]
rms = np.sqrt(np.mean(ful[575:995] ** 2))
print(f"  rms of the full stack in the window: {rms:.3f}")
print("   time ms | full  near   mid   far | synthetic (trend) near  mid   far")
for i in range(576, 994):
    if abs(ful[i]) > 1.6 * rms and abs(ful[i]) >= abs(ful[i - 1]) and abs(ful[i]) >= abs(ful[i + 1]):
        print(f"   {i*2:6d} | {ful[i]:+.2f} {real['nea'][i]:+.2f} {real['mid'][i]:+.2f} {real['far'][i]:+.2f} | {SY['trend']['nea'][i]:+18.3f} {SY['trend']['mid'][i]:+.3f} {SY['trend']['far'][i]:+.3f}")

print("\n== cube-wide mean trace (all traces): strongest flat events 1050-1990 ms")
mt = {k: np.asarray(S[k][::2, ::2, :], np.float64).mean((0, 1)) for k in KEYS}
r0 = np.sqrt(np.mean(mt["ful"][525:995] ** 2))
print("   time ms | full   near    mid    far   (mean over the whole cube)")
for i in range(526, 994):
    v = mt["ful"][i]
    if abs(v) > 1.5 * r0 and abs(v) >= abs(mt["ful"][i - 1]) and abs(v) >= abs(mt["ful"][i + 1]):
        print(f"   {i*2:6d} | {v:+.3f} {mt['nea'][i]:+.3f} {mt['mid'][i]:+.3f} {mt['far'][i]:+.3f}")
np.save(os.path.join(AVO, "synth_trend.npy"), np.stack([SY["trend"][k] for k in ("nea", "mid", "far")]))
