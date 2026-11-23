
# ---- fm_sl_1523.py: local shear velocity in the Hordaland Group from well 15/9-23 (FORCE 2020, measured DTS),
# about 12 km from 15/9-13. How far below a local Vs-Vp trend do Hordaland mudstones go, and what does the
# Sleipner soft-layer base give with this trend? ----
# Run:  cat sl_common.py fm_sl_1523.py | ssh linux-pc '<geoscience python> -I -'
LAS = R + "/force2020/las/Force_2020_all_wells_train_test_blind_hidden_final/15_9-23.las"
cols, rows = None, []
for line in open(LAS, errors="replace"):
    s = line.strip()
    if s.startswith("~C"):
        cols, mode = [], "c"
        continue
    if s.startswith("~A"):
        mode = "a"
        continue
    if s.startswith("~"):
        mode = None
        continue
    if cols is not None and mode == "c" and s and not s.startswith("#"):
        cols.append(s.split(".")[0].strip())
    elif cols is not None and mode == "a" and s and not s.startswith("#"):
        p = s.split()
        if len(p) == len(cols):
            rows.append([float(v) for v in p])
A = np.array(rows)
A[A == -999.25] = np.nan
C = {c: A[:, i] for i, c in enumerate(cols)}
z = C["DEPT"]
HO, RO, GC = 1133.0, 2323.0, 1962.0
vp = 304800.0 / C["DTC"]
vs = 304800.0 / C["DTS"]
rho = C["RHOB"]
gr, cal, bs, drho = C["GR"], C["CALI"], C["BS"], C["DRHO"]
lith = C["FORCE_2020_LITHOFACIES_LITHOLOGY"]
hz = (z >= HO) & (z < RO)
ok = hz & np.isfinite(vp) & np.isfinite(vs) & np.isfinite(rho)
print(f"15/9-23: Hordaland {HO:.0f}-{RO:.0f} m MD; DTS present {np.nanmin(z[np.isfinite(vs)]):.0f}-{np.nanmax(z[np.isfinite(vs)]):.0f} m; samples in Hordaland with Vp, Vs and rho: {ok.sum()} ({ok.sum()*0.152:.0f} m)")
wash = cal - bs
print(f"caliper minus bit size in the Hordaland with Vs: median {np.nanmedian(wash[ok]):+.2f} in, p90 {np.nanpercentile(wash[ok], 90):+.2f} in; |DRHO| > 0.05 in {np.mean(np.abs(drho[ok]) > 0.05)*100:.0f} % of samples")
codes, counts = np.unique(lith[ok][np.isfinite(lith[ok])], return_counts=True)
print("lithology codes in that interval:", {int(c): int(n) for c, n in zip(codes, counts)})

print("\n  depth m MD | Vp    Vs    Vp/Vs  rho   GR   caliper-bit | Vs vs Nini trend | lithology")
for a in range(1480, 2340, 40):
    m = ok & (z >= a) & (z < a + 40)
    if m.sum() < 20:
        continue
    v, s_, r_ = np.nanmedian(vp[m]), np.nanmedian(vs[m]), np.nanmedian(rho[m])
    lc = np.unique(lith[m][np.isfinite(lith[m])], return_counts=True)
    dom = int(lc[0][np.argmax(lc[1])]) if len(lc[0]) else -1
    print(f"  {a:4d}-{a+40:4d} | {v:5.0f} {s_:5.0f} {v/s_:5.2f} {r_:5.3f} {np.nanmedian(gr[m]):5.0f} {np.nanmedian(wash[m]):+6.2f} | {100*(s_/vs_trend(v)-1):+6.1f} % | {dom}")

# local mudstone trend, good hole only
sh = ok & (lith == 65000) & (np.abs(wash) < 1.0) & (np.abs(drho) < 0.05)
print(f"\nshale samples with good hole: {sh.sum()} ({sh.sum()*0.152:.0f} m)")
# 1 m medians to reduce sample noise
zb = np.arange(HO, RO, 1.0)
P = []
for a in zb:
    m = sh & (z >= a) & (z < a + 1)
    if m.sum() >= 4:
        P.append((a, np.median(vp[m]), np.median(vs[m]), np.median(rho[m]), np.median(gr[m])))
P = np.array(P)
b1, b0 = np.polyfit(P[:, 1], P[:, 2], 1)
res = P[:, 2] - (b0 + b1 * P[:, 1])
rel = P[:, 2] / (b0 + b1 * P[:, 1]) - 1
print(f"local trend (1 m medians, n = {len(P)}): Vs = {b0:.0f} + {b1:.3f} Vp; scatter {np.std(res):.0f} m/s; relative scatter {100*np.std(rel):.1f} %")
print(f"   relative residual percentiles p1/p5/p50/p95/p99: " + " ".join(f"{100*np.percentile(rel, q):+.1f}" for q in (1, 5, 50, 95, 99)) + " %")
for L in (5, 10, 20, 30):
    # smooth over L m: what do layers of that thickness reach?
    k = np.ones(L) / L
    rs = np.convolve(rel, k, mode="valid")
    print(f"   {L:2d} m running mean of the relative residual: min {100*rs.min():+.1f} %, p1 {100*np.percentile(rs, 1):+.1f} %, sd {100*np.std(rs):.1f} %")
for lab, fn in (("Nini Lark trend", vs_trend), ("Castagna line", vs_castagna)):
    print(f"   15/9-23 shales against the {lab}: Vs median {100*(np.median(P[:, 2]/fn(P[:, 1])) - 1):+.1f} %")
gcm = P[:, 0] >= GC
print(f"predicted Green Clay part {GC:.0f}-{RO:.0f} m: Vp/Vs median {np.median(P[gcm, 1]/P[gcm, 2]):.2f}, residual {100*np.median(rel[gcm]):+.1f} %; above it: Vp/Vs {np.median(P[~gcm, 1]/P[~gcm, 2]):.2f}, residual {100*np.median(rel[~gcm]):+.1f} %")
# do low-density (soft) mudstones have low Vs relative to the trend?
lowr = P[:, 3] < np.percentile(P[:, 3], 20)
print(f"softest 20 % by density (rho < {np.percentile(P[:, 3], 20):.3f}): residual median {100*np.median(rel[lowr]):+.1f} %; others {100*np.median(rel[~lowr]):+.1f} %")


def vs_local(v):
    return b0 + b1 * v


# the Sleipner soft-layer base with the local trend
w, t, L = well()
vel = Vel(w)
g = stack_geometry(vel, 1.45)
a1, r1, a2, r2 = 2041.0, 1.885, 2212.0, 2.026
print(f"\n== base of the softest layer in 15/9-13 (Vp {a1:.0f}, rho {r1:.3f} over Vp {a2:.0f}, rho {r2:.3f}) at 1.45 s")
print(f"   local trend gives Vs {vs_local(a1):.0f} over {vs_local(a2):.0f} (Vp/Vs {a1/vs_local(a1):.2f} and {a2/vs_local(a2):.2f}); Nini trend {vs_trend(a1):.0f} over {vs_trend(a2):.0f}")
for lab, fn in (("local 15/9-23 trend", vs_local), ("Nini Lark trend", vs_trend)):
    base = None
    print(f"   -- background: {lab}")
    for red in (0, 5, 7, 10, 15, 20, 25):
        b1v = fn(a1) * (1 - red / 100)
        rr = {k: float(np.mean(zoep(a1, b1v, r1, a2, fn(a2), r2, g[k][0]))) for k in g}
        fn_ = rr["far"] / rr["nea"]
        base = fn_ if red == 0 else base
        print(f"      soft-layer Vs {red:2d} % below: R near {rr['nea']:+.4f} mid {rr['mid']:+.4f} far {rr['far']:+.4f}; far/near {fn_:.3f} ({100*(fn_/base-1):+.0f} % against no anomaly)")

# reference interfaces in 15/9-13 (20 ms blocks, |R0| > 0.015) with the local trend
w2, t2, L2 = well()
vp2, rho2 = L2["vp"], L2["rho_used"]
tb = np.arange(1270, 1960, 20)
bv = [(np.nanmedian(vp2[(t2 >= a) & (t2 < a + 20)]), np.nanmedian(rho2[(t2 >= a) & (t2 < a + 20)])) for a in tb]
ref = []
for i in range(len(tb) - 1):
    (p1, q1), (p2, q2) = bv[i], bv[i + 1]
    if abs(0.5 * np.log(p2 * q2 / (p1 * q1))) < 0.015:
        continue
    gg = stack_geometry(vel, (tb[i] + 20) / 1000.0)
    rr = {k: float(np.mean(zoep(p1, vs_local(p1), q1, p2, vs_local(p2), q2, gg[k][0]))) for k in gg}
    ref.append((rr["mid"] / rr["nea"], rr["far"] / rr["nea"]))
ref = np.array(ref)
print(f"\n   {len(ref)} reference interfaces in 15/9-13 with the local trend: mid/near median {np.median(ref[:,0]):.3f}, far/near median {np.median(ref[:,1]):.3f} (p10 {np.percentile(ref[:,1],10):.3f}, p90 {np.percentile(ref[:,1],90):.3f})")
