"""Scratch, after review: (1) is the low Vs of the Horda clay unusual, or what any softer mudrock on the local Vs-Vp trend has?
(2) how sharp is the top of the clay in the logs? (3) log quality across the boundary (caliper, density correction)."""
import sys, numpy as np
LAS = sys.argv[1]
lines = open(LAS, errors="replace").read().split("\n")
c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
g = lambda n: arr[:, names.index(n)]
md = g("DEPT"); dtc = np.where((g("DTCO") > 40) & (g("DTCO") < 250), g("DTCO"), np.nan); dts = np.where((g("DTSM") > 80) & (g("DTSM") < 800), g("DTSM"), np.nan)
rho = np.where((g("RHOZ") > 1.5) & (g("RHOZ") < 3.1), g("RHOZ"), np.nan); vp = 304800 / dtc; vs = 304800 / dts
ok = np.isfinite(vp) & np.isfinite(vs) & np.isfinite(rho)
print("=== 1. local Vs-Vp trend of the Lark mudstone (808-1700 m), and where the Horda clay (1702-1762 m) falls")
L = ok & (md >= 808) & (md < 1700); H = ok & (md >= 1702) & (md < 1762)
b, a = np.polyfit(vp[L], vs[L], 1); res = vs[L] - (a + b * vp[L])
print("  Lark: Vs = %.0f + %.3f Vp  (n=%d, scatter %.0f m/s)   [Castagna mudrock line: Vs = -1172 + 0.862 Vp]" % (a, b, L.sum(), res.std()))
for lab, m in (("Lark L3 808-1220", ok & (md >= 808) & (md < 1220)), ("Lark L2 1220-1545", ok & (md >= 1220) & (md < 1545)), ("Lark L1 1545-1700", ok & (md >= 1545) & (md < 1700)), ("Horda 1702-1762", H), ("sand 1766-1790", ok & (md >= 1766) & (md < 1790)), ("Sele+ 1792-1860", ok & (md >= 1792) & (md < 1860))):
    r = vs[m] - (a + b * vp[m]); rc = vs[m] - (-1172 + 0.862 * vp[m])
    print("  %-18s Vp %4.0f Vs %4.0f | Vs minus local Lark trend: median %+4.0f m/s (16-84 %%: %+4.0f..%+4.0f) | minus mudrock line: %+4.0f m/s" % (lab, np.median(vp[m]), np.median(vs[m]), np.median(r), *np.percentile(r, [16, 84]), np.median(rc)))
vpA, vsA, rA, vpB, vsB, rB = 2438.0, 1012.0, 2.295, 2215.0, 811.0, 2.171
def zoep(a1, b1, r1, a2, b2, r2, th):
    th = np.radians(np.asarray(th, float)); p = np.sin(th) / a1
    c1 = np.cos(th); c2 = np.sqrt(1 - (p * a2) ** 2 + 0j); d1 = np.sqrt(1 - (p * b1) ** 2 + 0j); d2 = np.sqrt(1 - (p * b2) ** 2 + 0j)
    aa = r2 * (1 - 2 * b2 ** 2 * p ** 2) - r1 * (1 - 2 * b1 ** 2 * p ** 2); bb = r2 * (1 - 2 * b2 ** 2 * p ** 2) + 2 * r1 * b1 ** 2 * p ** 2
    cc = r1 * (1 - 2 * b1 ** 2 * p ** 2) + 2 * r2 * b2 ** 2 * p ** 2; dd = 2 * (r2 * b2 ** 2 - r1 * b1 ** 2)
    E = bb * c1 / a1 + cc * c2 / a2; F = bb * d1 / b1 + cc * d2 / b2; G = aa - dd * c1 / a1 * d2 / b2; Hh = aa - dd * c2 / a2 * d1 / b1; D = E * F + G * Hh * p ** 2
    return (((bb * c1 / a1 - cc * c2 / a2) * F - (aa + dd * c1 / a1 * d2 / b2) * Hh * p ** 2) / D).real
def grad(a1, b1, r1, a2, b2, r2):
    am, bm, rm = (a1 + a2) / 2, (b1 + b2) / 2, (r1 + r2) / 2
    return 0.5 * (a2 - a1) / am - 2 * (bm / am) ** 2 * ((r2 - r1) / rm + 2 * (b2 - b1) / bm)
def fn(vs2):
    near = zoep(vpA, vsA, rA, vpB, vs2, rB, np.linspace(5, 15, 101)).mean(); far = zoep(vpA, vsA, rA, vpB, vs2, rB, np.linspace(25, 35, 101)).mean()
    full = zoep(vpA, vsA, rA, vpB, vs2, rB, np.linspace(0, 30, 301)).mean(); return grad(vpA, vsA, rA, vpB, vs2, rB), far / near, full
cases = (("measured Vs of the clay", vsB), ("clay on the local Lark Vs-Vp trend", a + b * vpB), ("clay with the same offset from the Lark trend as the layer above", a + b * vpB + (vsA - (a + b * vpA))),
         ("clay with the Vp/Vs of the layer above (2.43)", vpB / 2.43), ("clay on the mudrock line shifted through the layer above", vsA + 0.862 * (vpB - vpA)))
print("  null model for the clay's Vs -> Vs, gradient (mean Vs/Vp form), far/near (25-35 over 5-15 deg), full stack 0-30 deg")
for lab, v in cases:
    gq, f, full = fn(v); print("   %-66s Vs %4.0f  gradient %+.3f  far/near %.2f  full stack %+.4f" % (lab, v, gq, f, full))
print("\n=== 2. how sharp is the top? impedance in 1.5 m bins and the 10-90 % width of the step")
ai = vp * rho
up = np.nanmedian(ai[ok & (md >= 1670) & (md < 1696)]); dn = np.nanmedian(ai[ok & (md >= 1704) & (md < 1730)])
print("  median impedance 1670-1696 m: %.0f, 1704-1730 m: %.0f (step %.1f %%); 1600-1700 m: %.0f" % (up, dn, 100 * (dn - up) / up, np.nanmedian(ai[ok & (md >= 1600) & (md < 1700)])))
zz = np.arange(1686, 1716, 1.5)
print("  " + "  ".join("%.0f:%.0f" % (z, np.nanmedian(ai[ok & (md >= z) & (md < z + 1.5)])) for z in zz))
k = int(round(2.0 / 0.1524)) | 1; pad = np.pad(np.where(ok, ai, np.nan), k // 2, mode="edge"); sm = np.nanmedian(np.lib.stride_tricks.sliding_window_view(pad, k), axis=1)
w = (md >= 1690) & (md <= 1712); z90 = md[w][np.argmax(sm[w] < up - 0.1 * (up - dn))]; z10 = md[w][np.argmax(sm[w] < up - 0.9 * (up - dn))]
print("  2 m running median: falls through 10 %% of the step at %.1f m and 90 %% at %.1f m -> width %.1f m" % (z90, z10, z10 - z90))
print("\n=== 3. log quality: caliper (bit size from the well report is not in the LAS), density correction, by interval")
cal, hdra = g("HCAL"), g("HDRA")
for lab, lo, hi in (("Lark L1 1600-1700", 1600, 1700), ("Horda 1702-1762", 1702, 1762), ("sand 1766-1790", 1766, 1790), ("Sele+ 1792-1860", 1792, 1860)):
    m = (md >= lo) & (md < hi)
    print("  %-18s caliper median %.2f in (5-95 %%: %.2f-%.2f)  density correction median %+.3f (5-95 %%: %+.3f..%+.3f) g/cm3" % (lab, np.nanmedian(cal[m]), *np.nanpercentile(cal[m], [5, 95]), np.nanmedian(hdra[m]), *np.nanpercentile(hdra[m], [5, 95])))
print("  whole logged section 802-1865 m: caliper 5-95 %%: %.2f-%.2f in" % tuple(np.nanpercentile(cal[(md > 802) & (md < 1865)], [5, 95])))
