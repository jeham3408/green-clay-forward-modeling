"""Scratch: full-stack averaging of the AVA response, a wedge (tuning) model of the soft clay, and the VSP numbers."""
import numpy as np
# Nini-4 medians (from fm_nini.py): A = Lark L1 1600-1700 m, B = Horda 1700-1760 m, C = Sele and below 1792-1860 m
A = dict(vp=2438.0, vs=1012.0, rho=2.295); B = dict(vp=2215.0, vs=811.0, rho=2.171); C = dict(vp=2497.0, vs=1087.0, rho=2.280)
def zoep(m1, m2, th):
    a1, b1, r1, a2, b2, r2 = m1["vp"], m1["vs"], m1["rho"], m2["vp"], m2["vs"], m2["rho"]
    th = np.radians(np.asarray(th, float)); p = np.sin(th) / a1
    c1 = np.cos(th); c2 = np.sqrt(1 - (p * a2) ** 2 + 0j); d1 = np.sqrt(1 - (p * b1) ** 2 + 0j); d2 = np.sqrt(1 - (p * b2) ** 2 + 0j)
    a = r2 * (1 - 2 * b2 ** 2 * p ** 2) - r1 * (1 - 2 * b1 ** 2 * p ** 2); b = r2 * (1 - 2 * b2 ** 2 * p ** 2) + 2 * r1 * b1 ** 2 * p ** 2
    c = r1 * (1 - 2 * b1 ** 2 * p ** 2) + 2 * r2 * b2 ** 2 * p ** 2; d = 2 * (r2 * b2 ** 2 - r1 * b1 ** 2)
    E = b * c1 / a1 + c * c2 / a2; F = b * d1 / b1 + c * d2 / b2; G = a - d * c1 / a1 * d2 / b2; H = a - d * c2 / a2 * d1 / b1; D = E * F + G * H * p ** 2
    return (((b * c1 / a1 - c * c2 / a2) * F - (a + d * c1 / a1 * d2 / b2) * H * p ** 2) / D).real
Bn = dict(B); Bn["vs"] = B["vp"] / 2.43          # the Vp/Vs of the layer above as the median of the ratio, as in fm_nini.py and fm_null.py
print("=== what a FULL STACK sees of the top of the clay: mean of R over the angle range (equal weight per degree)")
for lo, hi in ((0, 25), (0, 30), (0, 35), (5, 35)):
    th = np.linspace(lo, hi, 301)
    print("  %2d-%2d deg: measured Vs %+.4f | same clay with the Vp/Vs of the mudstone above %+.4f | difference %.1f %% of the amplitude" % (lo, hi, zoep(A, B, th).mean(), zoep(A, Bn, th).mean(), 100 * (zoep(A, B, th).mean() / zoep(A, Bn, th).mean() - 1)))
print("=== what PARTIAL STACKS would see: near 5-15, mid 15-25, far 25-35 deg")
for nm, m2 in (("measured Vs", B), ("Vp/Vs of the mudstone above", Bn)):
    v = [zoep(A, m2, np.linspace(a, b, 101)).mean() for a, b in ((5, 15), (15, 25), (25, 35))]
    print("  %-28s near %+.4f mid %+.4f far %+.4f | far/near %.2f" % (nm, v[0], v[1], v[2], v[2] / v[0]))
print("=== wedge: soft clay (B) between the mudstone above (A) and the unit below (C); peak-to-trough amplitude of the composite, relative to a thick layer")
def ricker(fp, dt):
    half = int(np.ceil(2.0 / fp / dt)); t = np.arange(-half, half + 1) * dt; a = (np.pi * fp * t) ** 2
    return (1 - 2 * a) * np.exp(-a)
dt = 0.0005; rt = (B["vp"] * B["rho"] - A["vp"] * A["rho"]) / (B["vp"] * B["rho"] + A["vp"] * A["rho"]); rb = (C["vp"] * C["rho"] - B["vp"] * B["rho"]) / (C["vp"] * C["rho"] + B["vp"] * B["rho"])
print("  R top %+.4f, R base %+.4f" % (rt, rb))
for fp in (25, 40, 60):
    w = ricker(fp, dt); rows = []
    ref = None
    for h in (300, 100, 60, 40, 30, 22, 15, 10, 5):
        n = int(round(2.0 * h / B["vp"] / dt)); r = np.zeros(4000); r[1000] = rt; r[1000 + n] = rb
        s = np.convolve(r, w, mode="same"); amp = s.max() - s.min()
        if ref is None: ref = amp
        rows.append("%3d m: %.2f" % (h, amp / ref))
    print("  Ricker %2d Hz (tuning near %.0f m): " % (fp, B["vp"] / (4 * fp * 1.3)) + " | ".join(rows))
print("=== VSP at Nini-4: one-way times (ms below sea level, from the CGG listing) and what Q would do to the spectral ratio between two levels")
lv = {1490: 707.1, 1700: 795.8, 1760: 822.7, 1880: 869.1, 800: 415.8, 320: 165.9}
for z1, z2, lab in ((320, 1880, "whole profile"), (800, 1700, "Lark Formation"), (1490, 1700, "lower Lark, 10 m levels"), (1700, 1760, "Horda Formation"), (1490, 1880, "all 40 close levels")):
    dts = (lv[z2] - lv[z1]) / 1000.0
    print("  %4d-%4d m (%s): one-way %.1f ms; slope of ln(spectral ratio) for Q 30 / 50 / 100: %.4f / %.4f / %.4f per Hz -> change over 10-90 Hz: %.2f / %.2f / %.2f (ln units), %.1f / %.1f / %.1f dB"
          % (z1, z2, lab, dts * 1e3, np.pi * dts / 30, np.pi * dts / 50, np.pi * dts / 100, 80 * np.pi * dts / 30, 80 * np.pi * dts / 50, 80 * np.pi * dts / 100, 8.686 * 80 * np.pi * dts / 30, 8.686 * 80 * np.pi * dts / 50, 8.686 * 80 * np.pi * dts / 100))
