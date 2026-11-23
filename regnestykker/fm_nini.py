"""Scratch calculation for the forward-modelling assessment (not part of the app). Reads the open Nini-4 LAS read-only."""
import sys, numpy as np
np.set_printoptions(linewidth=200, suppress=True)
LAS = sys.argv[1]
lines = open(LAS, errors="replace").read().split("\n")
c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
g = lambda n: arr[:, names.index(n)]
md = g("DEPT"); dtc = np.where((g("DTCO") > 40) & (g("DTCO") < 250), g("DTCO"), np.nan); dts = np.where((g("DTSM") > 80) & (g("DTSM") < 800), g("DTSM"), np.nan)
rho = np.where((g("RHOZ") > 1.5) & (g("RHOZ") < 3.1), g("RHOZ"), np.nan); gr = g("GR"); cal = g("HCAL")
vp = 304800.0 / dtc; vs = 304800.0 / dts

def stats(a, b, label):
    m = (md >= a) & (md < b) & np.isfinite(vp) & np.isfinite(vs) & np.isfinite(rho)
    q = lambda x: np.nanpercentile(x[m], [16, 50, 84])
    P, S, R = q(vp), q(vs), q(rho); VV = q(vp / vs); AI = q(vp * rho); SI = q(vs * rho)
    pr = q(0.5 * ((vp / vs) ** 2 - 2) / ((vp / vs) ** 2 - 1)); G = q(rho * vs ** 2 * 1e-6); M = q(rho * vp ** 2 * 1e-6)
    K = q((rho * vp ** 2 - 4.0 / 3 * rho * vs ** 2) * 1e-6); Gm, Km = rho * vs ** 2 * 1e-6, (rho * vp ** 2 - 4.0 / 3 * rho * vs ** 2) * 1e-6
    E = q(9 * Km * Gm / (3 * Km + Gm))
    print("%-34s n=%4d  Vp %4.0f [%4.0f-%4.0f]  Vs %4.0f [%4.0f-%4.0f]  rho %.3f [%.2f-%.2f]  Vp/Vs %.2f [%.2f-%.2f]  AI %4.0f  SI %4.0f  PR %.3f  G %.2f GPa  M %.1f  K %.1f  E %.2f GPa  GR %.0f  cal %.2f in"
          % (label, m.sum(), P[1], P[0], P[2], S[1], S[0], S[2], R[1], R[0], R[2], VV[1], VV[0], VV[2], AI[1], SI[1], pr[1], G[1], M[1], K[1], E[1], np.nanmedian(gr[m]), np.nanmedian(cal[m])))
    return dict(vp=P[1], vs=S[1], rho=R[1], vpvs=VV[1], n=int(m.sum()), E=E[1], G=G[1], pr=pr[1])

print("=== interval statistics (median [16th-84th percentile]); depths m below KB")
A1 = stats(1600, 1700, "A1 Lark L1 1600-1700"); B1 = stats(1700, 1760, "B1 Horda 1700-1760")
A2 = stats(1650, 1698, "A2 Lark L1 1650-1698"); B2 = stats(1704, 1756, "B2 Horda 1704-1756")
A3 = stats(1545, 1700, "A3 Lark L1 1545-1700"); B3 = stats(1702, 1762, "B3 Horda 1702-1762")
L2 = stats(1220, 1545, "Lark L2 1220-1545"); L3 = stats(808, 1220, "Lark L3 808-1220")
SD = stats(1766, 1790, "sand below Horda 1766-1790"); SE = stats(1792, 1860, "Sele and below 1792-1860")

def zoep(a1, b1, r1, a2, b2, r2, th):
    th = np.radians(np.asarray(th, float)); p = np.sin(th) / a1
    c1 = np.cos(th); c2 = np.sqrt(1 - (p * a2) ** 2 + 0j); d1 = np.sqrt(1 - (p * b1) ** 2 + 0j); d2 = np.sqrt(1 - (p * b2) ** 2 + 0j)
    a = r2 * (1 - 2 * b2 ** 2 * p ** 2) - r1 * (1 - 2 * b1 ** 2 * p ** 2); b = r2 * (1 - 2 * b2 ** 2 * p ** 2) + 2 * r1 * b1 ** 2 * p ** 2
    c = r1 * (1 - 2 * b1 ** 2 * p ** 2) + 2 * r2 * b2 ** 2 * p ** 2; d = 2 * (r2 * b2 ** 2 - r1 * b1 ** 2)
    E = b * c1 / a1 + c * c2 / a2; F = b * d1 / b1 + c * d2 / b2; G = a - d * c1 / a1 * d2 / b2; H = a - d * c2 / a2 * d1 / b1; D = E * F + G * H * p ** 2
    return (((b * c1 / a1 - c * c2 / a2) * F - (a + d * c1 / a1 * d2 / b2) * H * p ** 2) / D).real

def shuey(a1, b1, r1, a2, b2, r2, dd=0.0, de=0.0):
    a, b, r = (a1 + a2) / 2, (b1 + b2) / 2, (r1 + r2) / 2; da, db, dr = a2 - a1, b2 - b1, r2 - r1
    Z1, Z2 = r1 * a1, r2 * a2; G1, G2 = r1 * b1 ** 2, r2 * b2 ** 2
    R0 = (Z2 - Z1) / (Z2 + Z1)
    grad = 0.5 * (da / a - (2 * b / a) ** 2 * (G2 - G1) / ((G1 + G2) / 2) + dd)
    curv = 0.5 * (da / a + de)
    return R0, grad, curv

ang = [0, 10, 20, 30, 40]
def face(name, A, B, extra=""):
    R0, G, C = shuey(A["vp"], A["vs"], A["rho"], B["vp"], B["vs"], B["rho"])
    z = zoep(A["vp"], A["vs"], A["rho"], B["vp"], B["vs"], B["rho"], ang)
    print("%-52s R0 %+.4f  gradient %+.4f  curvature %+.4f | exact R(0,10,20,30,40 deg) = %s | R(30)/R(0) = %.2f %s" % (name, R0, G, C, " ".join("%+.4f" % v for v in z), z[3] / z[0], extra))
    return R0, G
print("\n=== interface: stiffer mudstone (A, above) over smectite-rich clay (B, below)")
for nm, A, B in (("A1 over B1 (1600-1700 over 1700-1760)", A1, B1), ("A2 over B2 (tight around the boundary)", A2, B2), ("A3 over B3 (whole L1 over whole Horda)", A3, B3)):
    face(nm, A, B)
# what part of the response is the low Vs?
Bn = dict(B1); Bn["vs"] = B1["vp"] / A1["vpvs"]
face("B1 with the Vp/Vs of A (same Vp, rho; 'normal' Vs)", A1, Bn)
Bv = dict(A1); Bv["vs"] = A1["vp"] / B1["vpvs"]
face("only Vs changes (A with the Vp/Vs of B)", A1, Bv)
Bp = dict(A1); Bp["vp"] = B1["vp"]; Bp["vs"] = B1["vp"] / A1["vpvs"]
face("only Vp changes (Vp/Vs kept)", A1, Bp)
Br = dict(A1); Br["rho"] = B1["rho"]
face("only density changes", A1, Br)
print("\n=== base of the clay: B over the units below")
face("B1 over the sand 1766-1790", B1, SD); face("B1 over Sele and below 1792-1860", B1, SE)
print("\n=== anisotropy as a confounder (Rueger, weak VTI): gradient and R(30 deg) for a contrast in Thomsen delta and epsilon (lower minus upper)")
for dd, de in ((0, 0), (0.05, 0.0), (-0.05, 0.0), (0.10, 0.0), (-0.10, 0.0), (0.0, 0.10), (0.0, -0.10), (0.05, 0.10), (-0.05, -0.10)):
    R0, G, C = shuey(A1["vp"], A1["vs"], A1["rho"], B1["vp"], B1["vs"], B1["rho"], dd, de)
    s2 = np.sin(np.radians(30)) ** 2; t2 = np.tan(np.radians(30)) ** 2
    print("  d_delta %+.2f d_epsilon %+.2f : gradient %+.4f  R(30) = %+.4f  (isotropic 3-term R(30) = %+.4f)" % (dd, de, G, R0 + G * s2 + C * s2 * t2, shuey(A1["vp"], A1["vs"], A1["rho"], B1["vp"], B1["vs"], B1["rho"])[0] + shuey(A1["vp"], A1["vs"], A1["rho"], B1["vp"], B1["vs"], B1["rho"])[1] * s2 + shuey(A1["vp"], A1["vs"], A1["rho"], B1["vp"], B1["vs"], B1["rho"])[2] * s2 * t2))

print("\n=== tuning: quarter wavelength in the clay (Vp of B1) and two-way time thickness")
for f in (20, 25, 30, 40, 50, 60, 75):
    print("  %2d Hz: wavelength %5.1f m, tuning thickness (1/4) %5.1f m, resolvable above about (1/8..1/4) %4.1f-%4.1f m" % (f, B1["vp"] / f, B1["vp"] / f / 4, B1["vp"] / f / 8, B1["vp"] / f / 4))
for h in (60, 100, 141, 200, 300, 430):
    print("  thickness %3d m: two-way time %5.1f ms, one-way %5.1f ms" % (h, 2000.0 * h / B1["vp"], 1000.0 * h / B1["vp"]))

print("\n=== attenuation: reflection from below a clay of thickness h; t* = two-way time / Q")
f = np.linspace(0.0, 150.0, 3001)
def ricker_spec(fp): return f ** 2 * np.exp(-(f / fp) ** 2)
def centroid(s): return float(np.sum(f * s) / np.sum(s))
for fp in (25, 30, 40):
    s0 = ricker_spec(fp); c0 = centroid(s0); var = float(np.sum((f - c0) ** 2 * s0) / np.sum(s0))
    print(" Ricker peak %d Hz: centroid %.1f Hz, std %.1f Hz" % (fp, c0, var ** 0.5))
    for h in (100, 200, 400):
        twt = 2.0 * h / B1["vp"]
        row = []
        for Q in (20, 30, 50, 80, 100, 150, 200):
            ts = twt / Q; s = s0 * np.exp(-np.pi * f * ts)
            row.append("Q%3d: t* %4.1f ms, A(peak) x%.2f, A(2xpeak) x%.2f, centroid %+.1f Hz" % (Q, ts * 1e3, np.exp(-np.pi * fp * ts), np.exp(-np.pi * 2 * fp * ts), centroid(s) - c0))
        print("   h %3d m (%.0f ms two-way): " % (h, twt * 1e3) + " | ".join(row))
print("\n=== the question that matters: a soft clay with Q_B against a stiffer mudstone with Q_A of the same thickness (difference only)")
for fp in (30,):
    s0 = ricker_spec(fp)
    for QA, QB in ((100, 50), (80, 40), (60, 30), (50, 25), (200, 50), (40, 20)):
        for h in (100, 200, 400):
            twt = 2.0 * h / B1["vp"]; dts_ = twt * (1.0 / QB - 1.0 / QA)
            ca = centroid(s0 * np.exp(-np.pi * f * twt / QA)); cb = centroid(s0 * np.exp(-np.pi * f * twt / QB))
            print("   Q %3d vs %3d, h %3d m: dt* %4.1f ms; amplitude at 30 Hz %.2f (%.1f dB), at 60 Hz %.2f (%.1f dB); centroid difference %.2f Hz" % (QA, QB, h, dts_ * 1e3, np.exp(-np.pi * 30 * dts_), 20 * np.log10(np.exp(-np.pi * 30 * dts_)), np.exp(-np.pi * 60 * dts_), 20 * np.log10(np.exp(-np.pi * 60 * dts_)), cb - ca))
