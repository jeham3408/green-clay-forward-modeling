"""Scratch: check-shot drift at Nini-4 (VSP one-way times against the integrated sonic log), as a rough size-of-effect check on
dispersion and hence Q, after the idea of Carter et al. 2020.  Not a Q measurement: borehole effects on the sonic log in soft clay,
picking conventions and anisotropy are not handled."""
import sys, numpy as np
las, vspf = sys.argv[1], sys.argv[2]
lines = open(las, errors="replace").read().split("\n")
c0 = next(i for i, l in enumerate(lines) if l.startswith("~C")); a0 = next(i for i, l in enumerate(lines) if l.startswith("~A"))
end = next(i for i in range(c0 + 1, len(lines)) if lines[i].startswith("~"))
names = [l.split(".")[0].strip() for l in lines[c0 + 1:end] if l.strip() and not l.startswith("#")]
arr = np.array([l.split() for l in lines[a0 + 1:] if l.strip()], float); arr[arr == -999.25] = np.nan
md = arr[:, 0]; dt = arr[:, names.index("DTCO")]; dt = np.where((dt > 40) & (dt < 250), dt, np.nan)
ok = np.isfinite(dt); z0, z1 = md[ok][0], md[ok][-1]
gaps = np.diff(md[ok]); print("sonic %.1f-%.1f m, %d samples, largest gap %.1f m" % (z0, z1, ok.sum(), gaps.max()))
zz = np.arange(z0, z1, 0.1524); slow = np.interp(zz, md[ok], dt[ok]) / 304800.0          # s/m
tson = np.concatenate([[0.0], np.cumsum(0.5 * (slow[1:] + slow[:-1]) * np.diff(zz))]) * 1000.0   # one-way ms from z0
rows = []
for line in open(vspf, errors="replace"):
    p = line.split()
    if len(p) == 9:
        try: v = [float(q) for q in p]
        except ValueError: continue
        rows.append((v[0], v[7]))
rows.sort(); z = np.array([r[0] for r in rows]); t = np.array([r[1] for r in rows])
m = (z >= z0 + 5) & (z <= z1 - 1); z, t = z[m], t[m]
ts = np.interp(z, zz, tson)
drift = (t - t[0]) - (ts - ts[0])                 # positive = the VSP (seismic frequency) is slower than the sonic
print("levels used: %d (%.0f-%.0f m). VSP one-way %.1f ms, integrated sonic %.1f ms, total drift %+.1f ms (%+.2f %%)" % (len(z), z[0], z[-1], t[-1] - t[0], ts[-1] - ts[0], drift[-1], 100 * drift[-1] / (t[-1] - t[0])))
print("depth  VSP-ms  sonic-ms  drift-ms")
for k in range(0, len(z), max(1, len(z) // 16)): print("%6.0f %7.1f %8.1f %+7.2f" % (z[k], t[k] - t[0], ts[k] - ts[0], drift[k]))
print("\ninterval: drift gradient and the Q a constant-Q medium would need (sonic 10 kHz against 50 Hz; range for 5-15 kHz against 30-80 Hz)")
for a, b, lab in ((820, 1220, "Lark L3"), (1220, 1540, "Lark L2"), (1540, 1700, "Lark L1"), (1700, 1760, "Horda"), (1760, 1870, "below Horda"), (820, 1700, "whole Lark"), (1490, 1700, "lower Lark, 10 m levels"), (1490, 1870, "all 10 m levels")):
    k = (z >= a) & (z <= b)
    if k.sum() < 3: print("  %-26s too few levels" % lab); continue
    dv = (t[k][-1] - t[k][0]); ds = (ts[k][-1] - ts[k][0])
    g = np.polyfit(t[k] - t[k][0], drift[k], 1)[0]                     # drift per unit VSP time (fraction), least squares over the levels
    resid = np.std(drift[k] - np.polyval(np.polyfit(t[k], drift[k], 1), t[k]))
    def q(frac, lo=np.log(10000 / 50.0)): return lo / (np.pi * frac) if frac > 0 else float("inf")
    print("  %-26s %2d levels, VSP %.1f ms, sonic %.1f ms, end-to-end %+.2f ms (%+.2f %%), fitted gradient %+.2f %%, scatter of levels %.2f ms -> Q about %s (range %s-%s)" % (
        lab, k.sum(), dv, ds, dv - ds, 100 * (dv - ds) / dv, 100 * g, resid,
        ("%.0f" % q(g)) if g > 0 else "none (negative drift)", ("%.0f" % q(g, np.log(5000 / 80.0))) if g > 0 else "-", ("%.0f" % q(g, np.log(15000 / 30.0))) if g > 0 else "-"))
