#!/usr/bin/env python3
"""Sleipner partial stacks: which offsets and angles does each stack hold?

The SEG-Y headers say nothing, so this is worked out from the data:
 - the full stack as a weighted sum of near, mid and far gives the live fold of each range against time,
 - the time-depth curve of 15/9-13 gives ray angles for any offset (exact ray tracing in a layered model),
 - plus polarity at the sea bed and the caliper in the Hordaland clay (is the density log to be trusted?).
Read-only.
"""
import json
import os
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026"
OUT = R + "/sleipner/avo"
DT = 0.002
S = {k: np.load(os.path.join(OUT, k + ".npy"), mmap_mode="r") for k in ("nea", "mid", "far", "ful")}
W = json.load(open(R + "/sleipner/cubes/sleipner94/wells.json"))
w = W["wells"][0]
meta = json.load(open(R + "/sleipner/cubes/sleipner94/meta.json"))
print("polarity notes in app files:", {k: v for k, v in w.items() if "pol" in k.lower()}, {k: v for k, v in meta.items() if "pol" in k.lower()}, {k: v for k, v in meta.get("segy", {}).items() if "pol" in k.lower()})

# ---------- fold weights against time
sub = {k: np.asarray(v[::3, ::3, :], np.float64) for k, v in S.items()}
print("\n== full ~ a*near + b*mid + c*far in 100 ms windows: live fold of mid and far relative to near")
print(" centre ms |    a      b      c    sum |  mid/near  far/near | R2")
rows = []
for c0 in range(450, 1960, 50):
    i0, i1 = (c0 - 50) // 2, (c0 + 50) // 2
    X = np.stack([sub[k][..., i0:i1].ravel() for k in ("nea", "mid", "far")], 1)
    y = sub["ful"][..., i0:i1].ravel()
    co, *_ = np.linalg.lstsq(X, y, rcond=None)
    r2 = 1 - np.sum((y - X @ co) ** 2) / np.sum(y ** 2)
    rows.append((c0, co, r2))
    print(f" {c0:8d} | {co[0]:6.3f} {co[1]:6.3f} {co[2]:6.3f} {co.sum():6.3f} | {co[1]/co[0]:8.2f} {co[2]/co[0]:9.2f} | {r2:.3f}")

# ---------- sea bed polarity
print("\n== mean near-stack trace around the sea bed (all traces), ms: amplitude")
m = np.asarray(S["nea"][:, :, 40:75], np.float64).mean((0, 1))
print(" ".join(f"{(40+i)*2}:{v:+.2f}" for i, v in enumerate(m)))
mf = np.asarray(S["ful"][:, :, 40:75], np.float64).mean((0, 1))
print("full: " + " ".join(f"{(40+i)*2}:{v:+.2f}" for i, v in enumerate(mf)))

# ---------- velocity model and ray angles
md = np.array(w["td"]["md"], float)
twt = np.array(w["td"]["twt"], float) / 1000.0
z = md - w["kb"]
ok = np.isfinite(z) & np.isfinite(twt)
z, twt = z[ok], twt[ok]
tl = np.arange(0.0, twt.max() + 1e-9, 0.02)
zl = np.interp(tl, np.r_[0.0, twt], np.r_[0.0, z])
vint = 2 * np.diff(zl) / np.diff(tl)
h = np.diff(zl)
print("\n== velocity model from the 15/9-13 time-depth curve")
print(" t0 s | depth m bsl | v_avg | v_rms | v_int just above")
for t0 in (0.6, 0.9, 1.1, 1.3, 1.45, 1.6, 1.85, 1.98):
    n = int(round(t0 / 0.02))
    vr = np.sqrt(np.sum(vint[:n] ** 2 * 0.02) / t0)
    print(f" {t0:4.2f} | {zl[n]:8.0f} | {2*zl[n]/t0:5.0f} | {vr:5.0f} | {vint[n-1]:5.0f}")


def angle(x, t0):
    """Incidence angle (deg) at a flat reflector at zero-offset time t0 for source-receiver offset x. Layered ray tracing."""
    n = int(round(t0 / 0.02))
    v, hh = vint[:n], h[:n]
    lo, hi = 0.0, 0.999999 / v.max()
    for _ in range(60):
        p = 0.5 * (lo + hi)
        xx = np.sum(2 * hh * p * v / np.sqrt(1 - (p * v) ** 2))
        if xx < x:
            lo = p
        else:
            hi = p
    return np.degrees(np.arcsin(min(1.0, p * v[-1])))


def offset_for_angle(th, t0):
    n = int(round(t0 / 0.02))
    v, hh = vint[:n], h[:n]
    p = np.sin(np.radians(th)) / v[-1]
    if p * v.max() >= 1:
        return np.nan
    return float(np.sum(2 * hh * p * v / np.sqrt(1 - (p * v) ** 2)))


print("\n== offset (m) that gives a certain incidence angle at the reflector")
print(" t0 s |   10     20     30     35     40     45     50 deg")
for t0 in (0.6, 0.9, 1.1, 1.3, 1.45, 1.6, 1.85):
    print(f" {t0:4.2f} | " + " ".join(f"{offset_for_angle(a, t0):6.0f}" for a in (10, 20, 30, 35, 40, 45, 50)))

# ---------- hypotheses for the offset ranges (width ratio near:mid:far from the late-time fold weights)
late = [r for r in rows if r[0] >= 1650]
rm = np.mean([r[1][1] / r[1][0] for r in late])
rf = np.mean([r[1][2] / r[1][0] for r in late])
print(f"\nlate-time fold ratios near:mid:far = 1 : {rm:.2f} : {rf:.2f}")
for cable, x0 in ((3000, 150), (3000, 250), (3600, 150)):
    wn = cable / (1 + rm + rf)
    x1, x2, x3 = x0 + wn, x0 + wn * (1 + rm), x0 + cable
    print(f"\n-- hypothesis: streamer {cable} m, nearest offset {x0} m -> near {x0:.0f}-{x1:.0f}, mid {x1:.0f}-{x2:.0f}, far {x2:.0f}-{x3:.0f} m")
    print("    t0 s | live far fold | outer live offset m | near deg       mid deg        far deg (to outer live offset)")
    for c0, co, r2 in rows:
        if c0 not in (900, 1000, 1100, 1200, 1300, 1450, 1600, 1750, 1850, 1950):
            continue
        t0 = c0 / 1000
        ff = np.clip((co[2] / co[0]) / rf, 0, 1)
        xo = x2 + ff * (x3 - x2)
        a = [angle(x, t0) for x in (x0, x1, x2, xo)]
        print(f"    {t0:4.2f} | {ff:12.2f} | {xo:18.0f} | {a[0]:4.1f}-{a[1]:4.1f}     {a[1]:4.1f}-{a[2]:4.1f}     {a[2]:4.1f}-{a[3]:4.1f}")

# ---------- caliper and logs from the LAS in the Hordaland clay
las = R + "/sleipner/raw/well_data/well_data/data/Wells_released_2011/Well 159_13/159-13 Logs.las"
cols, data = None, []
for line in open(las, errors="replace"):
    if line.startswith("~A"):
        cols = line.split()[1:]
        continue
    if cols is not None and line.strip():
        p = line.split()
        if len(p) == len(cols):
            data.append([float(v) for v in p])
A = np.array(data)
A[A == -999.25] = np.nan
col = {c: A[:, i] for i, c in enumerate(cols)}
print("\n== 15/9-13 LAS, 20 m bins: caliper (in), density, sonic, neutron, GR")
print("  MD m | CALI in | RHOB  | DT us/ft | DT-CKS | NPHI |  GR  | VSH")
for a in list(range(820, 1500, 20)) + list(range(1500, 2000, 40)):
    mm = (col["DEPTH"] >= a) & (col["DEPTH"] < a + 20)
    g = lambda k: np.nanmedian(col[k][mm]) if np.isfinite(col[k][mm]).any() else np.nan
    print(f" {a:5d} | {g('CALI'):7.2f} | {g('RHOB'):5.3f} | {g('DT'):8.1f} | {g('DT-CKS'):6.1f} | {g('NPHI'):4.2f} | {g('GR'):5.1f} | {g('VSH'):4.2f}")
