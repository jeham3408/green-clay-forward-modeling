#!/usr/bin/env python3
"""FORCE 2020 wells near Sleipner with measured shear sonic (DTS) in the upper Hordaland Group, at the level of the
soft layer in 15/9-13 (240-272 m below the top of the Hordaland Group). Is there a layer with low Vs there?
Read-only. Run:  ssh linux-pc '<geoscience python> -I -' < fm_sl_dtsnear.py
"""
import glob
import json
import math
import os
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026"
LASDIR = R + "/force2020/las/Force_2020_all_wells_train_test_blind_hidden_final"
G = json.load(open(R + "/natures-plug/data/greenclay.json"))
ex = G["ex"]
TOPS = {ex["n"][i]: (ex["la"][i], ex["lo"][i], ex["ho"][i], ex["ro"][i]) for i in range(len(ex["n"]))}
LAT0, LON0 = 58.37388, 1.93413     # 15/9-13


def dist_km(la, lo):
    return math.hypot((la - LAT0) * 111.2, (lo - LON0) * 111.2 * math.cos(math.radians(LAT0)))


def read_las(p):
    cols, rows, mode, well = [], [], None, None
    for line in open(p, errors="replace"):
        s = line.strip()
        if s.startswith("~C"):
            mode = "c"; continue
        if s.startswith("~A"):
            mode = "a"; continue
        if s.startswith("~"):
            mode = "w" if s.startswith("~W") else None; continue
        if mode == "w" and s.startswith("WELL"):
            well = s.split(".", 1)[1].split(":")[0].strip()
        elif mode == "c" and s and not s.startswith("#"):
            cols.append(s.split(".")[0].strip())
        elif mode == "a" and s and not s.startswith("#"):
            q = s.split()
            if len(q) == len(cols):
                rows.append([float(v) for v in q])
    A = np.array(rows) if rows else np.zeros((0, len(cols)))
    if A.size:
        A[A == -999.25] = np.nan
    return well, {c: A[:, i] for i, c in enumerate(cols)}


rows = []
for p in sorted(glob.glob(LASDIR + "/*.las")):
    well, C = read_las(p)
    if well is None or "DTS" not in C or "DTC" not in C:
        continue
    name = well.split()[0]
    key = next((k for k in TOPS if k == name or k.replace(" ", "") == name.replace(" ", "")), None)
    if key is None:
        continue
    la, lo, ho, ro = TOPS[key]
    d = dist_km(la, lo)
    if d > 60:
        continue
    z, dts, dtc = C["DEPT"], C["DTS"], C["DTC"]
    rho = C.get("RHOB", np.full_like(z, np.nan))
    ok = np.isfinite(dts) & np.isfinite(dtc) & (dts > 0) & (dtc > 0)
    if not ok.any():
        continue
    zs = z[ok]
    a, b = ho + 240, ho + 272          # level of the soft layer below the Hordaland top in 15/9-13
    lev = ok & (z >= ho + 150) & (z < ho + 400)
    cover = np.mean(ok[(z >= a) & (z < b)]) if np.any((z >= a) & (z < b)) else 0
    rows.append((d, key, ho, ro, zs.min(), zs.max(), cover, p, lev.sum()))
rows.sort()
print("FORCE wells within 60 km of 15/9-13 with DTS and DTC (distance, well, Hordaland top-base, DTS range, coverage of top+240..272 m):")
for r in rows:
    print(f"   {r[0]:5.1f} km  {r[1]:12s}  Hordaland {r[2]:6.0f}-{r[3]:6.0f} m  DTS {r[4]:6.0f}-{r[5]:6.0f} m  level covered {100*r[6]:3.0f} %")

print("\nwells that cover the level: Vs, Vp, Vp/Vs and density in 20 m bins from top+150 to top+400 m, and Vs relative to a")
print("Vs-Vp line fitted to that well's upper Hordaland shale (GR above the well's median)")
for d, key, ho, ro, zmin, zmax, cover, p, n in rows:
    if cover < 0.5:
        continue
    _, C = read_las(p)
    z = C["DEPT"]
    vp = 304800.0 / C["DTC"]
    vs = 304800.0 / C["DTS"]
    rho = C.get("RHOB", np.full_like(z, np.nan))
    gr = C.get("GR", np.full_like(z, np.nan))
    m = (z >= ho + 100) & (z < ho + 450) & np.isfinite(vp) & np.isfinite(vs)
    sh = m & (gr > np.nanmedian(gr[m]))
    b1, b0 = np.polyfit(vp[sh], vs[sh], 1)
    print(f"\n   {key} ({d:.1f} km), Hordaland top {ho:.0f} m; shale line Vs = {b0:.0f} + {b1:.3f} Vp")
    print("     depth below top | Vp    Vs    Vp/Vs  rho    GR  | Vs vs line")
    for off in range(150, 400, 20):
        mm = m & (z >= ho + off) & (z < ho + off + 20)
        if mm.sum() < 20:
            continue
        v, s_, r_ = np.nanmedian(vp[mm]), np.nanmedian(vs[mm]), np.nanmedian(rho[mm])
        res = np.nanmedian(vs[mm] / (b0 + b1 * vp[mm]) - 1)
        flag = "  <-- level of the soft layer in 15/9-13" if ho + off <= ho + 272 and ho + off + 20 > ho + 240 else ""
        print(f"     {off:4d}-{off+20:4d} m      | {v:5.0f} {s_:5.0f} {v/s_:5.2f} {r_:5.3f} {np.nanmedian(gr[mm]):5.0f} | {100*res:+5.1f} %{flag}")
    k = np.ones(66) / 66            # 10 m running mean at 0.152 m sampling
    rel = vs[m] / (b0 + b1 * vp[m]) - 1
    rel = rel[np.isfinite(rel)]
    if len(rel) > 200:
        rs = np.convolve(rel, k, mode="valid")
        print(f"     10 m running mean of Vs relative to the line, top+100..450 m: min {100*rs.min():+.1f} %, p1 {100*np.percentile(rs, 1):+.1f} %")
