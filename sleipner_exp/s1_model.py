"""Step 1: 1D elastic model of the Sleipner section at well 15/9-13, from open data only.

Inputs (open):  FORCE 2020 LAS of 15/9-13 (DTC, RHOB, CALI, z_loc), 16/7-6 and 15/9-23 (measured DTS for the Vs trend);
                the checkshot-calibrated time-depth curve of 15/9-13 from the Sleipner 2019 reference dataset
                (CO2 DataShare), as stored in the hackathon app's wells.json.
Output:         model_1d.npz in SLEIPNER_EXP_ROOT (everything the later steps need; they do not read the logs).

Run:  SLEIPNER_EXP_ROOT=... FORCE_LAS_DIR=... SLEIPNER_WELLS_JSON=... python -B s1_model.py
"""
import glob
import json
import os

import numpy as np

ROOT = os.environ.get("SLEIPNER_EXP_ROOT", os.path.dirname(os.path.abspath(__file__)))
LAS_DIR = os.environ.get("FORCE_LAS_DIR", "/mnt/ssd/Data/force-hackathon-2026/force2020/las")
WELLS_JSON = os.environ.get("SLEIPNER_WELLS_JSON", "/mnt/ssd/Data/force-hackathon-2026/sleipner/cubes/sleipner94/wells.json")

ZMAX, DZ = 2300.0, 1.0
WATER_Z, VW, RW = 81.0, 1480.0, 1.03
SOFT_TOP, SOFT_BASE = 1358.0, 1390.0        # softest clay layer, m TVDSS (1383-1415 m MD, KB 25 m)
VS_FLOOR = 600.0
VP_CLIP = 3400.0


def read_las(name):
    f = glob.glob(os.path.join(LAS_DIR, "*", name))[0]
    names, rows, sec = [], [], None
    for line in open(f, errors="replace"):
        if line.startswith("~"):
            sec = line[1:2].upper()
            continue
        if line.startswith("#") or not line.strip():
            continue
        if sec == "C":
            names.append(line.split(".")[0].strip())
        elif sec == "A":
            rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    a[a == -999.25] = np.nan
    return {n: a[:, i] for i, n in enumerate(names)}


def block(z, v, ze, slow=False):
    out = np.full(len(ze) - 1, np.nan)
    ok = np.isfinite(z) & np.isfinite(v)
    idx = np.searchsorted(ze, z[ok], side="right") - 1
    g = (idx >= 0) & (idx < len(out))
    idx, vv = idx[g], v[ok][g]
    w = 1.0 / vv if slow else vv
    s = np.bincount(idx, weights=w, minlength=len(out))
    n = np.bincount(idx, minlength=len(out))
    m = n > 0
    out[m] = (n[m] / s[m]) if slow else (s[m] / n[m])
    return out


def fill(x):
    i = np.arange(len(x))
    ok = np.isfinite(x)
    return np.interp(i, i[ok], x[ok])


def medfilt(x, k):
    r = k // 2
    p = np.pad(x, r, mode="edge")
    return np.median(np.stack([p[i:i + len(x)] for i in range(k)]), 0)


# ------------------------------------------------------------------ Vs trend from measured shear logs (open FORCE wells)
trend = {}
pts = []
for nm, las, zlo, zhi in (("16/7-6", "16_7-6.las", 1081.0, 2005.0), ("15/9-23", "15_9-23.las", 1486.0, 2300.0)):
    L = read_las(las)
    z = -L["z_loc"]
    vp, vs = 304800.0 / L["DTC"], 304800.0 / L["DTS"]
    lith = L.get("FORCE_2020_LITHOFACIES_LITHOLOGY")
    m = np.isfinite(vp) & np.isfinite(vs) & (z >= zlo) & (z <= zhi) & (vp > 1800) & (vp < 2700) & (vp / vs > 1.6) & (vp / vs < 4.0)
    if lith is not None:
        m &= (lith == 65000)                                # shale
    A = np.c_[np.ones(m.sum()), vp[m]]
    c, *_ = np.linalg.lstsq(A, vs[m], rcond=None)
    res = vs[m] - A @ c
    trend[nm] = dict(a=float(c[0]), b=float(c[1]), n=int(m.sum()), sd=float(res.std()), z=[zlo, zhi],
                     vp_p10_p90=[float(np.percentile(vp[m], 10)), float(np.percentile(vp[m], 90))])
    pts.append((vp[m], vs[m]))
    # upper Hordaland of 16/7-6 separately (the level of the soft layer in 15/9-13: 240-270 m below top Hordaland)
    if nm == "16/7-6":
        mu = m & (z < 1081.0 + 450.0)
        A = np.c_[np.ones(mu.sum()), vp[mu]]
        c, *_ = np.linalg.lstsq(A, vs[mu], rcond=None)
        trend["16/7-6 upper Hordaland"] = dict(a=float(c[0]), b=float(c[1]), n=int(mu.sum()), sd=float((vs[mu] - A @ c).std()), z=[1081.0, 1531.0],
                                               vp_p10_p90=[float(np.percentile(vp[mu], 10)), float(np.percentile(vp[mu], 90))])
TR = trend["16/7-6 upper Hordaland"]                       # the trend used in the model
vs_trend = lambda v: TR["a"] + TR["b"] * v

# ------------------------------------------------------------------ 15/9-13 logs on 1 m cells
L = read_las("15_9-13.las")
z = -L["z_loc"]
kb = float(np.nanmedian((L["DEPT"] + L["z_loc"])[np.isfinite(L["z_loc"])][:200]))
ze = np.arange(0.0, ZMAX + DZ / 2, DZ)
zc = ze[:-1] + DZ / 2
vp_log = block(z, 304800.0 / L["DTC"], ze, slow=True)
rh_log = block(z, L["RHOB"], ze)
cali = block(z, L["CALI"], ze)
vp = fill(np.where((vp_log > 1450) & (vp_log < 6000), vp_log, np.nan))
vp = np.minimum(vp, VP_CLIP)
sea = zc < WATER_Z
vp[sea] = VW

# checkshot calibration: scale the slowness in 100 m segments so that the model's two-way time equals the
# checkshot-corrected time-depth curve of the well (Sleipner 2019 reference dataset) down to its end
W = json.load(open(WELLS_JSON))["wells"][0]
td_md, td_t = np.array(W["td"]["md"], float), np.array(W["td"]["twt"], float) / 1000.0
okt = np.isfinite(td_md) & np.isfinite(td_t)
td_z, td_t = td_md[okt] - W["kb"], td_t[okt]
keep = np.r_[True, np.diff(td_t) > 1e-6]
td_z, td_t = td_z[keep], td_t[keep]
twt_raw = np.r_[0.0, np.cumsum(2 * DZ / vp)]
z_end = float(td_z.max())
knots = np.r_[WATER_Z, np.arange(200.0, z_end, 100.0), z_end]
fac = np.ones(len(zc))
drift = []
for a, b in zip(knots[:-1], knots[1:]):
    t_td = np.interp(b, td_z, td_t) - np.interp(a, td_z, td_t)
    t_mod = np.interp(b, ze, twt_raw) - np.interp(a, ze, twt_raw)
    f = t_td / t_mod
    fac[(zc >= a) & (zc < b)] = f
    drift.append((float(a), float(b), float(f)))
vp_cal = vp / fac
vp_cal[sea] = VW
vp_cal = np.minimum(vp_cal, VP_CLIP)

# density: the log where the hole is good, else a Gardner-type trend fitted to the good samples
good = np.isfinite(rh_log) & (rh_log > 1.70) & (rh_log < 2.80) & np.isfinite(cali) & (cali < 21.0) & ~sea
gcoef = float(np.exp(np.mean(np.log(rh_log[good]) - 0.25 * np.log(vp_cal[good]))))
rho = np.where(good, rh_log, gcoef * vp_cal ** 0.25)
rho = medfilt(rho, 3)
# sea-bed sediment (no density log above about 470 m): Gardner trend, limited to 1.95-2.05 g/cm3 in the top 300 m
shallow = (~sea) & (zc < 470.0)
rho[shallow] = np.clip(gcoef * vp_cal[shallow] ** 0.25, 1.95, 2.05)
rho[sea] = RW

vs = np.maximum(vs_trend(vp_cal), VS_FLOOR)
vs = np.minimum(vs, vp_cal / 1.5)
vs[sea] = 0.0
soft = (zc >= SOFT_TOP) & (zc < SOFT_BASE)
twt = np.r_[0.0, np.cumsum(2 * DZ / vp_cal)]            # two-way time to cell edges, s

i_sb = int(np.argmax(~sea))
R_sb = (vp_cal[i_sb] * rho[i_sb] - VW * RW) / (vp_cal[i_sb] * rho[i_sb] + VW * RW)
info = dict(kb=kb, water_depth=WATER_Z, R_seabed=float(R_sb), vp_seabed=float(vp_cal[i_sb]), rho_seabed=float(rho[i_sb]),
            vs_trend_used="16/7-6 upper Hordaland", trends=trend, gardner_coef=gcoef,
            twt_seabed_ms=float(1000 * np.interp(WATER_Z, ze, twt)),
            twt_soft_ms=[float(1000 * np.interp(SOFT_TOP, ze, twt)), float(1000 * np.interp(SOFT_BASE, ze, twt))],
            twt_soft_td_ms=[float(1000 * np.interp(SOFT_TOP, td_z, td_t)), float(1000 * np.interp(SOFT_BASE, td_z, td_t))],
            soft_layer=dict(vp=float(1 / np.mean(1 / vp_cal[soft])), rho=float(rho[soft].mean()), vs_trend=float(vs[soft].mean())),
            below_soft_30m=dict(vp=float(1 / np.mean(1 / vp_cal[(zc >= SOFT_BASE) & (zc < SOFT_BASE + 30)])), rho=float(rho[(zc >= SOFT_BASE) & (zc < SOFT_BASE + 30)].mean())),
            above_soft_30m=dict(vp=float(1 / np.mean(1 / vp_cal[(zc >= SOFT_TOP - 30) & (zc < SOFT_TOP)])), rho=float(rho[(zc >= SOFT_TOP - 30) & (zc < SOFT_TOP)].mean())),
            vp_max=float(vp_cal.max()), vp_min_sed=float(vp_cal[~sea].min()), vs_min_sed=float(vs[~sea].min()), vs_max=float(vs.max()),
            cells_with_gardner_density=int((~good & ~sea).sum()), cells_vs_at_floor=int(((vs <= VS_FLOOR + 1e-6) & ~sea).sum()),
            drift_factor_minmax=[float(min(d[2] for d in drift)), float(max(d[2] for d in drift))], td_end_z=z_end,
            twt_at_zmax_ms=float(1000 * twt[-1]), twt_drift_before_cal_ms_at_soft_top=float(1000 * (np.interp(SOFT_TOP, ze, twt_raw) - np.interp(SOFT_TOP, td_z, td_t))))
np.savez_compressed(os.path.join(ROOT, "model_1d.npz"), z=zc, ze=ze, vp=vp_cal, vs=vs, rho=rho, twt=twt, soft=soft,
                    vp_uncal=vp, good_rho=good, info=json.dumps(info))
json.dump(info, open(os.path.join(ROOT, "model_info.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
for a in range(0, 2300, 100):
    m = (zc >= a) & (zc < a + 100)
    print("%4d-%4d m  twt %6.0f ms  vp %4.0f  vs %4.0f  rho %.3f  gardner %3d %%" % (a, a + 100, 1000 * np.interp(a, ze, twt), 1 / np.mean(1 / vp_cal[m]), vs[m].mean(), rho[m].mean(), 100 * (~good[m]).mean()))
