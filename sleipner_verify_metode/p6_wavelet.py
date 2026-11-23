"""Point 1, robustness: how much do the baseline and the low-Vs scenarios depend on the wavelet and on NMO stretch?

(A) check of my NMO version against the experiment's primaries-only Zoeppritz stacks (C_trend_main in stacks.npz), with
    the same ghosts (6 m source, 8 m receiver) put into my wavelet;
(B) baseline zone/reference against wavelet frequency, with and without NMO stretch, with and without ghosts;
(C) the key low-Vs scenarios with the ghosted wavelet and NMO (the set-up closest to the experiment);
(D) other Vs laws for the whole section (bed-scale Vs following the mudrock slope instead of the flat local trend);
(E) low Vs in the soft beds everywhere, not only in the zone.
Writes p6_results.json.
"""
import json, os
import numpy as np
import slp, vz
HERE = os.path.dirname(os.path.abspath(__file__))
m = vz.load_model()
z = m["z"]; n = len(z)
F = vz.Forward(m)
Z = np.load(os.path.join(HERE, "stacks.npz"))
OUT = {}

def ricker(f, nn=401):
    t = (np.arange(nn) - nn // 2) * vz.DT
    a = (np.pi * f * t) ** 2
    return (1 - 2 * a) * np.exp(-a)

def ghosted(w, zs=6.0, zr=8.0, c=1480.0):
    ts, tr = 2 * zs / c / vz.DT, 2 * zr / c / vz.DT           # delays in samples (vertical incidence)
    nn = len(w)
    f = np.fft.rfftfreq(nn)
    W = np.fft.rfft(w)
    sh = lambda s: np.exp(-2j * np.pi * f * s)
    # centred on the ghost-free arrival time shifted by half the total delay, as in the experiment's synthetic
    H = (1 - sh(ts)) * (1 - sh(tr)) * sh(-(ts + tr) / 2)
    return np.fft.irfft(W * H, nn)

def peak_freq(w):
    f = np.fft.rfftfreq(4096, vz.DT)
    return float(f[np.argmax(np.abs(np.fft.rfft(w, 4096)))])

def setw(w):
    F.wavelet = lambda f=None, n=301: w

def run(fac=None, stretch=True, vs=None, events=None):
    if vs is not None:
        R = F._refl(vs, np.arange(F.n - 1))
    else:
        R = F.refl(fac)
    S = F.stacks(R, stretch=stretch)
    r = vz.rms_measure(S)
    e = vz.event_measure(S, events)
    return S, dict(rms_fn=r["fn"], ev_fn=e["fn"], pk_fn=e["pk_fn"], rms_mn=r["mn"], ev_mn=e["mn"], n=[e["n_zone"], e["n_ref"]], zone_values=e["fn_zone_values"])

# ---- (A)
wg = ghosted(ricker(25.0))
OUT["peak_freq_hz"] = dict(ricker25=peak_freq(ricker(25.0)), ricker25_ghosted=peak_freq(wg), ricker35=peak_freq(ricker(35.0)), ricker35_ghosted=peak_freq(ghosted(ricker(35.0))))
setw(wg)
S, r = run(None, True)
chk = {}
w19 = (vz.T_OUT >= 1.1) & (vz.T_OUT <= 1.9)
for k in ("nea", "mid", "far"):
    c = Z["C_trend_main__" + k].astype(float)
    best = max(((float(np.corrcoef(c[w19], np.roll(S[k], L)[w19])[0, 1]), L) for L in range(-20, 21)), key=lambda v: v[0])
    chk[k] = dict(best_corr=best[0], lag_ms=best[1])
OUT["A_my_nmo_ghosted_vs_experiment_C_trend"] = dict(corr=chk, my_zone_over_ref=r,
    experiment_C=dict(rms_fn=0.96, ev_fn=0.89, pk_fn=0.86))
# their C stack measured with my functions (same events as theirs = picked on their G; here: picked on their C full stack)
Sc = {k: Z["C_trend_main__" + k].astype(float) for k in ("nea", "mid", "far")}
Sc["ful"] = sum(slp.W_FULL[k] * Sc[k] for k in ("nea", "mid", "far"))
rc = vz.rms_measure(Sc); ec = vz.event_measure(Sc, None)
OUT["A_experiment_C_measured_here_own_events"] = dict(rms_fn=rc["fn"], ev_fn=ec["fn"], pk_fn=ec["pk_fn"], n=[ec["n_zone"], ec["n_ref"]])

# ---- (B)
B = {}
for nm, w in [("ricker%d" % f, ricker(float(f))) for f in (15, 20, 25, 30, 35, 40, 50)] + [("ricker%d_ghosted" % f, ghosted(ricker(float(f)))) for f in (15, 20, 25, 30, 35)]:
    setw(w)
    row = dict(peak_hz=peak_freq(w))
    for st in (False, True):
        _, r = run(None, st)
        row["nmo" if st else "nmofree"] = {q: r[q] for q in ("rms_fn", "ev_fn", "pk_fn", "rms_mn", "ev_mn", "n")}
    B[nm] = row
OUT["B_baseline_by_wavelet"] = B

# ---- (C) scenarios
SOFT = (1358.0, 1390.0)
zt, zb = [float(np.interp(t, m["twt"], m["ze"])) for t in vz.ZONE]
ai = m["vp"] * m["rho"]
def runmean(x, L):
    k = np.ones(L) / L
    return np.convolve(np.pad(x, L // 2, mode="edge"), k, mode="valid")[:len(x)]
def box(a, b, X):
    f = np.ones(n); f[(z >= a) & (z < b)] = 1 - X; return f
def softbeds(a, b, X, L=61):
    d = runmean(ai, L) - ai
    f = np.ones(n); f[(z >= a) & (z < b) & (d > 0)] = 1 - X; return f
SCN = {"baseline": None}
for X in (0.07, 0.15, 0.20):
    p = "%02d" % round(100 * X)
    SCN["box32_softlayer_" + p] = box(SOFT[0], SOFT[1], X)
    SCN["box_wholezone_" + p] = box(zt, zb, X)
    SCN["box150_centred_" + p] = box(1374 - 75, 1374 + 75, X)
    SCN["softbeds_zone_L61_" + p] = softbeds(zt, zb, X)
    SCN["softbeds_zone_L31_" + p] = softbeds(zt, zb, X, 31)
    SCN["softbeds_everywhere_1000_2000m_L61_" + p] = softbeds(1000.0, 2000.0, X)
Cc = {}
for wn, w in (("ricker25", ricker(25.0)), ("ricker25_ghosted", wg), ("ricker35", ricker(35.0)), ("ricker18", ricker(18.0))):
    setw(w)
    for st in (False, True):
        tag = "%s|%s" % (wn, "nmo" if st else "nmofree")
        S0, r0 = run(None, st)
        ev0 = slp.pick_events(S0["ful"])
        row = {}
        for nm, fac in SCN.items():
            _, r = run(fac, st, events=ev0)
            _, r2 = run(fac, st, events=None)
            row[nm] = dict(rms_fn=r["rms_fn"], ev_fixed_fn=r["ev_fn"], ev_own_fn=r2["ev_fn"], pk_fixed_fn=r["pk_fn"], rms_mn=r["rms_mn"], ev_fixed_mn=r["ev_mn"], zone_values=r["zone_values"], n_own=r2["n"])
        Cc[tag] = row
OUT["C_scenarios_by_wavelet"] = Cc

# ---- (D) other Vs laws for the whole section
vs0 = m["vs"].copy()
vp = m["vp"]
sed = vs0 > 0
laws = {}
vps = runmean(vp, 101)
for slope in (0.343, 0.6, 0.77, 0.862):
    v = vs0.copy()
    v[sed] = (46.24 + 0.3433 * vps[sed]) + slope * (vp[sed] - vps[sed])       # background as the experiment, bed scale with another slope
    laws["bedscale_slope_%.3f" % slope] = np.clip(v, 300.0, vp / 1.5)
v = vs0.copy(); v[sed] = np.clip(0.862 * vp[sed] - 1172.0, 300.0, None); laws["castagna_mudrock_line"] = v
D = {}
for wn, w in (("ricker25", ricker(25.0)), ("ricker25_ghosted", wg)):
    setw(w)
    for st in (False, True):
        tag = "%s|%s" % (wn, "nmo" if st else "nmofree")
        row = {}
        for nm, v in laws.items():
            _, r = run(vs=v, stretch=st)
            c = (z >= SOFT[0]) & (z < SOFT[1])
            row[nm] = dict(rms_fn=r["rms_fn"], ev_fn=r["ev_fn"], pk_fn=r["pk_fn"], rms_mn=r["rms_mn"], ev_mn=r["ev_mn"], n=r["n"], vs_soft_layer=float(v[c].mean()), vpvs_soft_layer=float((vp[c] / v[c]).mean()))
        D[tag] = row
OUT["D_vs_laws_whole_section"] = D
json.dump(OUT, open(os.path.join(HERE, "p6_results.json"), "w"), indent=1, default=float)

f = lambda v: "  nan" if v is None else "%5.2f" % v
print("peak frequencies:", OUT["peak_freq_hz"])
print("(A)", json.dumps(OUT["A_my_nmo_ghosted_vs_experiment_C_trend"], default=lambda v: round(float(v), 3)))
print("    experiment's C measured here:", json.dumps(OUT["A_experiment_C_measured_here_own_events"], default=lambda v: round(float(v), 3)))
print("(B) baseline by wavelet: peak Hz | nmofree RMS ev pk | nmo RMS ev pk | mid/near nmo RMS ev")
for nm, r in B.items():
    print("  %-18s %5.1f | %s %s %s | %s %s %s | %s %s  n %s" % (nm, r["peak_hz"], f(r["nmofree"]["rms_fn"]), f(r["nmofree"]["ev_fn"]), f(r["nmofree"]["pk_fn"]), f(r["nmo"]["rms_fn"]), f(r["nmo"]["ev_fn"]), f(r["nmo"]["pk_fn"]), f(r["nmo"]["rms_mn"]), f(r["nmo"]["ev_mn"]), r["nmo"]["n"]))
print("(C) scenarios: far/near RMS evFix evOwn pkFix | mid/near RMS evFix")
for tag, row in Cc.items():
    print(" ", tag)
    for nm, r in row.items():
        print("    %-40s %s %s %s %s | %s %s  %s" % (nm, f(r["rms_fn"]), f(r["ev_fixed_fn"]), f(r["ev_own_fn"]), f(r["pk_fixed_fn"]), f(r["rms_mn"]), f(r["ev_fixed_mn"]), r["zone_values"]))
print("(D) Vs laws: far/near RMS ev pk | mid/near RMS ev | Vs and Vp/Vs in the soft layer")
for tag, row in D.items():
    print(" ", tag)
    for nm, r in row.items():
        print("    %-26s %s %s %s | %s %s | %4.0f %.2f  n %s" % (nm, f(r["rms_fn"]), f(r["ev_fn"]), f(r["pk_fn"]), f(r["rms_mn"]), f(r["ev_mn"]), r["vs_soft_layer"], r["vpvs_soft_layer"], r["n"]))
