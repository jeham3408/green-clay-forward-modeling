"""Scratch: one 'what-if' at Nini-4. Same arithmetic as the app's tie (reflection coefficients from impedance on the 2 ms grid,
convolved with the zero-phase wavelet that has the spectrum of the seismic, positive number = increase in impedance, no shift).
The impedance of the Horda interval is replaced and the synthetic is compared with the real trace at the well."""
import json, sys, numpy as np
d = json.load(open(sys.argv[1])); w = d["wells"][0]; L = w["logs"]; tr = np.array(w["seis"]["trace"], float); dt_s, t0_s = d["dt_ms"], d["t0_ms"]
ai0 = np.array([np.nan if v is None else v for v in L["ai"]], float); n = L["n"]; tg = L["t0"] + L["dt"] * np.arange(n); md = np.array(L["md"], float)
tops = {t["key"]: t for t in w["tops"]}; t_hd, t_ro = tops["hd"]["twt"], tops["ro"]["twt"]
wv = np.array(w["wavelet"]["w"], float)
def ricker(f, dt):
    half = int(np.ceil(1500 / f / dt)); k = np.arange(-half, half + 1); a = np.pi * f * k * dt / 1000
    return (1 - 2 * a * a) * np.exp(-a * a)
def synth(ai, wav):
    rc = np.zeros(n); fin = np.isfinite(ai); k0, k1 = np.flatnonzero(fin)[[0, -1]]; prev = np.nan
    for k in range(k0, k1 + 1):
        a = ai[k] if np.isfinite(ai[k]) else prev
        if np.isfinite(prev) and np.isfinite(a) and a + prev > 0: rc[k - 1] = (a - prev) / (a + prev)
        if np.isfinite(a): prev = a
    return np.convolve(rc, wav, mode="same"), k0, k1
def seis_at(t):
    x = (t - t0_s) / dt_s; k = np.floor(x).astype(int); ok = (k >= 0) & (k < len(tr) - 1); out = np.full(len(x), np.nan)
    out[ok] = tr[k[ok]] + (tr[k[ok] + 1] - tr[k[ok]]) * (x[ok] - k[ok]); return out
def corr(syn, k0, k1, ta, tb):
    ka, kb = max(k0, int(np.ceil((ta - L["t0"]) / L["dt"]))), min(k1, int(np.floor((tb - L["t0"]) / L["dt"])))
    ks = np.arange(ka, kb + 1, 2); a = syn[ks]; b = seis_at(L["t0"] + ks * L["dt"]); m = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m], b[m])[0, 1]), int(m.sum())
inH = (md >= 1702.0) & (md < 1762.0) & np.isfinite(ai0)       # the Horda interval that has samples (top to the last Horda core)
above = (md >= 1600.0) & (md < 1700.0) & np.isfinite(ai0)
print("Horda on the time grid: %d samples, %.0f-%.0f ms; impedance median %.0f; mudstone above (1600-1700 m): %.0f" % (inH.sum(), tg[inH][0], tg[inH][-1], np.median(ai0[inH]), np.median(ai0[above])))
cases = {"as logged": ai0}
a = ai0.copy(); a[inH] = ai0[inH] * (np.median(ai0[above]) / np.median(ai0[inH])); cases["Horda scaled up to the impedance of the mudstone above (no soft layer)"] = a
a = ai0.copy(); a[inH] = np.median(ai0[above]); cases["Horda replaced by a constant = mudstone above"] = a
a = ai0.copy(); a[inH] = ai0[inH] * 0.90; cases["Horda 10 % softer than logged"] = a
a = ai0.copy(); a[inH] = ai0[inH] * 1.07; cases["Horda 7 % stiffer (half-way to the mudstone above)"] = a
for wname, wav in (("wavelet with the spectrum of the seismic (peak %.0f Hz)" % w["wavelet"]["peak_hz"], wv), ("Ricker 40 Hz", ricker(40, L["dt"]))):
    print("\n" + wname)
    for name, ai in cases.items():
        s, k0, k1 = synth(ai, wav)
        loc, nl = corr(s, k0, k1, t_hd - 60, t_ro + 40); allc, na = corr(s, k0, k1, 779, 1724)
        print("  %-72s local %4.0f-%4.0f ms: r = %+.2f (n=%d) | whole tie window: r = %+.2f" % (name, t_hd - 60, t_ro + 40, loc, nl, allc))
# how big is the change in the synthetic compared with the misfit that is already there?
s0, k0, k1 = synth(ai0, wv); s1, _, _ = synth(cases["Horda scaled up to the impedance of the mudstone above (no soft layer)"], wv)
ka, kb = int((t_hd - 60 - L["t0"]) / L["dt"]), int((t_ro + 40 - L["t0"]) / L["dt"]); ks = np.arange(ka, kb + 1, 2)
obs = seis_at(L["t0"] + ks * L["dt"]); sc = np.dot(obs, s0[ks]) / np.dot(s0[ks], s0[ks])
print("\nlocal window, synthetic scaled to the seismic: rms of seismic %.0f, rms of misfit (as logged) %.0f, rms of the change made by removing the soft layer %.0f" % (np.sqrt(np.mean(obs ** 2)), np.sqrt(np.mean((obs - sc * s0[ks]) ** 2)), np.sqrt(np.mean((sc * (s1[ks] - s0[ks])) ** 2))))
