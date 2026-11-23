
# ---- fm_sl_robust.py: is the far-offset dimming at the soft clay layer real, or a measuring artefact? ----
# (Multiples are tested in fm_sl_mult.py and removed in fm_sl_demult.py.)
# Run:  cat sl_common.py sl_events_lib.py fm_sl_robust.py | ssh linux-pc '<geoscience python> -I -'
import sys

S = stacks()
S3 = {k: S[k] for k in ("nea", "mid", "far")}
KEY = [1276, 1316, 1332, 1386, 1456, 1536, 1746]

# ---------- 0. shallow section: are the deep events at twice the time of strong shallow reflectors?
mt = {k: np.asarray(S[k][::2, ::2, :], np.float64).mean((0, 1)) for k in KEYS}
print("== strongest flat events 500-1000 ms on the cube-wide mean near trace (possible multiple generators)")
r0 = np.sqrt(np.mean(mt["nea"][250:500] ** 2))
for i in range(251, 499):
    v = mt["nea"][i]
    if abs(v) > 1.3 * r0 and abs(v) >= abs(mt["nea"][i - 1]) and abs(v) >= abs(mt["nea"][i + 1]):
        print(f"   {i*2:5d} ms  near {v:+.3f}  mid {mt['mid'][i]:+.3f}   (twice the time: {i*4} ms; plus one water layer: {i*2+110} ms)")

# ---------- 1. all events, lower threshold: is the dimming confined to the soft layer?
r1 = np.sqrt(np.mean(mt["ful"][540:985] ** 2))
ev = []
for i in range(541, 951):          # stop at 1900 ms: the tracking window must stay inside the 2 s record
    v = mt["ful"][i]
    if abs(v) > 0.8 * r1 and abs(v) >= abs(mt["ful"][i - 1]) and abs(v) >= abs(mt["ful"][i + 1]):
        if not ev or i * 2 - ev[-1] >= 12:
            ev.append(i * 2)
print(f"\n== {len(ev)} events from the cube-wide mean full-stack trace (threshold 0.8 x rms)")
print("  time | pol | kept % | far/near (95 %)       mid/near (95 %)      | no-alignment peak far/near  mid/near | lag far-near ms | mean wavelet amp near mid far")
ALL = {}
HOR = {}
for T0 in ev:
    hz = track(S["ful"], T0 // 2)
    HOR[T0] = hz
    m = measure(S3, S["ful"], T0, horizon=hz)
    ALL[T0] = m
    c = H
    a = {k: m["mw"][k][c - 2:c + 3][np.argmax(np.abs(m["mw"][k][c - 2:c + 3]))] for k in m["mw"]}
    lagfn = np.nanmedian((m["lf"] - m["ln"])[m["good"]]) * 2
    print(f"  {T0:4d} | {int(m['pol']):+d}  | {100*m['kept']:5.1f} | {m['fn'][0]:5.2f} ({m['fn'][1]:.2f}-{m['fn'][2]:.2f})   {m['mn'][0]:5.2f} ({m['mn'][1]:.2f}-{m['mn'][2]:.2f}) | "
          f"{m['pk_fn'][0]:10.2f} ({m['pk_fn'][1]:.2f}-{m['pk_fn'][2]:.2f}) {m['pk_mn'][0]:5.2f} ({m['pk_mn'][1]:.2f}-{m['pk_mn'][2]:.2f}) | {lagfn:+6.1f} | {a['nea']:+.3f} {a['mid']:+.3f} {a['far']:+.3f}")
    sys.stdout.flush()

# ---------- 2. frequency bands and quadrants for the key events
print("\n== key events: low band (8-22 Hz) against high band (25-50 Hz), and the four quadrants of the cube (far/near)")
print("  time | low far/near  mid/near | high far/near  mid/near | quadrants far/near SW SE NW NE | quadrants mid/near")
bi = (UB // 1000) < NBI / 2
bx = (UB % 1000) < NBX / 2
for T0 in KEY:
    T = min(ev, key=lambda e: abs(e - T0))
    lo_ = measure(S3, S["ful"], T, band=(8, 22), horizon=HOR[T])
    hi_ = measure(S3, S["ful"], T, band=(25, 50), horizon=HOR[T])
    m = ALL[T]
    q, q2 = [], []
    for sel in (bi & bx, ~bi & bx, bi & ~bx, ~bi & ~bx):
        q.append(boot(m["P"][("far", "mid")], m["P"][("nea", "mid")], sel=sel)[0])
        q2.append(boot(m["P"][("far", "mid")], m["P"][("nea", "far")], sel=sel)[0])
    print(f"  {T:4d} | {lo_['fn'][0]:5.2f} {lo_['mn'][0]:9.2f}     | {hi_['fn'][0]:5.2f} {hi_['mn'][0]:10.2f}     | " + " ".join(f"{v:.2f}" for v in q) + "  | " + " ".join(f"{v:.2f}" for v in q2))

# ---------- 3. mean wavelets (flattened on the full-stack horizon, no alignment between stacks)
print("\n== mean wavelets around two events, every 4 ms from -32 to +32 ms (flattened on the full stack, stacks not aligned)")
for T0 in (1316, 1456):
    T = min(ev, key=lambda e: abs(e - T0))
    m = ALL[T]
    for k in ("nea", "mid", "far"):
        print(f"  {T} {k}: " + " ".join(f"{v:+.2f}" for v in m["mw_raw"][k][H - 16:H + 17:2]))

# ---------- 4. multiple interference: does the amplitude of the event 110 ms above control the ratio?
print("\n== does the event one water-layer period (about 110 ms) above control the far/near ratio? Correlation over blocks")
for T0 in (1386, 1456, 1536):
    T = min(ev, key=lambda e: abs(e - T0))
    m = ALL[T]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = m["P"][("far", "mid")] / m["P"][("nea", "mid")]
    for Tp in ev:
        if 96 <= T - Tp <= 126:
            p = ALL[Tp]
            amp_parent = np.sqrt(np.abs(p["P"][("nea", "mid")]))
            okb = np.isfinite(ratio) & np.isfinite(amp_parent) & (ratio > 0) & (ratio < 3)
            cc = np.corrcoef(np.log(ratio[okb]), amp_parent[okb])[0, 1]
            an = np.sqrt(np.abs(m["P"][("nea", "mid")]))
            cc2 = np.corrcoef(an[okb], amp_parent[okb])[0, 1]
            print(f"   event {T} with parent {Tp} (gap {T-Tp} ms): corr(log far/near, parent amplitude) = {cc:+.2f}; corr(own near-mid amplitude, parent amplitude) = {cc2:+.2f}; blocks {okb.sum()}")

# ---------- 5. horizon times and amplitudes at the well
print("\n== at the well (inline index 57, crossline index 143): horizon time and whether the trace was kept")
for T0 in ev:
    tau_s, keep, ccb, lo = HOR[T0]
    print(f"   {T0}: horizon at {T0 + 2*int(tau_s[WELL_IL, WELL_XL])} ms, kept {bool(keep[WELL_IL, WELL_XL])}, correlation with pilot {ccb[WELL_IL, WELL_XL]:.2f}")
