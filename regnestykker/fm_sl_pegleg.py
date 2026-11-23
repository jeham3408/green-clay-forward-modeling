
# ---- fm_sl_pegleg.py: expected delay and stack amplitude of a water-layer peg-leg on near, mid and far ----
# Run:  cat sl_common.py fm_sl_pegleg.py | ssh linux-pc '<geoscience python> -I -'
# Primary at t0p = 1.316 s. First-order peg-leg: one extra round trip in the water layer (sea bed at 0.109 s,
# water 1480 m/s), on the source side or the receiver side. After NMO with the primary velocity at the multiple's
# zero-offset time, the peg-leg is under-corrected. Here: exact ray tracing of both paths through the layered model,
# NMO with the primary model, then stack over each offset range with a Ricker wavelet and NMO stretch.
w, t, L = well()
vel = Vel(w)
TW, VW = 0.109, 1480.0


def tx_primary(x, t0):
    """Two-way time of a primary at zero-offset time t0 for offset x (layered model, exact)."""
    n = vel._n(t0)
    v, h = vel.v[:n], vel.h[:n]
    lo, hi = 0.0, 0.999999 / v.max()
    for _ in range(60):
        p = 0.5 * (lo + hi)
        xx = np.sum(2 * h * p * v / np.sqrt(1 - (p * v) ** 2))
        lo, hi = (p, hi) if xx < x else (lo, p)
    return float(np.sum(2 * h / (v * np.sqrt(1 - (p * v) ** 2)))), p


def tx_pegleg(x, t0):
    """Peg-leg: primary path plus one extra round trip in the water layer at the same ray parameter (exact)."""
    n = vel._n(t0)
    v, h = vel.v[:n], vel.h[:n]
    hw = TW * VW / 2

    def X(p):
        return np.sum(2 * h * p * v / np.sqrt(1 - (p * v) ** 2)) + 2 * hw * p * VW / np.sqrt(1 - (p * VW) ** 2)

    lo, hi = 0.0, 0.999999 / max(v.max(), VW)
    for _ in range(60):
        p = 0.5 * (lo + hi)
        lo, hi = (p, hi) if X(p) < x else (lo, p)
    tt = np.sum(2 * h / (v * np.sqrt(1 - (p * v) ** 2))) + 2 * hw / (VW * np.sqrt(1 - (p * VW) ** 2))
    return float(tt), p


t0p = 1.32                                  # on the 20 ms layer grid of the model (reflector at 1316 ms)
t0p = tx_primary(1e-3, t0p)[0]              # exact zero-offset time of the layered model
t0m = t0p + TW
vr = vel.vrms(t0m)
print(f"primary t0 {t0p:.4f} s, peg-leg t0 {t0m:.4f} s, rms velocity of the model at {t0m:.3f} s: {vr:.0f} m/s")
print(" offset m | primary t(x) | peg-leg t(x) | peg-leg after NMO with the primary velocity: t0' - t0m (ms)")
res = {}
for x in (150, 450, 761, 1000, 1263, 1500, 1763, 2000, 2250, 2500, 2750, 3000, 3150):
    tp, _ = tx_primary(x, t0p)
    tm, _ = tx_pegleg(x, t0m - TW)
    # NMO: map time tm at offset x back to zero offset with the primary hyperbola at t0m (rms velocity)
    t0c = np.sqrt(max(tm ** 2 - (x / vr) ** 2, 0))
    res[x] = (t0c - t0m) * 1000
    print(f" {x:7d}  |   {tp:7.4f}    |   {tm:7.4f}    |   {res[x]:+6.1f}")

# stack a Ricker 35 Hz over each offset range, with NMO stretch, primary and peg-leg; measure delay and amplitude
dt = 0.0005
tt = np.arange(-0.2, 0.4, dt)
geo = {"nea": (150, 761), "mid": (761, 1763), "far": (1763, 3000)}
print("\n stack | peg-leg delay after stacking (ms, relative to t0m) | stacked peg-leg peak / unstretched peak | primary stacked peak")
for k, (a, b) in geo.items():
    xs = np.linspace(a, b, 40)
    tr_m = np.zeros_like(tt)
    tr_p = np.zeros_like(tt)
    for x in xs:
        tp, _ = tx_primary(x, t0p)
        t0pc = np.sqrt(max(tp ** 2 - (x / vel.vrms(t0p)) ** 2, 0)) - t0p
        tm, _ = tx_pegleg(x, t0m - TW)
        t0mc = np.sqrt(max(tm ** 2 - (x / vr) ** 2, 0)) - t0m
        st = t0m / np.sqrt(t0m ** 2 + (x / vr) ** 2)
        tr_m += ricker(tt - t0mc, F0 * st)
        tr_p += ricker(tt - t0pc, F0 * t0p / np.sqrt(t0p ** 2 + (x / vel.vrms(t0p)) ** 2))
    tr_m /= len(xs)
    tr_p /= len(xs)
    i = int(np.argmax(tr_m))
    j = int(np.argmax(tr_p))
    print(f" {k:5s} | copy {TW*1000 + (tt[i] - tt[j])*1000:.1f} ms after the primary | {tr_m[i]:.2f} | {tr_p[j]:.2f}  -> peg-leg/primary on the stack, before the -2R factor: {tr_m[i]/tr_p[j]:.2f}"
          f"; with R = 0.375 the copy ratio would be {-0.75*tr_m[i]/tr_p[j]:+.2f}")
print("observed (fm_sl_mult.py): near -0.75 at 110 ms, mid -0.36 at 112 ms, far -0.27 at 120 ms")
