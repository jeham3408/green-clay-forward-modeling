
# ---- fm_sl_trend167.py: the base of the softest layer in 15/9-13 with the upper-Hordaland shale line of 16/7-6
# (Vs = -383 + 0.538 Vp, measured DTS, 14 km away), compared with the other backgrounds ----
# Run:  cat sl_common.py fm_sl_trend167.py | ssh linux-pc '<geoscience python> -I -'
w, t, L = well()
vel = Vel(w)
g = stack_geometry(vel, 1.45)
a1, r1, a2, r2 = 2041.0, 1.885, 2212.0, 2.026
TR = {
    "16/7-6 upper Hordaland": lambda v: -383.0 + 0.538 * v,
    "Nini-4 Lark": vs_trend,
    "15/9-23 deeper Hordaland": lambda v: 415.0 + 0.225 * v,
}
for lab, fn in TR.items():
    base = None
    print(f"\n== background {lab}: soft layer Vs {fn(a1):.0f} (Vp/Vs {a1/fn(a1):.2f}) over Vs {fn(a2):.0f} (Vp/Vs {a2/fn(a2):.2f})")
    for red in (0, 5, 7, 10, 15, 20, 25):
        b1v = fn(a1) * (1 - red / 100)
        rr = {k: float(np.mean(zoep(a1, b1v, r1, a2, fn(a2), r2, g[k][0]))) for k in g}
        fn_ = rr["far"] / rr["nea"]
        base = fn_ if red == 0 else base
        print(f"   soft-layer Vs {red:2d} % below: far/near {fn_:.3f} ({100*(fn_/base-1):+.0f} %), mid/near {rr['mid']/rr['nea']:.3f}")
    # Vs anomaly needed for -35 % and -40 %
    for target in (0.65, 0.60):
        lo, hi = 0.0, 0.6
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            b1v = fn(a1) * (1 - mid)
            rr = {k: float(np.mean(zoep(a1, b1v, r1, a2, fn(a2), r2, g[k][0]))) for k in g}
            if rr["far"] / rr["nea"] / base > target:
                lo = mid
            else:
                hi = mid
        print(f"   anomaly for far/near x{target:.2f}: {100*mid:.1f} % below the background")
