"""(1) energy conservation of the own fluid-solid solution; (2) FD primaries-only stacks against the own Zoeppritz
as-built synthetic (e5), and the Vs effect in both (a grid / physics cross-check that does not need the 1 m run)."""
import numpy as np

import e1_lib as L
import e3_zoep as Z

# (1) fluid-solid: solve the same 3x3 system here and add up the vertical energy fluxes
def fs_full(a1, r1, a2, b2, r2, p):
    p = complex(p)
    vert = lambda v: (lambda q: np.conj(q) if q.imag < 0 else q)(np.sqrt(complex(1.0 / v ** 2 - p ** 2)))
    q1, qa2, qb2 = vert(a1), vert(a2), vert(b2)
    lam2, mu2, K1 = r2 * a2 ** 2 - 2 * r2 * b2 ** 2, r2 * b2 ** 2, r1 * a1 ** 2
    inc = np.array([a1 * q1, K1 * (p * a1 * p + q1 * a1 * q1), 0.0], complex)
    up = np.array([-a1 * q1, K1 * (p * a1 * p + q1 * a1 * q1), 0.0], complex)
    def sol(ux, uz, kz):
        return np.array([uz, lam2 * (p * ux + kz * uz) + 2 * mu2 * kz * uz, mu2 * (kz * ux + p * uz)], complex)
    A = np.column_stack([up, -sol(a2 * p, a2 * qa2, qa2), -sol(b2 * qb2, -b2 * p, qb2)])
    x = np.linalg.solve(A, -inc)
    e = np.array([abs(x[0]) ** 2 * r1 * a1 ** 2 * q1.real, abs(x[1]) ** 2 * r2 * a2 ** 2 * qa2.real, abs(x[2]) ** 2 * r2 * b2 ** 2 * qb2.real]) / (r1 * a1 ** 2 * q1.real)
    return x, e
rng = np.random.default_rng(2)
worst, worst_diff = 0.0, 0.0
for _ in range(2000):
    a2 = rng.uniform(1600, 2500); b2 = a2 / rng.uniform(1.8, 4.0); r2 = rng.uniform(1.7, 2.3)
    th = np.radians(rng.uniform(0, 0.98 * np.degrees(np.arcsin(1480.0 / a2))))
    p = np.sin(th) / 1480.0
    x, e = fs_full(1480.0, 1.03, a2, b2, r2, p)
    worst = max(worst, abs(e.sum() - 1))
    worst_diff = max(worst_diff, abs(x[0].real - Z.fluid_solid_rpp(1480.0, 1.03, a2, b2, r2, p)))
print("fluid-solid energy conservation, 2000 random pre-critical cases: max |sum - 1| = %.2e ; same Rpp as e3_zoep.fluid_solid_rpp: max diff %.1e" % (worst, worst_diff))
R = dict(fluid_solid_energy_max_dev=worst)

# (2)
E = np.load(L.OUT + "/e2_stacks.npz")
Zs = np.load(L.OUT + "/e5_zoep_stacks.npz")
T = E["t"]
w = L.win(1.330, 1.560)
out = {}
for k in ("near", "mid", "far"):
    fd, zo = E["G_trend_w25_a48__" + k], Zs["trend__" + k]
    cc = [float(np.corrcoef(fd[w], np.roll(zo, s)[w])[0, 1]) for s in range(-6, 7)]
    out[k] = dict(corr_best=max(cc), shift_ms=int(np.argmax(cc)) - 6, corr_zero_shift=cc[6])
    for case in ("m07", "m15", "m20"):
        dfd = E["G_%s_w25_a48__%s" % (case, k)] - fd
        dzo = Zs["%s__%s" % (case, k)] - zo
        w2 = L.win(1.380, 1.500)
        sc_fd = L.rms(dfd, w2) / L.rms(fd, w)
        sc_zo = L.rms(dzo, w2) / L.rms(zo, w)
        out[k]["vs_effect_%s" % case] = dict(corr=float(np.corrcoef(dfd[w2], dzo[w2])[0, 1]), strength_fd_over_zoeppritz=float(sc_fd / sc_zo))
    print(k, out[k])
R["fd_vs_zoeppritz"] = out
# far/near RMS in 1330-1560 ms: FD against Zoeppritz (the difference is stretch, spreading, ghost angle dependence, transmission)
R["far_over_near_rms_1330_1560"] = dict(fd=L.rms(E["G_trend_w25_a48__far"], w) / L.rms(E["G_trend_w25_a48__near"], w), zoeppritz=L.rms(Zs["trend__far"], w) / L.rms(Zs["trend__near"], w))
print(R["far_over_near_rms_1330_1560"])
L.jdump(R, "e8_results.json")
print("saved")
