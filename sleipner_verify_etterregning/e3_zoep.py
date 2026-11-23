"""Own Zoeppritz: (1) solid-solid P-SV scattering by solving the 4x4 boundary-condition system numerically for an
incident P wave from above (displacement amplitudes), with an energy-flux check; (2) fluid over solid P reflection.
No closed-form Aki-Richards expression is used for the solid-solid case."""
import numpy as np


def solid_solid(a1, b1, r1, a2, b2, r2, p):
    """Incident P from medium 1 (above). Returns displacement amplitudes (Rpp, Rps, Tpp, Tps) and the energy-flux
    fractions of the four outgoing waves. Complex-safe (post-critical)."""
    p = complex(p)
    def vert(v):
        q = np.sqrt(complex(1.0 / v ** 2 - p ** 2))
        if q.imag < 0:
            q = np.conj(q)
        return q
    qa1, qb1, qa2, qb2 = vert(a1), vert(b1), vert(a2), vert(b2)
    mu1, mu2 = r1 * b1 ** 2, r2 * b2 ** 2
    # plane waves exp(i w (p x +- q z - t)), z down. displacement polarisations (ux, uz):
    #  P down:  (a p,  a q) ; P up: (a p, -a q) ; SV down: (b q, -b p) ; SV up: (b q, b p)
    # tractions on z = const: tzz = i w [lam (dux/dx) + (lam+2mu) duz/dz]/(i w) etc.  (common factor i w dropped)
    def pol(kind, a, b, r, qa, qb):
        lam = r * a ** 2 - 2 * r * b ** 2
        mu = r * b ** 2
        if kind == "Pd":
            ux, uz, kz = a * p, a * qa, qa
        elif kind == "Pu":
            ux, uz, kz = a * p, -a * qa, -qa
        elif kind == "Sd":
            ux, uz, kz = b * qb, -b * p, qb
        else:
            ux, uz, kz = b * qb, b * p, -qb
        tzz = lam * (p * ux + kz * uz) + 2 * mu * kz * uz
        txz = mu * (kz * ux + p * uz)
        return np.array([ux, uz, tzz, txz], complex)
    inc = pol("Pd", a1, b1, r1, qa1, qb1)
    A = np.column_stack([pol("Pu", a1, b1, r1, qa1, qb1), pol("Su", a1, b1, r1, qa1, qb1),
                         -pol("Pd", a2, b2, r2, qa2, qb2), -pol("Sd", a2, b2, r2, qa2, qb2)])
    x = np.linalg.solve(A, -inc)
    rpp, rps, tpp, tps = x
    # vertical energy flux of a unit-displacement-amplitude wave: 0.5 rho v^2 w^2 * Re(q) * v ... proportional to rho * v^2 * Re(q) (per unit amplitude^2, v = wave speed)
    fl = lambda r, v, q: r * v ** 2 * q.real
    e_in = fl(r1, a1, qa1)
    e = np.array([abs(rpp) ** 2 * fl(r1, a1, qa1), abs(rps) ** 2 * fl(r1, b1, qb1), abs(tpp) ** 2 * fl(r2, a2, qa2), abs(tps) ** 2 * fl(r2, b2, qb2)]) / e_in
    return x, e


def rpp(a1, b1, r1, a2, b2, r2, p):
    return solid_solid(a1, b1, r1, a2, b2, r2, p)[0][0]


def fluid_solid_rpp(a1, r1, a2, b2, r2, p):
    """Pressure reflection coefficient, fluid over solid, from the 3x3 boundary system (continuity of uz, tzz = -P, txz = 0)."""
    p = complex(p)
    def vert(v):
        q = np.sqrt(complex(1.0 / v ** 2 - p ** 2))
        return np.conj(q) if q.imag < 0 else q
    q1, qa2, qb2 = vert(a1), vert(a2), vert(b2)
    lam2, mu2 = r2 * a2 ** 2 - 2 * r2 * b2 ** 2, r2 * b2 ** 2
    # fluid: displacement potential; unit incident displacement amplitude along the ray. P down (a p, a q), P up (a p, -a q)
    # fluid tzz = lam1 * div u = r1 a1^2 (p ux + kz uz)
    K1 = r1 * a1 ** 2
    inc = np.array([a1 * q1, K1 * (p * a1 * p + q1 * a1 * q1), 0.0], complex)                # uz, tzz, txz
    up = np.array([-a1 * q1, K1 * (p * a1 * p + q1 * a1 * q1), 0.0], complex)
    def sol(kind):
        if kind == "Pd":
            ux, uz, kz = a2 * p, a2 * qa2, qa2
        else:
            ux, uz, kz = b2 * qb2, -b2 * p, qb2
        return np.array([uz, lam2 * (p * ux + kz * uz) + 2 * mu2 * kz * uz, mu2 * (kz * ux + p * uz)], complex)
    A = np.column_stack([up, -sol("Pd"), -sol("Sd")])
    x = np.linalg.solve(A, -inc)
    return float(x[0].real) if abs(x[0].imag) < 1e-9 else float(np.sign(x[0].real) * abs(x[0])) if abs(x[0].real) > abs(x[0].imag) else float(abs(x[0]))


if __name__ == "__main__":
    # tests: energy conservation, normal incidence, fluid-solid against the impedance form at p = 0
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(2000):
        a1, a2 = rng.uniform(1800, 3500, 2)
        b1, b2 = a1 / rng.uniform(1.6, 3.5), a2 / rng.uniform(1.6, 3.5)
        r1, r2 = rng.uniform(1.8, 2.6, 2)
        th = np.radians(rng.uniform(0, 89))
        x, e = solid_solid(a1, b1, r1, a2, b2, r2, np.sin(th) / a1)
        worst = max(worst, abs(e.sum() - 1.0))
    print("energy conservation, 2000 random interfaces and angles 0-89 deg: max |sum - 1| = %.2e" % worst)
    x, e = solid_solid(2000.0, 800.0, 1.9, 2200.0, 900.0, 2.0, 0.0)
    print("normal incidence Rpp %.6f, impedance formula %.6f (displacement convention: sign -)" % (x[0].real, (2200 * 2.0 - 2000 * 1.9) / (2200 * 2.0 + 2000 * 1.9)))
    print("fluid-solid p=0: %.6f, impedance formula %.6f" % (fluid_solid_rpp(1480.0, 1.03, 1840.6, 678.0, 1.95, 0.0), (1840.6 * 1.95 - 1480 * 1.03) / (1840.6 * 1.95 + 1480 * 1.03)))
