"""Which grid points does each discrete update read? (symbolic only, no run).  Looks for a stencil that is not
centred on the staggered position of the field being updated (that would break the mirror symmetry about z = 0)."""
import json
import re

import numpy as np
from devito import Eq, Grid, Function, TensorTimeFunction, VectorTimeFunction, diag, div, grad

SO = 8
grid = Grid(shape=(41, 41), extent=(80.0, 80.0))
v = VectorTimeFunction(name="v", grid=grid, space_order=SO, time_order=1)
tau = TensorTimeFunction(name="t", grid=grid, space_order=SO, time_order=1)
lam, mu, b = (Function(name=n, grid=grid, space_order=SO) for n in ("lam", "mu", "b"))
s = grid.stepping_dim.spacing
u_v = Eq(v.forward, v + s * b * div(tau))
u_t = Eq(tau.forward, tau + s * (lam * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).transpose(inner=False))))
ev, et = u_v.evaluate, u_t.evaluate
out = {}
def offs(st, name, dim):
    res = set()
    for a in re.findall(name + r"\(([^()]*)\)", st):
        parts = [p.strip() for p in a.split(",")]
        p = parts[-1] if dim == "y" else parts[-2]
        p = p.replace(" ", "")
        m_ = re.fullmatch(dim + r"(?:([+-])(?:(\d+)\*)?h_" + dim + r"(?:/(\d+))?)?", p)
        if m_ is None:
            res.add(p)
            continue
        if m_.group(1) is None:
            res.add(0.0)
        else:
            val = float(m_.group(2) or 1) / float(m_.group(3) or 1)
            res.add(val if m_.group(1) == "+" else -val)
    return sorted(res, key=lambda q: (isinstance(q, str), q))
for nm, lhs, e in (("v_x", v[0], ev.rhs[0]), ("v_y", v[1], ev.rhs[1]), ("t_xx", tau[0, 0], et.rhs[0, 0]), ("t_yy", tau[1, 1], et.rhs[1, 1]), ("t_xy", tau[0, 1], et.rhs[0, 1])):
    st = str(e)
    row = {"lhs_staggered": str(lhs.staggered)}
    for f in ("t_xx", "t_yy", "t_xy", "v_x", "v_y"):
        oy = offs(st, f, "y")
        if oy:
            row[f + "_offsets_in_z_units_of_h"] = oy
    out[nm] = row
    print(nm, row, flush=True)
json.dump(out, open("/workspace/sleipner_verify_etterregning/e11_results.json", "w"), indent=1, default=str)
