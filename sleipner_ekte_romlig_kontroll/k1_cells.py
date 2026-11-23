#!/usr/bin/env python3
"""Independent check, step 1: cell stacks (8 inlines x 6 crosslines = 100 m x 75 m) of the Sleipner 1994 offset
stacks, samples 500-975 (1000-1950 ms). Memory-mapped, one strip of 8 inlines at a time. Own code, own grid."""
import json
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig_kontroll"
NIL, NXL, NS = 249, 468, 1001
CI, CX = 8, 6
NCI, NCX = NIL // CI, NXL // CX
S0, S1 = 500, 976
out, info = {}, {}
for proc, d in (("p07", "avo"), ("p10", "avo10")):
    cov = None
    if proc == "p10":
        cov = np.load(R + "/avo10/coverage.npy")
        info["coverage"] = [list(cov.shape), str(cov.dtype), float(np.mean(cov != 0))]
        cov = cov != 0
    for key in ("nea", "mid", "far", "ful"):
        a = np.load(f"{R}/{d}/{key}.npy", mmap_mode="r")
        assert a.shape == (NIL, NXL, NS), a.shape
        S = np.zeros((NCI, NCX, S1 - S0))
        C = np.zeros((NCI, NCX))
        top = np.zeros(200)
        for ci in range(NCI):
            blk = np.array(a[ci * CI:(ci + 1) * CI, :NCX * CX, :], dtype=np.float32)
            live = (np.abs(blk[:, :, 550:950]).sum(-1) > 0) & np.isfinite(blk).all(-1)
            if cov is not None:
                live &= cov[ci * CI:(ci + 1) * CI, :NCX * CX]
            blk[~live] = 0
            S[ci] = blk[:, :, S0:S1].reshape(CI, NCX, CX, -1).sum((0, 2), dtype=np.float64)
            C[ci] = live.reshape(CI, NCX, CX).sum((0, 2))
            top += np.abs(blk[:, :, :200]).sum((0, 1), dtype=np.float64)
        out[f"{proc}_{key}"] = (S / np.maximum(C, 1)[..., None]).astype(np.float32)
        out[f"{proc}_{key}_n"] = C.astype(np.int16)
        out[f"{proc}_{key}_top"] = (top / max(C.sum(), 1)).astype(np.float32)
        info[f"{proc}_{key}_live"] = float(C.sum() / (NCI * CI * NCX * CX))
        print(proc, key, info[f"{proc}_{key}_live"], flush=True)
        del a
np.savez(OUT + "/cells.npz", **out)
json.dump(info, open(OUT + "/k1_info.json", "w"), indent=1)
print("done")
