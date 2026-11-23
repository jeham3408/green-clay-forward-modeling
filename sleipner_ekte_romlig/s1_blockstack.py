#!/usr/bin/env python3
"""Sleipner 94p07 / 94p10 near, mid, far, full: mean trace per cell of 8 x 9 traces (100 m x 112.5 m).

Read-only on the cubes (memory-mapped, read 8 inlines at a time). Writes blockstack.npz (small) in OUT.
Run on linux-pc:  python -B s1_blockstack.py
"""
import time
import numpy as np

R = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_romlig"
NIL, NXL, NS = 249, 468, 1001
FI, FX = 8, 9
NFI, NFX = (NIL + FI - 1) // FI, NXL // FX        # 32, 52
assert NFX * FX == NXL

res = {}
t0 = time.time()
for proc, d in (("p07", "avo"), ("p10", "avo10")):
    cov = np.load(R + "/avo10/coverage.npy") if proc == "p10" else None
    for key in ("nea", "mid", "far", "ful"):
        a = np.load(f"{R}/{d}/{key}.npy", mmap_mode="r")
        assert a.shape == (NIL, NXL, NS), a.shape
        S = np.zeros((NFI, NFX, NS), np.float64)
        C = np.zeros((NFI, NFX))
        for i0 in range(0, NIL, FI):
            blk = np.array(a[i0:i0 + FI], dtype=np.float32)
            live = np.abs(blk[:, :, 550:950]).sum(-1) > 0
            if cov is not None:
                live &= cov[i0:i0 + FI]
            blk *= live[..., None]
            ni = blk.shape[0]
            S[i0 // FI] = blk.reshape(ni, NFX, FX, NS).sum((0, 2), dtype=np.float64)
            C[i0 // FI] = live.reshape(ni, NFX, FX).sum((0, 2))
        res[f"{proc}_{key}"] = (S / np.maximum(C, 1)[..., None]).astype(np.float32)
        res[f"{proc}_{key}_n"] = C.astype(np.int16)
        del a
        print(proc, key, "live fraction", round(float(C.sum() / (NIL * NXL)), 4), "t", round(time.time() - t0, 1), flush=True)
np.savez(OUT + "/blockstack.npz", **res)
print("done")
