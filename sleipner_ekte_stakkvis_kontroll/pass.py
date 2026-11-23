#!/usr/bin/env python3
"""pass.py (kontroll): one pass over the Sleipner stacks, written independently of sleipner_ekte_stakkvis.
Per version it accumulates, per sub-area (24) and per time sample (800-1998 ms):
  sum a^2, sum |a|, neighbour products (crossline lag 4 and inline lag 4), cross products between stacks at lags
  -4..4 samples, lagged auto products at 96..132 ms (the "copy"), and block-mean traces (10 x 12 traces).
Reads memory-mapped cubes 12-16 inlines at a time; the map is reopened for each chunk.
Run: run_capped.sh 600 python -B pass.py ver [ver ...]
"""
import os, sys, time, resource
import numpy as np

SL = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
OUT = "/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
VERS = {"p07": SL + "/avo/{}.npy", "p10": SL + "/avo10/{}.npy", "p07_dmE": SL + "/avo/dmE_{}.npy",
        "p10_dmE": SL + "/avo10/dmE_{}.npy", "p07_dm": SL + "/avo/dm_{}.npy", "p07_dmA": SL + "/avo/dmA_{}.npy",
        "p07_dmB": SL + "/avo/dmB_{}.npy", "p07_dmD": SL + "/avo/dmD_{}.npy"}
KS = ("nea", "mid", "far")
NIL, NXL, NS = 249, 468, 1001
S0, S1 = 400, 1000
NT = S1 - S0
IL_EDGES = [0, 62, 124, 186, 249]            # 4 inline bands
XL_EDGES = [0, 78, 156, 234, 312, 390, 468]  # 6 crossline bands -> 24 sub-areas of about 0.78 km x 0.98 km
NB = 24
XLAGS = list(range(-4, 5))
ALAGS = list(range(48, 67))
PAIRS = (("nea", "mid"), ("nea", "far"), ("mid", "far"))
BI, BX = 10, 12                              # block-mean traces 125 m x 150 m
NBI, NBX = NIL // BI, NXL // BX              # 24 x 39 (last 9 inlines dropped)


def lagprod(a, b, L):
    """sum over inlines of a(t)*b(t+L), returned on the full time axis (zero where t+L is outside)."""
    out = np.zeros((a.shape[1], NT))
    if L >= 0:
        out[:, :NT - L] = np.einsum("ixt,ixt->xt", a[:, :, :NT - L], b[:, :, L:])
    else:
        out[:, -L:] = np.einsum("ixt,ixt->xt", a[:, :, -L:], b[:, :, :NT + L])
    return out


def xlsum(p):
    return np.add.reduceat(p, XL_EDGES[:-1], axis=0)       # (6, NT)


def run(ver):
    t0 = time.time()
    pat = VERS[ver]
    A = dict(ss={k: np.zeros((NB, NT)) for k in KS}, ab={k: np.zeros((NB, NT)) for k in KS},
             nx={k: np.zeros((NB, NT)) for k in KS}, ni={k: np.zeros((NB, NT)) for k in KS},
             cx={p: np.zeros((len(XLAGS), NB, NT)) for p in PAIRS}, ac={k: np.zeros((len(ALAGS), NB, NT)) for k in KS})
    cnt = np.zeros(NB)
    cnt_nx = np.zeros(NB)
    cnt_ni = np.zeros(NB)
    tr_ss = {k: np.zeros((NIL, NXL, 12), np.float32) for k in KS}     # per trace energy in 100 ms windows
    B = {k: np.zeros((NBI, NBX, NT), np.float32) for k in KS}
    live_map = np.zeros((NIL, NXL), bool)
    for ib in range(4):
        a0, a1 = IL_EDGES[ib], IL_EDGES[ib + 1]
        for c0 in range(a0, a1, 16):
            c1 = min(a1, c0 + 16)
            X = {}
            for k in KS:
                mm = np.load(pat.format(k), mmap_mode="r")
                assert mm.shape == (NIL, NXL, NS)
                X[k] = np.array(mm[c0:c1, :, S0:S1], dtype=np.float32)
                del mm
            live = np.ones((c1 - c0, NXL), bool)
            for k in KS:
                e = (X[k][:, :, 150:550].astype(np.float64) ** 2).sum(-1)
                live &= np.isfinite(e) & (e > 0)
            for k in KS:
                X[k][~live] = 0.0
                X[k] = np.nan_to_num(X[k])
            live_map[c0:c1] = live
            sl = slice(ib * 6, ib * 6 + 6)
            cnt[sl] += np.add.reduceat(live.sum(0), XL_EDGES[:-1])
            for k in KS:
                x = X[k]
                A["ss"][k][sl] += xlsum(np.einsum("ixt,ixt->xt", x, x))
                A["ab"][k][sl] += xlsum(np.abs(x).sum(0, dtype=np.float64))
                tr_ss[k][c0:c1] = (x.astype(np.float64) ** 2).reshape(c1 - c0, NXL, 12, 50).sum(-1)
                # neighbour products: crossline lag 4 (50 m), inline lag 4 (50 m, inside the chunk)
                p = np.zeros((NXL, NT))
                p[:-4] = np.einsum("ixt,ixt->xt", x[:, :-4], x[:, 4:])
                A["nx"][k][sl] += xlsum(p)
                if c1 - c0 > 4:
                    A["ni"][k][sl] += xlsum(np.einsum("ixt,ixt->xt", x[:-4], x[4:]))
                for j, L in enumerate(ALAGS):
                    A["ac"][k][j, sl] += xlsum(lagprod(x, x, L))
            pl = np.zeros((NXL,))
            pl[:-4] = (live[:, :-4] & live[:, 4:]).sum(0)
            cnt_nx[sl] += np.add.reduceat(pl, XL_EDGES[:-1])
            if c1 - c0 > 4:
                cnt_ni[sl] += np.add.reduceat((live[:-4] & live[4:]).sum(0), XL_EDGES[:-1])
            for (k1, k2) in PAIRS:
                for j, L in enumerate(XLAGS):
                    A["cx"][(k1, k2)][j, sl] += xlsum(lagprod(X[k1], X[k2], L))
            # block-mean traces
            for r in range(NBI):
                lo, hi = max(c0, r * BI), min(c1, (r + 1) * BI)
                if hi <= lo:
                    continue
                for k in KS:
                    B[k][r] += X[k][lo - c0:hi - c0, :NBX * BX].reshape(hi - lo, NBX, BX, NT).sum((0, 2)) / (BI * BX)
            del X
    np.savez_compressed(os.path.join(OUT, f"acc_{ver}.npz"), cnt=cnt, cnt_nx=cnt_nx, cnt_ni=cnt_ni,
                        live=live_map, **{f"ss_{k}": A["ss"][k] for k in KS}, **{f"ab_{k}": A["ab"][k] for k in KS},
                        **{f"nx_{k}": A["nx"][k] for k in KS}, **{f"ni_{k}": A["ni"][k] for k in KS},
                        **{f"cx_{a}_{b}": A["cx"][(a, b)] for a, b in PAIRS}, **{f"ac_{k}": A["ac"][k] for k in KS},
                        **{f"tr_{k}": tr_ss[k] for k in KS}, **{f"B_{k}": B[k] for k in KS})
    mx = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"{ver}: live {live_map.mean():.4f}  {time.time() - t0:.0f} s  maxrss {mx:.0f} MB", flush=True)


if __name__ == "__main__":
    for v in sys.argv[1:]:
        run(v)
