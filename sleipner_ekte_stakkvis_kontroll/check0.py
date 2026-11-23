#!/usr/bin/env python3
"""check0.py: time origin, polarity at the sea bed, well position, raw SEG-Y delay; tiny reads only."""
import json, struct, numpy as np
SL = "/mnt/ssd/Data/force-hackathon-2026/sleipner"
out = {}
for ver, pat in (("p07", SL + "/avo/{}.npy"), ("p10", SL + "/avo10/{}.npy")):
    for k in ("nea", "mid", "far"):
        mm = np.load(pat.format(k), mmap_mode="r")
        x = np.array(mm[::25, ::25, 0:120], dtype=np.float64)
        del mm
        m = x.reshape(-1, 120).mean(0)
        j = int(np.argmax(np.abs(m)))
        jp = int(np.argmax(m))
        out[f"{ver}_{k}"] = dict(strongest_sample_ms=2 * j, strongest_value=float(m[j]), max_positive_ms=2 * jp,
                                 max_positive=float(m[jp]), first_sample_above_10pct_ms=int(2 * np.argmax(np.abs(m) > 0.1 * np.abs(m).max())))
# raw SEG-Y header of the 2007 near stack: sample interval, delay recording time, first trace il/xl
p = SL + "/raw/seismic/2007 processing/data/94p07nea.sgy"
with open(p, "rb") as f:
    f.seek(3200); b = f.read(400); th = f.read(240)
out["segy07"] = dict(dt_us=struct.unpack(">H", b[16:18])[0], ns=struct.unpack(">H", b[20:22])[0],
                     delay_ms=struct.unpack(">h", th[108:110])[0], lag_a=struct.unpack(">h", th[104:106])[0],
                     il=struct.unpack(">i", th[188:192])[0], xl=struct.unpack(">i", th[192:196])[0])
w = json.load(open(SL + "/cubes/sleipner94/wells.json"))["wells"][0]
meta = json.load(open(SL + "/cubes/sleipner94/meta.json"))
ox, oy = meta["origin_xy"]; ia = np.array(meta["il_axis_m"]); xa = np.array(meta["xl_axis_m"])
M = np.array([ia, xa]).T
ij = np.linalg.solve(M, np.array([w["x"] - ox, w["y"] - oy]))
out["well"] = dict(il_json=w["il"], xl_json=w["xl"], idx_from_il=w["il"] - meta["il0"], idx_from_xl=w["xl"] - meta["xl0"],
                   idx_from_xy=[float(ij[0]), float(ij[1])], seabed_twt=w["seabed_twt"])
print(json.dumps(out, indent=1))
json.dump(out, open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll/check0.json", "w"), indent=1)
