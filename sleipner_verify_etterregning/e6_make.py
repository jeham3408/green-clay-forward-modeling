"""Builds e6_fs_surface.py: a copy of the experiment's s2_fd.py (read-only source) with extra receiver rows at
z = 0, 2, 4 and -8 m in the mirror (FS) configuration, to measure the pressure at the sea surface directly."""
src = open("/workspace/sleipner_exp/s2_fd.py").read()
a = 'zr = [ZREC] if cfg == "FS" else [ZREC, -ZREC]'
b = 'zr = [ZREC, 0.0, 2.0, 4.0, -ZREC] if cfg == "FS" else [ZREC, -ZREC]'
assert src.count(a) == 1
src = src.replace(a, b)
a = '    d = np.array(recs[0].data, np.float64)\n'
b = a + ('    if cfg == "FS":\n'
         '        raw = {str(zz): np.array(r.data, np.float64) for zz, r in zip(zr, recs)}\n'
         '        ref = np.abs(raw["8.0"]).max(axis=0)\n'
         '        chk = {"max_abs_p_z%s_over_max_abs_p_z8_per_trace_max" % k: float((np.abs(v).max(axis=0) / ref).max()) for k, v in raw.items()}\n'
         '        chk["max_abs_p_z8"] = float(ref.max())\n'
         '        chk["max_abs_of_p(z=-8)+p(z=8)_over_max_abs_p_z8"] = float(np.abs(raw["-8.0"] + raw["8.0"]).max() / ref.max())\n'
         '        chk["max_abs_p_z0_abs"] = float(np.abs(raw["0.0"]).max())\n'
         '        json.dump(chk, open(os.path.join(ROOT, "e6_surface_pressure_%s.json" % case), "w"), indent=1)\n'
         '        print("SURFACE", json.dumps(chk), flush=True)\n')
assert src.count(a) == 1
src = src.replace(a, b)
open("/workspace/sleipner_verify_etterregning/e6_fs_surface.py", "w").write(src)
print("written")
