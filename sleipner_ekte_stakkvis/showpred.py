import json
P = json.load(open("/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_exp/prediksjoner.json"))
for k in P:
    if k != "prediksjoner":
        print(k, ":", json.dumps(P[k], ensure_ascii=False)[:1500])
for s, d in P["prediksjoner"].items():
    print("==", s, "|", d.get("offsetomraader"))
    for q, e in d["forklaringer"].items():
        print(f"  {q:38s} near {e['near_alene']:.3f} mid {e['mid_alene']:.3f} far {e['far_alene']:.3f} | mid/near rms {e['mid_over_near']['rms']:.3f} refl {e['mid_over_near']['per_reflektor']:.3f} topp {e['mid_over_near']['per_reflektor_topp']:.3f} | far/near rms {e['far_over_near']['rms']:.3f} refl {e['far_over_near']['per_reflektor']:.3f} topp {e['far_over_near']['per_reflektor_topp']:.3f}")
