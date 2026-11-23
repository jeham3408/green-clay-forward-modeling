import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
r=json.load(open(OUT+"/results_ekstra.json")); g=json.load(open(OUT+"/results.json"))
fig,ax=plt.subplots(2,2,figsize=(10,6.5),sharex=True)
for j,(v,t) in enumerate([("p07","94p07 ufiltrert"),("p07_dmE","94p07 etter multippelfilter E")]):
    p=np.array(r["profil_"+v]); a=ax[0,j]
    for c,n,col in [(1,"near","C3"),(2,"mid","C2"),(3,"far","C0")]: a.plot(p[:,0]+18,p[:,c],label=n,color=col)
    a.axvspan(1408,1492,color="0.9"); a.set_title(t); a.set_ylabel("RMS rel. referanse (36 ms vinduer)"); a.legend(fontsize=8)
    a=ax[1,j]; a.plot(p[:,0]+18,p[:,4],label="mid/near",color="C1"); a.plot(p[:,0]+18,p[:,5],label="far/near",color="C4")
    a.axhline(1,color="k",lw=.5); a.axvspan(1408,1492,color="0.9"); a.set_xlabel("tid (ms)"); a.legend(fontsize=8)
plt.tight_layout(); plt.savefig(OUT+"/fig_kontroll_profil.png",dpi=80)
gg=g["p07"]["gain_rms_per_96ms"]; print(gg)
