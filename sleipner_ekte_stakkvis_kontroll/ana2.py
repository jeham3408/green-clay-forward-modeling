import numpy as np, json
SL="/mnt/ssd/Data/force-hackathon-2026/sleipner"
OUT="/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
res={}
# 1. profile of ratios in 36 ms windows (3 bins) relative to ref
R=np.r_[25:47,63:91]; tb=808+12*np.arange(97)
for v in ["p07","p10","p07_dmE","p10_dmE"]:
    sq=np.load(OUT+f"/bins_{v}.npz")["sq"]
    E=sq.sum((0,1),dtype=np.float64)
    rr=np.sqrt(E[:,R].mean(1))
    rows=[]
    for s in range(24,92,3):
        W=np.arange(s,s+3)
        q=np.sqrt(E[:,W].mean(1))/rr
        rows.append([int(tb[s]),round(q[0],2),round(q[1],2),round(q[2],2),round(q[1]/q[0],2),round(q[2]/q[0],2)])
    res["profil_"+v]=rows  # t_start, near,mid,far (rel ref), m/n, f/n
# 2. copy regression
S=[np.load(f"{SL}/avo/{k}.npy",mmap_mode="r") for k in ("nea","mid","far")]
tg=[(1366,1466),(1200,1300),(1500,1600),(1600,1700),(1700,1800)]
lagsms=list(range(80,142,2))
acc={}
for s0 in range(3): 
    for (a,b) in tg:
        acc[(s0,a)]=np.zeros((len(lagsms),3))
for i0 in range(0,249,4):
    i1=min(i0+4,249)
    for s0 in range(3):
        c=np.array(S[s0][i0:i1,:,500:900],dtype=np.float32)  # samples 500..899 = 1000..1798ms
        for (a,b) in tg:
            ta,tb_=a//2-500,b//2-500
            x=c[:,:,ta:tb_]
            for li,l in enumerate(lagsms):
                y=c[:,:,ta-l//2:tb_-l//2] if ta-l//2>=0 else None
                if y is None: continue
                acc[(s0,a)][li]+=[ (x*y).sum(), (y*y).sum(), (x*x).sum()]
for s0,nm in enumerate(["near","mid","far"]):
    for (a,b) in tg:
        A=acc[(s0,a)]
        coef=A[:,0]/np.where(A[:,1]>0,A[:,1],1); rho2=A[:,0]**2/np.where(A[:,1]*A[:,2]>0,A[:,1]*A[:,2],1)
        k=int(np.argmin(coef))
        res[f"kopi_{nm}_{a}-{b}"]=dict(lag_ms=lagsms[k],coef=round(float(coef[k]),3),energi_forklart=round(float(rho2[k]*np.sign(-A[k,0])),3) if True else 0)
json.dump(res,open(OUT+"/results_ekstra.json","w"),indent=1)
for k,v in res.items():
    print(k)
    if isinstance(v,list):
        for r in v: print("  ",r)
    else: print("  ",v)
