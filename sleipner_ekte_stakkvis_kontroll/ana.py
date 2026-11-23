import numpy as np, json
OUT="/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
LAGS=list(range(-4,5))
Z=np.arange(50,57); R=np.r_[25:47,63:91]
tb=808+12*np.arange(97)
res={}
def r3(x): return round(float(x),3)
for v in ["p07","p10","p07_dmE","p10_dmE"]:
    d=np.load(OUT+f"/bins_{v}.npz"); sq=d["sq"]; ab=d["ab"]; cpg=d["cpg"]
    o={}
    # whole cube pooled
    E=sq.sum((0,1),dtype=np.float64); A=ab.sum((0,1),dtype=np.float64)   # [3,NB]
    rms=lambda W,E=E: np.sqrt(E[:,W].sum(1)/len(W))
    mab=lambda W,A=A: A[:,W].sum(1)/len(W)
    o["rms_zr"]=[r3(a/b) for a,b in zip(rms(Z),rms(R))]
    o["mabs_zr"]=[r3(a/b) for a,b in zip(mab(Z),mab(R))]
    # per trace median
    tz=np.sqrt(sq[:,:,:,Z].mean(-1)); tr=np.sqrt(sq[:,:,:,R].mean(-1))
    ok=(tr>0).all(-1)&(tz>0).all(-1)
    pr=(tz/np.where(tr>0,tr,1))[ok]
    o["trace_median_zr"]=[r3(np.median(pr[:,s])) for s in range(3)]
    o["trace_p16_p84"]=[[r3(np.percentile(pr[:,s],16)),r3(np.percentile(pr[:,s],84))] for s in range(3)]
    mn=pr[:,1]/pr[:,0]; fn=pr[:,2]/pr[:,0]
    o["trace_median_mid_near"]=r3(np.median(mn)); o["trace_median_far_near"]=r3(np.median(fn))
    # sub-areas 18 : 83 il x 78 xl
    sub=[]
    for gi in range(3):
        for gx in range(6):
            e=sq[gi*83:(gi+1)*83,gx*78:(gx+1)*78].sum((0,1),dtype=np.float64)
            sub.append(np.sqrt(e[:,Z].mean(1))/np.sqrt(e[:,R].mean(1)))
    sub=np.array(sub)  # 18,3
    for nm,val in [("near",sub[:,0]),("mid",sub[:,1]),("far",sub[:,2]),("mid/near",sub[:,1]/sub[:,0]),("far/near",sub[:,2]/sub[:,0]),("far/mid",sub[:,2]/sub[:,1])]:
        o["sub_"+nm]=[r3(np.median(val)),r3(val.std(ddof=1)/np.sqrt(18)),r3(val.min()),r3(val.max())]
    # sliding window placebo: all 7-bin windows with start 0.. in 1100-1900 excluding overlap with 1380-1560 region
    def ratio_for(W,Rb):
        e=E[:,W].mean(1); er=E[:,Rb].mean(1)
        return np.sqrt(e/er)
    pl=[]
    starts=[s for s in range(24,90-6) if (s+6<47 or s>=63)]
    for s in starts:
        W=np.arange(s,s+7); Rb=np.array([b for b in R if b<s-1 or b>s+7])
        q=ratio_for(W,Rb); pl.append([s,q[0],q[1],q[2],q[1]/q[0],q[2]/q[0],q[2]/q[1]])
    pl=np.array(pl)
    names=["near","mid","far","mid/near","far/near","far/mid"]
    o["placebo_n"]=len(pl)
    zq=ratio_for(Z,R); zv=[zq[0],zq[1],zq[2],zq[1]/zq[0],zq[2]/zq[0],zq[2]/zq[1]]
    for i,n in enumerate(names):
        col=pl[:,1+i]
        o["placebo_"+n]=dict(zone=r3(zv[i]),mean=r3(col.mean()),sd=r3(col.std(ddof=1)),min=r3(col.min()),max=r3(col.max()),z=r3((zv[i]-col.mean())/col.std(ddof=1)))
    # cross-product estimators, whole cube per lag: cpg [il,6,3pairs,NB,lag]
    C=cpg.sum((0,1))  # pairs,NB,lag
    l0=LAGS.index(0)
    def est(W):
        out={}
        nn,mm,ff=E[0,W].sum(),E[1,W].sum(),E[2,W].sum()
        nm=C[0,W,l0].sum(); nf=C[1,W,l0].sum(); mf=C[2,W,l0].sum()
        out["corr_nm"]=nm/np.sqrt(nn*mm); out["corr_nf"]=nf/np.sqrt(nn*ff); out["corr_mf"]=mf/np.sqrt(mm*ff)
        out["rms_mn"]=np.sqrt(mm/nn); out["rms_fn"]=np.sqrt(ff/nn)
        out["ls_mn"]=nm/nn      # regression of mid on near
        out["ls_fn"]=nf/nn
        out["ls_fm"]=mf/mm
        # instrument est using other as aux: mid/near = (mid*far)/(near*far)
        out["iv_mn"]=mf/nf
        out["iv_fn_viamid"]=mf/nm
        # best lag
        for pi,n in enumerate(["nm","nf","mf"]):
            cl=C[pi,W,:].sum(0); out["bestlag_"+n]=int(np.array(LAGS)[np.argmax(cl)])*2
        return out
    ez,er=est(Z),est(R)
    o["est_zone"]={k:r3(x) for k,x in ez.items()}; o["est_ref"]={k:r3(x) for k,x in er.items()}
    o["est_zr"]={k:r3(ez[k]/er[k]) for k in ["rms_mn","rms_fn","ls_mn","ls_fn","ls_fm","iv_mn","iv_fn_viamid"]}
    # gain: RMS per 96 ms
    g=[]
    for w0 in range(0,96,8):
        W=np.arange(w0,min(w0+8,97)); 
        if len(W)<8: break
        g.append([int(tb[w0]),*[r3(x) for x in np.sqrt(E[:,W].mean(1)/(249*468*6))]])
    o["gain_rms_per_96ms"]=g
    del sq,ab,tz,tr,pr
    res[v]=o
json.dump(res,open(OUT+"/results.json","w"),indent=1)
for v,o in res.items():
    print("=====",v)
    for k,x in o.items():
        if k!="gain_rms_per_96ms": print(k,x)
