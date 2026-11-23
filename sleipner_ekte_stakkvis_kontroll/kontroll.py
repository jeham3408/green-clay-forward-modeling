import numpy as np, sys, json, time
SL="/mnt/ssd/Data/force-hackathon-2026/sleipner"
OUT="/mnt/ssd/Data/force-hackathon-2026/review-linux/sleipner_ekte_stakkvis_kontroll"
V={"p07":SL+"/avo/{}.npy","p10":SL+"/avo10/{}.npy","p07_dmE":SL+"/avo/dmE_{}.npy","p10_dmE":SL+"/avo10/dmE_{}.npy"}
S0=404; NB=97; BW=6   # bins of 12 ms from sample 404 (808 ms) to 998
LAGS=list(range(-4,5))
def run(v):
    m=[np.load(V[v].format(k),mmap_mode="r") for k in ("nea","mid","far")]
    NIL,NXL=249,468
    sq=np.zeros((NIL,NXL,3,NB),np.float32); ab=np.zeros((NIL,NXL,3,NB),np.float32)
    cp=np.zeros((NIL,NB,3,len(LAGS)),np.float64)  # pooled over xl per inline: pairs nm, nf, mf with lags
    # cross-products pooled per inline x xl-group of 78 -> do per (il, xlgroup)
    cpg=np.zeros((NIL,6,3,NB,len(LAGS)),np.float64)
    for i0 in range(0,NIL,4):
        i1=min(i0+4,NIL)
        c=[np.array(a[i0:i1,:,S0-6:S0+NB*BW+6],dtype=np.float32) for a in m]  # pad 6 each side
        core=[x[:,:,6:6+NB*BW] for x in c]
        for s in range(3):
            q=core[s].reshape(i1-i0,NXL,NB,BW)
            sq[i0:i1,:,s]=(q**2).sum(-1); ab[i0:i1,:,s]=np.abs(q).sum(-1)
        for pi,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
            for li,l in enumerate(LAGS):
                sh=c[b][:,:,6+l:6+l+NB*BW]
                p=(core[a]*sh).reshape(i1-i0,NXL,NB,BW).sum(-1)  # il,xl,bin
                for g in range(6):
                    cpg[i0:i1,g,pi,:,li]=p[:,g*78:(g+1)*78].sum(1)
        del c,core
    np.savez(OUT+f"/bins_{v}.npz",sq=sq,ab=ab,cpg=cpg)
t=time.time()
for v in sys.argv[1:]:
    run(v); print(v,time.time()-t,flush=True)
