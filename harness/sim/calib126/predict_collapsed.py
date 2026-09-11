# Predicted cell MAC under the engine's rule (LEDGER #118 units): kept markers =
# annotation in group, 1 <= MAC, MAF <= maxMAF (and MAC <= t for MAC-threshold
# sets, high-only for lof sets); MAC = sum over non-ultra-rare kept markers of
# their MAC + sum over samples of max over ultra-rare (MAC <= 10) kept markers.
import sys, numpy as np, collections
D='/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k'
import os; S=os.environ['CALIB_DIR']
N=50000; B=(N+3)//4; UR=10
idx={}
with open(f'{D}/chr20_sim.bim') as f:
    for i,l in enumerate(f): idx[l.split()[1]]=i
mac={}
with open(f'{D}/allcounts.acount') as f:
    next(f)
    for l in f:
        c=l.split(); a=int(c[5]); n=int(c[6]); mac[c[1]]=min(a,n-a)
lut=np.zeros((256,4),dtype=np.int8)
for b in range(256):
    for k in range(4):
        v=(b>>(2*k))&3; lut[b,k]={0:2,1:-1,2:1,3:0}[v]
bed=open(f'{D}/chr20_sim.bed','rb')
gcache={}
def geno(x):
    if x in gcache: return gcache[x]
    bed.seek(3+idx[x]*B); raw=np.frombuffer(bed.read(B),dtype=np.uint8)
    g=lut[raw].reshape(-1)[:N].astype(np.int8); g=np.where(g<0,0,g)
    if g.sum()!=mac[x]: g=2-g
    assert g.sum()==mac[x]
    gcache[x]=g; return g
def cellmac(v,a,annos,maxmaf,maxmac=None):
    kept=[x for x,an in zip(v,a) if an in annos and mac[x]>=1 and mac[x]/(2*N)<=maxmaf and (maxmac is None or mac[x]<=maxmac)]
    if not kept: return 0,0
    s=sum(mac[x] for x in kept if mac[x]>UR)
    ur=[x for x in kept if mac[x]<=UR]
    if ur:
        m=np.zeros(N,dtype=np.int8)
        for x in ur: m=np.maximum(m,geno(x))
        s+=int(m.sum())
    return s,len(kept)
def load(fn):
    g=collections.OrderedDict()
    for l in open(fn):
        c=l.split(); g.setdefault(c[0],{})[c[1]]=c[2:]
    return g
groups={'high':{'high'},'high:medium':{'high','medium'},'high:medium:low':{'high','medium','low'}}
with open(f'{S}/predicted_cell_mac_collapsed.tsv','w') as o:
    o.write("gene\tgroup\tmaxMAF\tMAC\tn\n")
    for g,d in load(f'{D}/groupFile.txt').items():
        v,a=d['var'],d['anno']
        for gn,annos in groups.items():
            for M in (1e-4,1e-3,1e-2):
                s,nm=cellmac(v,a,annos,M); o.write(f"{g}\t{gn}\t{M:g}\t{s}\t{nm}\n")
        sets={'lof_unw':({'high'},None),'unw_unw':({'high','medium'},None)}
        for t in (1,3,5,10,25,50,100,200,500):
            sets[f'lof_MAC{t}']=({'high'},t); sets[f'unw_MAC{t}']=({'high','medium'},t)
        for nm,(annos,mm) in sets.items():
            s,n2=cellmac(v,a,annos,1e-3,mm); o.write(f"{g}\tflex:{nm}\t0.001\t{s}\t{n2}\n")
    gcache.clear()
with open(f'{S}/predicted_w25_collapsed.tsv','w') as o:
    o.write("gene\tgroup\tmaxMAF\tMAC\tn\n")
    for g,d in load(f'{S}/groupFile_w25.txt').items():
        v,a=d['var'],d['anno']
        for gn in ('high','high:medium'):
            for M in (1e-4,1e-3):
                s,nm=cellmac(v,a,groups[gn],M); o.write(f"{g}\t{gn}\t{M:g}\t{s}\t{nm}\n")
print("done")
