# Windowed group file: each gene's var/anno lists cut into consecutive windows
# of W markers (region GENE_wNN). Same markers, same labels; only the grouping
# changes, which a null phenotype cannot see. Prints predicted cell counts by
# unweighted MAC stratum for the No_Weights arm.
import collections, sys
D='/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k'
W=int(sys.argv[1]); out=sys.argv[2]; pred=sys.argv[3]   # e.g. 25 $CALIB_DIR/groupFile_w25.txt $CALIB_DIR/predicted_w25.tsv
mac={}
with open(f'{D}/allcounts.acount') as f:
    next(f)
    for l in f:
        c=l.split(); a=int(c[5]); n=int(c[6]); mac[c[1]]=min(a,n-a)
genes=collections.OrderedDict()
with open(f'{D}/groupFile.txt') as f:
    for l in f:
        c=l.split(); genes.setdefault(c[0],{})[c[1]]=c[2:]
N=50000
groups={'high':{'high'},'high:medium':{'high','medium'}}
mafs=[1e-4,1e-3]
cnt={(g,M):collections.Counter() for g in groups for M in mafs}
nreg=0
with open(out,'w') as o, open(pred,'w') as p:
    p.write("region\tgroup\tmaxMAF\tMAC\tn\n")
    for g,d in genes.items():
        v,a=d['var'],d['anno']; assert len(v)==len(a)
        for w in range(0,len(v),W):
            r=f"{g}_w{w//W+1:02d}"; vv=v[w:w+W]; aa=a[w:w+W]; nreg+=1
            o.write(f"{r} var "+" ".join(vv)+"\n"); o.write(f"{r} anno "+" ".join(aa)+"\n")
            for gn,annos in groups.items():
                for M in mafs:
                    s=0; nm=0
                    for vid,an in zip(vv,aa):
                        m=mac[vid]
                        if m>=1 and an in annos and m/(2*N)<=M: s+=m; nm+=1
                    p.write(f"{r}\t{gn}\t{M:g}\t{s}\t{nm}\n")
                    cnt[(gn,M)]['0' if s==0 else str(s) if s<5 else '>=5']+=1
print("regions:",nreg, "window:",W)
for k,c in cnt.items():
    print(f"{k[0]:12s} maxMAF={k[1]:g}  "+"  ".join(f"{s}:{c[s]}" for s in ['0','1','2','3','4','>=5']))
# The reduced file (regions whose high:medium x 1e-3 cell has raw MAC >= 5), the
# `wgridr` arm: the SKAT-O path dies at a raised gate on a region whose every
# cell is refused (LEDGER #128), so its gate-5 control runs on this subset.
keep=set()
with open(pred) as f:
    for l in f:
        c=l.rstrip("\n").split("\t")
        if c[0]!="region" and c[1]=="high:medium" and c[2]=="0.001" and int(c[3])>=5: keep.add(c[0])
with open(out.replace(".txt","r.txt"),"w") as o:
    for l in open(out):
        if l.split()[0] in keep: o.write(l)
print("reduced file:",out.replace(".txt","r.txt"),"regions",len(keep))
