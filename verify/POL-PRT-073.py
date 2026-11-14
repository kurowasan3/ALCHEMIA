# POL-PRT-073 二本鎖DNAの塩基組成・水素結合の数・熱安定性
data={'P':(1000,0.30),'Q':(1000,0.18),'R':(1000,0.25),'S':(1200,0.35)}
res={}
for k,(bp,a) in data.items():
    nb=2*bp; A=round(nb*a); AT=A; GC=bp-AT
    H=2*AT+3*GC; gc=GC/bp
    res[k]=(H,gc); print(k,'AT',AT,'GC',GC,'H',H,'GC frac',gc,'G%',(1-2*a)/2*100,'H per bp',H/bp)
print('Tm order by GC fraction:',sorted(res,key=lambda k:-res[k][1]))
