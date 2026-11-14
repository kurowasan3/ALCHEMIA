# POL-RUB-032 検証：選択性の異なる2種の陽イオンを同時に流したカラムの流出曲線（段モデルによる数値計算）
# 樹脂の選択性 Ca2+ >> K+ >> Na+ を「完全選択」として，各段で樹脂は Ca→K→Na の順に電荷を満たす．
def run(cK, cCa, Q=4.0e-3, N=400, dV=0.5, Vmax=500):
    q=Q/N
    resin=[{'Na':q,'K':0.0,'Ca':0.0} for _ in range(N)]   # 電荷当量（mol of charge）
    out=[]
    V=0.0
    while V<Vmax-1e-9:
        sol={'Na':0.0,'K':cK*dV/1000,'Ca':2*cCa*dV/1000}   # 電荷当量
        for r in resin:
            tot={k:r[k]+sol[k] for k in r}
            cap=q; new={}
            for k in ('Ca','K','Na'):
                take=min(cap,tot[k]); new[k]=take; cap-=take
            for k in r: sol[k]=tot[k]-new[k]
            r.update(new)
        V+=dV
        out.append((V,sol['Na']/(dV/1000),sol['K']/(dV/1000),sol['Ca']/2/(dV/1000)))
    return out
def at(out,v):
    return min(out,key=lambda t:abs(t[0]-v))
o=run(0.020,0.010)
for v in (50,99,101,150,199,201,300):
    print(v,['%.4f'%x for x in at(o,v)[1:]])
assert abs(at(o,50)[1]-0.040)<1e-6 and at(o,50)[2]<1e-9
assert abs(at(o,150)[2]-0.040)<1e-6 and at(o,150)[1]<1e-9 and at(o,150)[3]<1e-9
assert abs(at(o,300)[2]-0.020)<1e-6 and abs(at(o,300)[3]-0.010)<1e-6
dV=0.5
Na=sum(t[1] for t in o)*dV/1000; K=sum(t[2] for t in o if t[0]<=200)*dV/1000
print('流出 Na',Na,'0-200mL の流出 K',K)
assert abs(Na-4.0e-3)<1e-9 and abs(K-0.020*0.2)<1e-9
o2=run(0.020,0.005,Vmax=600)
for v in (100,140,200,390,410):
    print(v,['%.4f'%x for x in at(o2,v)[1:]])
assert abs(at(o2,100)[1]-0.030)<1e-6
assert abs(at(o2,200)[2]-0.030)<1e-6
assert abs(at(o2,410)[2]-0.020)<1e-6 and abs(at(o2,410)[3]-0.005)<1e-6
print('front', 4.0/(0.020+0.010), 4.0/0.010)
print('OK')
