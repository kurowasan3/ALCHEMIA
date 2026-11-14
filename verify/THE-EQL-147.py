# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/THE-EQL-147/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; depth = 1
        while depth:
            if s[j] == '{': depth += 1
            elif s[j] == '}': depth -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import math, numpy as np
Vinf=50.0; K=2.0e-5
ps=np.array([2.0e4,5.0e4,1.0e5]); Vs=Vinf*K*ps/(1+K*ps)
tab=np.array([14.3,25.0,33.3]); assert np.allclose(np.round(Vs,1),tab)
sl,ic=np.polyfit(ps,ps/tab,1); Vfit=1/sl; Kfit=1/(ic*Vfit)
print(Vfit,Kfit); assert sig(Vfit,2)==50 and sig(Kfit,2)==2.0e-5
nmol=50e-3/22.4; area=nmol*6.02e23*2.0e-19; assert sig(area/2,2)==130
th1=K*1e5/(1+K*1e5); th2=K*2e5/(1+K*2e5); assert sig(th2/th1,2)==1.2
assert sig(9.0/K,2)==4.5e5
N=1e-4/50e-18; lat=N*2*math.pi*3.0e-9*2.0e-6; faces=2e-4*(1-math.pi*(3.0e-9)**2/50e-18)
r=(lat+faces)/2e-4; print(N,lat,faces,r); assert sig(r,2)==380

print('ALL OK')
