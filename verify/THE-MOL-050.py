# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/out/items/THE-MOL-050/item.json'
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

from sympy import symbols, solve, Rational as R
x,y=symbols('x y')
Ag,S,O=108,32,16
mAg2S=2*Ag+S; mAg2SO4=2*Ag+S+4*O
n0=R(1240,100)/mAg2S
sol=solve([x+y-n0, R(1240,100)-(mAg2S-2*Ag)*x+(mAg2SO4-mAg2S)*y-R(1224,100)],[x,y])
print(sol); assert abs(float(sol[x])-0.0350)<1e-9 and abs(float(sol[y])-0.0150)<1e-9
final=2*n0*Ag; assert sig(float(final),3)==10.8
so2=float(n0*R(224,10)); o2=float((sol[x]+2*sol[y]-sol[y])*R(224,10))
assert sig(so2,3)==1.12 and sig(o2,3)==1.12
from collections import Counter
L=Counter({'Ag2S':1,'O2':2,'Ag2SO4':1}); Rr=Counter({'Ag2SO4':1,'Ag':2,'SO2':1,'O2':1})
netL=L-Rr; netR=Rr-L
print('Q+R 正味:',dict(netL),'->',dict(netR))
assert dict(netL)=={'Ag2S':1,'O2':1} and dict(netR)=={'Ag':2,'SO2':1}

print('ALL OK')
