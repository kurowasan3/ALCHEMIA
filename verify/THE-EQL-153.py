# 検証スクリプト（グループ g05）  python3 -I このファイル
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemk/outh/items/THE-EQL-153/item.json'
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

from sympy import symbols, solve, diff, simplify, Rational as R
tA,tB,KA,KB,pA,pB,x,a=symbols('tA tB KA KB pA pB x a', positive=True)
s=solve([tA-KA*pA*(1-tA-tB), tB-KB*pB*(1-tA-tB)],[tA,tB])
assert simplify(s[tA]-KA*pA/(1+KA*pA+KB*pB))==0 and simplify(s[tB]-KB*pB/(1+KA*pA+KB*pB))==0
crit=solve(diff(x/(a+x)**2,x),x); assert crit==[a]
k_A,k_B=symbols('k_A k_B')
sol=solve([3*k_B-1-2*k_A, 5*k_B-1-6*k_A],[k_A,k_B]); assert sol=={k_A:R(1,4),k_B:R(1,2)}
import numpy as np
for pa,expect in [(2.0,3.0),(6.0,5.0)]:
    pb=np.linspace(0.01,40,400001); v=0.25*pa*0.5*pb/(1+0.25*pa+0.5*pb)**2
    assert abs(pb[v.argmax()]-expect)<1e-3
vm=lambda pa:0.25*pa/(4*(1+0.25*pa)); assert sig(vm(6.0)/vm(2.0),2)==1.8
f=lambda pb:0.5*pb/(1.5+0.5*pb)**2; assert sig(f(12)/f(3),2)==0.64
thA=lambda pb:0.5/(1.5+0.5*pb); thB=lambda pb:0.5*pb/(1.5+0.5*pb)
print(thA(3),thA(12),thB(3),thB(12))
assert abs(thA(3)-1/6)<1e-12 and abs(thA(12)-1/15)<1e-12 and thB(3)==0.5 and thB(12)==0.8
# 機構 I：dv/dpB>0
vI=KA*pA/(1+KA*pA)*pB; assert simplify(diff(vI,pB)-KA*pA/(1+KA*pA))==0
print('数値OK')

print('ALL OK')
