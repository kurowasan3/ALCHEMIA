# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/outh/items/THE-EQL-141/item.json'
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

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans <= 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

k1,km1,k2,km2=1.0,0.5,0.05,0.005
assert k1/k2==20
nA=1/13; assert sig(nA,2)==0.077 and sig(2*nA,2)==0.15 and sig(10*nA,2)==0.77
A,B,C=1.0,0.0,0.0; dt=1e-4; t=0; mx=(0,0); th=None
while t<400:
    dA=-(k1+k2)*A+km1*B+km2*C; dB=k1*A-km1*B; dC=k2*A-km2*C
    A+=dA*dt; B+=dB*dt; C+=dC*dt; t+=dt
    if B>mx[0]: mx=(B,t)
    if th is None and C>=0.5*10/13: th=t
assert abs(C-10/13)<1e-3
assert sig(mx[0],1)==0.6 and 2.5<mx[1]<3.5, mx
ta=0.69/(k2/3+km2); assert sig(ta,2)==32 and abs(ta-th)/th<0.03, (ta,th)
print('ALL OK')
