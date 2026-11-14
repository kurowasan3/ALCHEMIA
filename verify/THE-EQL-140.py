# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/outh/items/THE-EQL-140/item.json'
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

def f(pH,p1,p2):
    h=10**-pH; K1=10**-p1; K2=10**-p2
    return K1*h/((h+K1)*(h+K2))
fm=f(6.0,4,8); assert sig(fm,3)==0.980
assert sig(9.8/fm,2)==10
assert sig(9.8/fm*f(4.0,4,8),2)==5.0
fq=f(6.5,6,7); assert sig(fq,2)==0.58
assert abs(fq-1/(1+10**-0.5)**2)<1e-12
# 最大の位置
import numpy as np
xs=np.linspace(5,8,30001); ys=[f(x,6,7) for x in xs]; assert abs(xs[int(np.argmax(ys))]-6.5)<1e-3
K1,K2=1e-6,1e-7; ff=fq/2
roots=sorted(-np.log10(np.roots([ff, ff*(K1+K2)-K1, ff*K1*K2])))
assert round(roots[0],1)==5.6 and round(roots[1],1)==7.4, roots
for x in (5.2,5.9,6.1): assert abs(f(x,6,7)-f(13-x,6,7))<1e-12
print('ALL OK')
