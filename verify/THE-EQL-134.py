# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/out/items/THE-EQL-134/item.json'
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

K1=1/2; K2=1/0.25; s=1+K1+K2
assert sig(11/s,2)==2.0 and sig(11/s*K2,2)==8.0
a=19/s; assert sig(a,2)==3.5 and sig(a*K1,2)==1.7 and sig(a*K2,2)==14
# 数値積分で確認（k1=1, k-1=2, k2=1, k-2=0.25）
A,B,C=11.0,0.0,0.0; dt=1e-3
for _ in range(40000):
    r1=A-2*B; r2=A-0.25*C
    A,B,C=A-dt*(r1+r2),B+dt*r1,C+dt*r2
assert abs(A-2)<1e-6 and abs(C-8)<1e-6
# 触媒（正逆とも10倍）でも同じ
A,B,C=11.0,0.0,0.0; dt=1e-4
for _ in range(40000):
    r1=A-2*B; r2=10*(A-0.25*C)
    A,B,C=A-dt*(r1+r2),B+dt*r1,C+dt*r2
assert abs(C-8)<1e-4, C
print('ALL OK')
