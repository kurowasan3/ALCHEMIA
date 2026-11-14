# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/outh/items/THE-GAS-068/item.json'
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

T=[273,283,293,298,303]; L=[0.928,0.861,0.773,0.710,0.593]; G=[0.098,0.135,0.194,0.242,0.335]
x=lambda r,i:(r-G[i])/(L[i]-G[i])
assert sig(x(0.30,0),2)==0.24 and sig(x(0.30,3),2)==0.12
assert sig(x(0.70,0),2)==0.73 and sig(x(0.70,3),2)==0.98
assert x(0.30,4)<0 and x(0.30,3)>0
assert x(0.70,4)>1 and x(0.70,3)<1
assert all(0.44<x(0.47,i)<0.53 for i in range(5))
assert all(L[i]>L[i+1] and G[i]<G[i+1] for i in range(4))
assert sig(0.47*10.0,2)==4.7
print('ALL OK')
