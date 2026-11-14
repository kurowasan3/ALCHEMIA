# 検証スクリプト（x2 g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x2/out/items/THE-ION-112/item.json'
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

nprob = it['problem'].count('\\item')
nans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert nprob == nans or (nprob == 0 and nans == 1), (nprob, nans)
print('小問数', nprob, '答え', nans)

from math import log10
def frac(A,a,b): return (a-A)/(a-b)
x=0.5; A2=0.64*(1-x)+0.12*x; print('(2)',A2); assert abs(A2-0.38)<1e-9
x3=frac(0.25,0.64,0.12); pH3=4.0+log10(x3/(1-x3)); print('(3)',x3,pH3); assert abs(x3-0.75)<1e-9 and round(pH3,1)==4.5
s=0.28/0.35; x4=frac(0.408,0.64*s,0.12*s); pH4=4.0+log10(x4/(1-x4)); print('(4)',s,x4,pH4); assert abs(x4-0.25)<1e-9 and round(pH4,1)==3.5
import math
a,w1,b,w2=0.651,math.sqrt(5800),0.385,math.sqrt(4204)
acid=lambda L:a*math.exp(-((L-520)/w1)**2); base=lambda L:b*math.exp(-((L-440)/w2)**2)
assert abs(acid(460)-0.35)<0.005 and abs(base(460)-0.35)<0.005 and abs(acid(510)-0.64)<0.005 and abs(base(510)-0.12)<0.005
svg=open(ITEM.replace('item.json','problem.svg')).read(); assert 'viewBox' in svg and 'script' not in svg

print('ALL OK')
