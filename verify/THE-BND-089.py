# 検証スクリプト（x1 g03）  python3 -I この.py
import sys, json, re
from math import log10, floor, sqrt, pi, exp, log
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chemx/x1/out/items/THE-BND-089/item.json'
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

d={n:c*z/(a*1e-7)**3 for n,a,z,c in [('Li',0.351,1,2),('Na',0.429,1,2),('K',0.533,1,2),('Al',0.405,3,4)]}
assert sig(d['Na'],2)==2.5e22 and sig(d['K'],2)==1.3e22 and sig(d['Al'],2)==1.8e23
assert sig(d['Al']/6.0e23*9.65e4,2)==2.9e4
mp={'Li':181,'Na':98,'K':64,'Al':660}
assert sorted(d,key=d.get)==sorted(mp,key=mp.get)
print('ALL OK')
