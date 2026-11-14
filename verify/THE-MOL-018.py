# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-MOL-018/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

nAgCl = 2.65 / 143.5; mCl = nAgCl * 35.5; mM = 1.00 - mCl; per = mM / nAgCl
print(nAgCl, mCl, mM, per)
assert sig(nAgCl, 3) == 0.0185 and sig(mM, 3) == 0.344 and sig(per, 3) == 18.7
cand = {'Na': 23, 'Mg': 24, 'Al': 27, 'K': 39, 'Ca': 40, 'Fe': 56, 'Cu': 63.5, 'Zn': 65.4}
best = min(((n, s, abs(per * n - m) / m) for n in (1, 2, 3) for s, m in cand.items()), key=lambda t: t[2])
print(best); assert best[:2] == (3, 'Fe') and best[2] < 0.01
print('FeCl3 → AgCl', 1.00 / 162.5 * 3 * 143.5)
m2 = 1.00 / 127 * 2 * 143.5; print('FeCl2 → AgCl', m2); assert sig(m2, 3) == 2.26

print('全検証OK:', ITEM)
