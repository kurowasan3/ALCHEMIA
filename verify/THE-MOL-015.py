# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-MOL-015/item.json'
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

import sympy as sp
x, y = sp.symbols('x y')
s = sp.solve([24 * x + 27 * y - sp.Rational(228, 100), x + sp.Rational(3, 2) * y - sp.Rational(110, 1000)], [x, y])
print({k: float(v) for k, v in s.items()})
assert s[x] == sp.Rational(1, 20) and s[y] == sp.Rational(1, 25)
pct = 24 * 0.05 / 2.28 * 100; print('Mg%', pct); assert round(pct) == 53
V = 1.5 * 0.040 * 22.4; print('H2(NaOH)', V); assert sig(V, 3) == 1.34

print('全検証OK:', ITEM)
