# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-026/item.json'
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
x = sp.symbols('x'); assert sp.solve(2 * x + 3 * (3 - x) - 8, x) == [1]
y = sp.symbols('y'); s = sp.solve(3 * y + 2 * (sp.Rational(94, 100) - y) - 2, y)[0]; print('Fe3+', s, float(s / sp.Rational(94, 100)))
assert s == sp.Rational(12, 100) and round(float(s / sp.Rational(94, 100)) * 100) == 13
vac = 4 - 4 * 0.94; print('空位', vac); assert sig(vac, 2) == 0.24
rho = 4 * (0.94 * 56 + 16) / (6.0e23 * 0.0795e-21); print('密度', rho); assert sig(rho, 2) == 5.8

print('全検証OK:', ITEM)
