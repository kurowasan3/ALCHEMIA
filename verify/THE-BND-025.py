# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-025/item.json'
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
NA = 6.0e23
n = sp.Rational(12, 6) + sp.Rational(2, 2) + 3; print('原子数', n); assert n == 6
V = 1.5 * 1.73 * 0.295 ** 2 * 0.468; V2 = 1.5 * 3 ** 0.5 * 0.295 ** 2 * 0.468; print('体積 nm3', V, V2)
assert sig(V * 1e-21, 2) == 1.1e-22 and sig(V2 * 1e-21, 2) == 1.1e-22
rho = 6 * 47.9 / (NA * V * 1e-21); rho_r = 6 * 47.9 / (NA * 1.06e-22); print('密度', rho, rho_r)
assert sig(rho, 2) == 4.5 and sig(rho_r, 2) == 4.5
a = sp.symbols('a', positive=True); r = a / 2; c = sp.sqrt(sp.Rational(8, 3)) * a
phcp = sp.simplify(6 * sp.Rational(4, 3) * sp.pi * r ** 3 / (3 * sp.sqrt(3) / 2 * a ** 2 * c))
pfcc = sp.sqrt(2) * sp.pi / 6
print('hcp', phcp, 'fcc', pfcc); assert sp.simplify(phcp - pfcc) == 0

print('全検証OK:', ITEM)
