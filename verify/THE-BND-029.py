# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-029/item.json'
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
assert 1 + sp.Rational(12, 4) == 4
R = sp.symbols('R', positive=True); a = 2 * sp.sqrt(2) * R
r_oct = sp.simplify(a / 2 - R); r_tet = sp.simplify(sp.sqrt(3) / 4 * a - R)
print(r_oct / R, float(r_oct / R), r_tet / R, float(r_tet / R))
assert sp.simplify(r_oct / R - (sp.sqrt(2) - 1)) == 0 and sp.simplify(r_tet / R - (sp.sqrt(6) / 2 - 1)) == 0
assert round(1.414 - 1, 3) == 0.414 and round(2.4495 / 2 - 1, 3) == 0.225 and sig(float(r_tet / R), 3) == 0.225
# 組成
assert 8 / 4 == 2                       # M2X
frac = sp.Rational(8, 3) / 4; print('占有率', frac); assert frac == sp.Rational(2, 3)
assert 2 * 4 + 4 * (-2) == 0            # MO2

print('全検証OK:', ITEM)
