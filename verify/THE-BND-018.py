# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-018/item.json'
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
pf = sp.sqrt(2) * sp.pi / 6; pb = sp.sqrt(3) * sp.pi / 8
ratio = sp.nsimplify(pb / pf)
print('fcc', float(pf), 'bcc', float(pb), 'ratio', ratio, float(ratio))
assert sp.simplify(ratio - 3 * sp.sqrt(3) / (4 * sp.sqrt(2))) == 0
pfa = 1.41 * 3.14 / 6; pba = 1.73 * 3.14 / 8
print('近似', pfa, pba, pba / pfa, pfa / pba)
assert round(pfa * 100) == 74 and round(pba * 100) == 68
assert sig(float(ratio), 2) == 0.92 and sig(pba / pfa, 2) == 0.92
assert sig(1 / float(ratio), 2) == 1.1 and sig(pfa / pba, 2) == 1.1
r = sp.symbols('r', positive=True)
af = 2 * sp.sqrt(2) * r; ab = 4 * r / sp.sqrt(3)
assert sp.simplify(sp.sqrt(2) * af - 4 * r) == 0 and sp.simplify(sp.sqrt(3) * ab - 4 * r) == 0

print('全検証OK:', ITEM)
