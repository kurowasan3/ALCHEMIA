# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-ATM-011/item.json'
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

from fractions import Fraction as F
from collections import defaultdict
import sympy as sp
Br = {79: F(1, 2), 81: F(1, 2)}; Cl = {35: F(3, 4), 37: F(1, 4)}
def dist(A, B):
    d = defaultdict(F)
    for m1, p1 in A.items():
        for m2, p2 in B.items():
            d[m1 + m2] += p1 * p2
    return dict(sorted(d.items()))
d1 = dist(Br, Br); d2 = dist(Br, Cl); print(d1, d2)
assert list(d1) == [158, 160, 162] and [v * 4 for v in d1.values()] == [1, 2, 1]
assert list(d2) == [114, 116, 118] and [v * 8 for v in d2.values()] == [3, 4, 1]
avg = sum(m * p for m, p in d2.items()); print('平均', float(avg)); assert avg == F(231, 2)
assert float(sum(m * p for m, p in Br.items()) + sum(m * p for m, p in Cl.items())) == 115.5
p = sp.symbols('p', positive=True)
sol = [s for s in sp.solve(sp.Eq(p ** 2 / (2 * p * (1 - p)), sp.Rational(9, 6)), p) if 0 < s < 1]
print(sol); assert sol == [sp.Rational(3, 4)]
assert sp.simplify((2 * p * (1 - p)) / (1 - p) ** 2 - 6).subs(p, sp.Rational(3, 4)) == 0

print('全検証OK:', ITEM)
