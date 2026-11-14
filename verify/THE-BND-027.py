# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-027/item.json'
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

assert 1 + 12 / 4 == 4
rho = 4 * 106 / (6.0e23 * 0.0589e-21); rho2 = 4 * 106 / (6.0e23 * (0.389e-7) ** 3); print('密度', rho, rho2)
assert sig(rho, 3) == 12.0 and sig(rho2, 3) == 12.0
x = 4 * 0.60 / 4; assert abs(x - 0.60) < 1e-12
V = 12.0 / 106 * 0.60 / 2 * 22.4e3; print('H2 cm3', V); assert sig(V, 2) == 760
assert 8 / 4 == 2

print('全検証OK:', ITEM)
