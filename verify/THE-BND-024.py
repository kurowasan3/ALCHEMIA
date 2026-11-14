# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-024/item.json'
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

from math import sqrt
NA = 6.0e23
assert 6 * 1 / 3 == 2
S = 1.5 * 1.73 * 0.0202; S2 = 1.5 * sqrt(3) * 0.142 ** 2; print('面積', S, S2); assert sig(S, 2) == 0.052 and sig(S2, 2) == 0.052
v = S / 2 * 0.335 * 1e-21; rho = 12 / (NA * v); rho2 = 12 / (NA * S2 / 2 * 0.335e-21); print('密度', rho, rho2)
assert sig(rho, 2) == 2.3 and sig(rho2, 2) == 2.3
assert rho < 3.5

print('全検証OK:', ITEM)
