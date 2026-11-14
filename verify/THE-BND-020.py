# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-020/item.json'
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
a = 2 * (0.133 + 0.181); print('a', a); assert abs(a - 0.628) < 1e-9
dK = a / sqrt(2); dKa = 1.41 * a / 2; print('K-K', dK, dKa); assert sig(dK, 2) == 0.44 and sig(dKa, 2) == 0.44
rhoK = 4 * 74.5 / (NA * 0.248e-21); rhoK2 = 4 * 74.5 / (NA * (a * 1e-7) ** 3); print('KCl', rhoK, rhoK2)
assert sig(rhoK, 2) == 2.0 and sig(rhoK2, 2) == 2.0
aN = 2 * (0.101 + 0.181); rhoN = 4 * 58.5 / (NA * 0.179e-21); print('NaCl a', aN, rhoN)
assert abs(aN - 0.564) < 1e-9 and sig(rhoN, 2) == 2.2
print('質量比', 74.5 / 58.5, '体積比', 0.248 / 0.179)
assert 0.248 / 0.179 > 74.5 / 58.5

print('全検証OK:', ITEM)
