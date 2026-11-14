# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-MOL-016/item.json'
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

nH2O2 = 50.0 * 0.034 / 34; nO2 = nH2O2 / 2; print('O2', nO2); assert abs(nO2 - 0.025) < 1e-12
z = nO2 * 0.060; nO3 = 2 * z; print('O3', nO3); assert sig(nO3, 2) == 0.0030
tot = nO2 - z; Mavg = nO2 * 32 / tot; print('平均分子量', Mavg); assert sig(Mavg, 2) == 34
after = (nO2 - 3 * z) + nO3 - nO3 + nO3   # O3 1 mol → O2 1 mol
assert abs(after - tot) < 1e-12

print('全検証OK:', ITEM)
