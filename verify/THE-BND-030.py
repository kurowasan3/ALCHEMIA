# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-030/item.json'
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
lim = 1.73 - 1; print('CsCl 限界', lim, sqrt(3) - 1); assert sig(lim, 3) == 0.730
r = {'Cs': 0.174, 'Na': 0.102, 'Cu': 0.060, 'Br': 0.196, 'Cl': 0.181}
def typ(x):
    return 'CsCl' if x >= sqrt(3) - 1 else 'NaCl' if x >= sqrt(2) - 1 else 'ZnS' if x >= sqrt(6) / 2 - 1 else '?'
cases = {'CsBr': r['Cs'] / r['Br'], 'NaBr': r['Na'] / r['Br'], 'CuCl': r['Cu'] / r['Cl']}
for k, v in cases.items(): print(k, round(v, 3), typ(v))
assert [round(v, 3) for v in cases.values()] == [0.888, 0.520, 0.331]
assert [typ(v) for v in cases.values()] == ['CsCl', 'NaCl', 'ZnS']
a = 2 * (0.174 + 0.196) / 1.73; a2 = 2 * 0.370 / sqrt(3); print('a', a, a2); assert sig(a, 2) == 0.43 and sig(a2, 2) == 0.43
print('NaBr 仮想 CsCl 型 a=2r-', 2 * r['Br'])

print('全検証OK:', ITEM)
