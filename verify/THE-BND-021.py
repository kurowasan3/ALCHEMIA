# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-021/item.json'
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
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
NA = 6.0e23
print('C60 MW', MolWt(Chem.MolFromSmiles('[C]' * 1)) * 60)
assert 4 * 60 == 240
d = 1.42 / 1.41; d2 = 1.42 / sqrt(2); print('距離', d, d2); assert sig(d, 2) == 1.0 and sig(d2, 2) == 1.0
rho = 4 * 720 / (NA * 2.86e-21); print('密度', rho); assert sig(rho, 2) == 1.7
x = (4.19 / 3.60 - 1) * 720 / 39; print('x', x); assert round(x) == 3 and abs(x - 3) < 0.1
assert 4 * round(x) == 12

print('全検証OK:', ITEM)
