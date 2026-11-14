# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-022/item.json'
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
assert 8 * 1 / 8 + 6 / 2 + 4 == 8
bSi = 1.73 * 0.543 / 4; bSi2 = sqrt(3) * 0.543 / 4; print('Si-Si', bSi, bSi2); assert sig(bSi, 3) == 0.235 and sig(bSi2, 3) == 0.235
rho = 8 * 28.1 / (NA * 0.160e-21); print('密度', rho); assert sig(rho, 2) == 2.3
bC = 1.73 * 0.357 / 4; bC2 = sqrt(3) * 0.357 / 4; print('C-C', bC, bC2); assert sig(bC, 3) == 0.154
ref = {'エチレン': 0.134, 'エタン': 0.154, 'ベンゼン': 0.140, '黒鉛': 0.142}
best = min(ref, key=lambda k: abs(ref[k] - bC)); print(best); assert best == 'エタン'

print('全検証OK:', ITEM)
