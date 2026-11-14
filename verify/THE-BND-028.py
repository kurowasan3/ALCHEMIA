# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-BND-028/item.json'
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
n = (8 / 8, 12 / 4, 1); print(n); assert n == (1, 3, 1) and 2 * 1 + 4 * 1 - 2 * 3 == 0
a = 2 * (0.061 + 0.140); print('a', a); assert sig(a, 3) == 0.402
rmax = a / 1.41 - 0.140; rmax2 = a / sqrt(2) - 0.140; print('rmax', rmax, rmax2)
assert sig(rmax, 3) == 0.145 and 0.140 < rmax2 < rmax
aBa = 1.41 * (0.161 + 0.140); aBa2 = sqrt(2) * 0.301; print('a(Ba)', aBa, aBa2)
assert sig(aBa, 3) == 0.424
assert 0.161 > rmax


print('全検証OK:', ITEM)
