# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-MOL-049/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
from sympy import Rational as Q
Fe, O = 56, 16
n0 = Q(200, 1000)
m0 = n0 * (3*Fe + 4*O)
assert m0 == Q(464, 10)
mFe = 3*n0*Fe
assert mFe == Q(336, 10)
nH2_full = 3*n0*Q(4, 3)
assert nH2_full == Q(8, 10) and (m0 - mFe)/O == nH2_full
masses = {0: 46.4, 2.0: 42.4, 6.0: 37.6, 16.0: 34.2}
res = {}
for t in (2.0, 6.0, 16.0):
    h2 = (46.4 - masses[t]) / 16
    per_h = h2 * 60 / (t + 4.0)
    res[t] = per_h
    print(t, 'H2/cycle', round(h2, 4), 'H2/h', round(per_h, 4), sig(per_h, 2))
assert sig(res[2.0], 2) == 2.5 and sig(res[6.0], 2) == 3.3 and sig(res[16.0], 2) == 2.3
assert max(res, key=res.get) == 6.0
# 減少速度が単調に遅くなる（還元の進み方として妥当）
r1 = (46.4-42.4)/2; r2 = (42.4-37.6)/4; r3 = (37.6-34.2)/10
assert r1 > r2 > r3
# 34.4 g は完全還元(33.6 g)より大きい
assert masses[16.0] > float(mFe)
# 全体の反応 = 2式の和
from collections import Counter
print('OK')
