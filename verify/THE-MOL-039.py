ITEM = '/home/claude/chem4/out/items/THE-MOL-039/item.json'
# 検証スクリプト（グループ g07）  python3 -I この.py
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

# 小問数と答えの数の一致
n_items = it['problem'].count('\\item')
n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
assert n_items == n_ans or (n_items == 0 and n_ans == 0), (n_items, n_ans)
print('小問数', n_items, '答え', n_ans)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import Rational as R
x = R(500, 100); f = 25
M_CaCl2 = 40 + 2*R(355, 10); M_HCl = R(365, 10); M_H2O = 18
nCa = R(500, 10000) * R(240, 10000) * f
nCl = R(100, 1000) * R(248, 10000) * f
print('nCa', float(nCa), 'nCl', float(nCl))
z1 = (x - nCa*M_CaCl2) / x * 100
z2 = (x - nCl/2*M_CaCl2) / x * 100
mHCl = (nCl - 2*nCa) * M_HCl
w = x - nCa*M_CaCl2 - mHCl
z3 = w / x * 100
hyd = w / M_H2O / nCa
for name, v, n in [('z1', z1, 3), ('z2', z2, 3), ('mHCl', mHCl, 3), ('z3', z3, 3), ('水和数', hyd, 3)]:
    print(name, float(v), sig(float(v), n))
assert sig(float(z1), 3) == 33.4 and sig(float(z2), 3) == 31.2 and sig(float(mHCl), 3) == 0.073 and sig(float(z3), 3) == 31.9
assert sig(float(hyd), 3) == 2.96
assert z1 > z3 > z2
# 丸めの境界でないこと
for v in (z1, z2, z3):
    t = float(v) * 10 % 1
    assert abs(t - 0.5) > 0.02, v
print('すべて OK')
