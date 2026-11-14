# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/out/items/THE-MOL-036/item.json'
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

n_eq = 0
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if re.search(r'boxed|\\text|\\mathrm', b): continue
        if len(re.findall(r'<=>|->|<-', b)) != 1: continue
        if ' + ' not in b and not re.search(r'(^|\s)\d+\s*[A-Z(\[]', b): continue
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import sympy as sp
nCO2 = sp.Rational(440, 100) / 44; nH2O = sp.Rational(252, 100) / 18
ntot = nH2O - nCO2; nO2 = sp.Rational(3, 2) * nCO2
print('総物質量', float(ntot), 'O2', float(nO2), 'L', float(nO2 * sp.Rational(224, 10)))
assert ntot == sp.Rational(4, 100) and float(nO2 * sp.Rational(224, 10)) == 3.36
nbar = nCO2 / ntot; assert nbar == sp.Rational(5, 2)
a, b = sp.symbols('a b')
sol = sp.solve([a + b - ntot, 2*a + 3*b - nCO2], [a, b])
assert sol[a] == sol[b] == sp.Rational(2, 100)
mass = sol[a]*46 + sol[b]*60; assert mass == sp.Rational(212, 100)
assert sol[a]*3 + sol[b]*4 == nH2O
# O 原子の収支
assert ntot + 2*nO2 == 2*nCO2 + nH2O
for e in ['C2H5OH + 3O2 -> 2CO2 + 3H2O', '2C3H7OH + 9O2 -> 6CO2 + 8H2O', 'C2H5OC2H5 + 6O2 -> 4CO2 + 5H2O', 'C3H5OH + 4O2 -> 3CO2 + 3H2O']:
    assert check_equation(e) is None
print('OK')
