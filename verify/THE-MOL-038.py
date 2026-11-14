# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/out/items/THE-MOL-038/item.json'
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
nO = sp.Rational(240, 10) / 80
x = sp.symbols('x'); xs = sp.solve(x + 4*x - nO, x)[0]
nC = xs + 2*xs
print('C mol', nC, 'g', float(nC*12), 'L', float(nC*sp.Rational(224, 10)))
assert sig(float(nC*12), 2) == 2.2 and sig(float(nC*sp.Rational(224, 10)), 2) == 4.0
a, b = sp.symbols('a b')
s = sp.solve([a + b - sp.Rational(288, 100)/12, a + 2*b - nO], [a, b])
assert s[a] / s[b] == 3
lo, hi = nO/2*12, nO*12
print('範囲', float(lo), float(hi)); assert (float(lo), float(hi)) == (1.8, 3.6) and 1.5 < lo
for e in ['2CuO + C -> 2Cu + CO2', 'CuO + C -> Cu + CO']:
    assert check_equation(e) is None
print('OK')
