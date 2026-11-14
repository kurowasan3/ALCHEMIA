# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/POL-SUG-032/item.json'
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
        r = check_equation(b)
        assert r is None, (b, r)
        n_eq += 1
        print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from fractions import Fraction
AM = dict(H=1.0, C=12.0, N=14.0, O=16.0, Na=23.0, K=39.0, Cu=63.5, I=127.0)
from chemlib import parse_formula
def M(f):
    c, _ = parse_formula(f)
    return sum(AM[e] * float(n) for e, n in c.items())

import sympy as sp
x = sp.symbols('x')
u = lambda k: M('C6H10O5') + k*(M('NO2') - M('H'))
assert M('NO2') - M('H') == 45
assert sig(3*14/u(3)*100,3) == 14.1
xs = float(sp.solve(sp.Eq(14*x/(162+45*x), sp.Rational(1111,10000)), x)[0])
assert sig(xs,3) == 2.00 and sig(xs/3*100,3) == 66.7
n = 48.6/162
for xv in (xs, 2.0):
    assert sig(n*u(xv),3) == 75.6, n*u(xv)
    assert sig(n*xv*M('HNO3'),3) == 37.8, n*xv*M('HNO3')
print('POL-SUG-032 数値OK', xs)
