# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem4/out/items/THE-MOL-034/item.json'
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
x = sp.symbols('x')
M = sp.Rational(125, 100) * sp.Rational(224, 10)
assert M == 28
s = sp.solve(16*x + 44*(1-x) - 28, x)[0]
assert s == sp.Rational(4, 7)
print('CH4:C3H8 =', s, ':', 1-s, '→ 4:3')
b = sp.solve(16*x + 44*(1-x) - sp.Rational(288, 10), x)[0]
print('境界', float(b)); assert sig(float(b)*100, 2) == 54
pairs = {'ア': (2, 4), 'イ': (4, 32), 'ウ': (28, 28), 'エ': (44, 40), 'オ': (16, 32)}
ok = [k for k, (a, c) in pairs.items() if min(a, c) <= 28 <= max(a, c)]
assert ok == ['イ', 'ウ', 'オ'], ok
print('OK')
