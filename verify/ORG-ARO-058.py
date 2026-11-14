# 検証スクリプト（グループ g05）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ARO-058/item.json'
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
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
print('小問の数', np_)

from math import log10, floor, sqrt
def sig(x, n):
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

Kw = 1.0e-14; KbN = 1.8e-5; KbA = 4.0e-10; KbB = 2.2e-5
assert KbB > KbN > KbA
Ka_ = KbN / KbA; Kb_ = KbA / KbN
print('(a) K', Ka_, '(b) K', Kb_)
assert sig(Ka_, 2) == 4.5e4 and sig(Kb_, 2) == 2.2e-5
K3 = KbN / KbB
from sympy import symbols, solve, Rational
x = symbols('x', positive=True)
sol = [s for s in solve(x**2 - K3*(1-x)**2, x) if 0 < s < 1][0]
print('(3) K', K3, 'x', sol)
assert round(float(sol)*100) == 47
# (4) pH：近似と厳密（二次方程式）
def pH_salt(c, Kb):
    Ka = Kw / Kb
    h = (-Ka + sqrt(Ka*Ka + 4*Ka*c)) / 2
    return -log10(h), -log10(sqrt(c*Ka))
pa = pH_salt(0.10, KbA); pn = pH_salt(0.10, KbN)
print('アニリン塩酸塩 pH 厳密/近似', pa, ' NH4Cl', pn)
assert round(pa[1], 1) == 2.8 and round(pa[0], 1) == 2.8
assert round(pn[1], 1) == 5.1 and round(pn[0], 1) == 5.1
# log 値の与え方で計算した値
assert abs(0.5*(6 - (1 - 2*0.30)) - 2.8) < 1e-9
assert abs(0.5*(10 + (0.30 + 2*0.48 - 1)) - 5.13) < 1e-9
print('ORG-ARO-058 OK')
