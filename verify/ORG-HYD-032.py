# 検証スクリプト（グループ g02）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-HYD-032/item.json'
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
        r = check_equation(b); assert r is None, (b, r); n_eq += 1; print('式OK:', b)
print('検査した式の数:', n_eq)
for k in ('src_no','src_label','field','theme','difficulty','time_min','format','problem','answer','solution','core','tags','verified','verify_note','problem_fig','solution_fig'):
    assert k in it, k
np_ = it['problem'].count('\\item'); na = len(re.findall(r'(^|\s)\((\d+)\)', it['answer']))
assert np_ == na, ('小問の数', np_, na)
from math import log10, floor
def sig(x, n): return round(x, -int(floor(log10(abs(x)))) + (n - 1))

from sympy import Rational as R
ca = R(100,100)/64; al = 3*R(100,100)/144
assert sig(float(ca),3) == 0.0156 and sig(float(al),3) == 0.0208 and al > ca
print('CaC2 →', float(ca), ' Al4C3 →', float(al))
# 陰イオンの電荷 = 受け取る H+ の数 = 炭化水素の H の数
for ion, q, hc in [('C2', 2, 'C2H2'), ('C', 4, 'CH4'), ('C3', 4, 'C3H4')]:
    assert int(re.search(r'H(\d)', hc).group(1)) == q
print('OK')
