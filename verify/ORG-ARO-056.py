# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ARO-056/item.json'
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

# ---- 計算
from sympy import Rational as R, symbols, solve

# ---- 計算
from chemlib import parse_formula
from sympy import Rational as R
w = {'H': R(1), 'C': R(12), 'O': R(16), 'Na': R(23), 'Br': R(80)}
def M(f): return sum(w[e]*n for e, n in parse_formula(f)[0].items())
assert M('C6H5ONa') == 116 and M('C6H2Br3OH') == 331 and M('C6H4(OH)COOH') == 138
n0 = R(232,100)/116
nph = R(430,100)/331
conv = (n0 - nph)/n0*100
mB = (n0 - nph)*138
print('変化率 %', float(conv), ' B g', float(mB))
assert round(float(conv), 1) == 35.0 and round(float(mB), 3) == 0.967
# 途中で 0.00701 に丸めても同じ値か
assert round(0.00701/0.0200*100, 1) == 35.0 and round(0.00701*138, 3) == 0.967
print('OK')
