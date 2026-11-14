# 検証スクリプト（グループ g04）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem5/out/items/ORG-ARO-055/item.json'
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

# ---- 計算（NaOH の物質量は反応式の係数から）
from sympy import Rational as R
def M(f):
    from chemlib import parse_formula
    w = {'H': R(1), 'C': R(12), 'O': R(16), 'Na': R(23), 'S': R(32), 'Cl': R(355, 10)}
    return sum(w[e]*n for e, n in parse_formula(f)[0].items())
assert M('C6H6O') == 94 and M('NaOH') == 40 and M('Na2SO3') == 126 and M('C3H6O') == 58
n = R(9400)/94; assert n == 100
naoh_i = 1 + 2; naoh_u = 2
for e in ['C6H5SO3H + NaOH -> C6H5SO3Na + H2O', 'C6H5SO3Na + 2NaOH -> C6H5ONa + Na2SO3 + H2O', 'C6H5Cl + 2NaOH -> C6H5ONa + NaCl + H2O']:
    assert check_equation(e) is None
print('NaOH kg', float(n*naoh_i*40/1000), 'Na2SO3 kg', float(n*126/1000), 'アセトン kg', float(n*58/1000))
assert n*naoh_i*40 == 12000 and n*126 == 12600 and n*58 == 5800
print('OK')
