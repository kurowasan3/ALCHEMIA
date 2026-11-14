# 検証スクリプト INO-ALK-017  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/INO-ALK-017/item.json'
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
it = json.load(open(ITEM, encoding='utf-8'))
def ce_blocks(s):
    return re.findall(r'\\ce\{((?:[^{}]|\{[^{}]*\})*)\}', s)
def clean_ce(eq):
    eq = re.sub(r'(->|<=>|<-)\[[^\]]*\](\[[^\]]*\])?', r'\1', eq)
    eq = re.sub(r'\\[a-zA-Z]+\{[^{}]*\}|\\[,;:! ]|\\quad', ' ', eq)
    eq = re.sub(r'\{(\d*[+-])\}', r'\1', eq)
    return eq
def should_check(eq):
    # 矢印が1つの完全な式はすべて（validate ツールが一般式とみなして飛ばす Zn・Sn なども含めて）検査する．
    # 一般式・同位体の式は矢印に空ラベル [\\text{}] を付けてあるので除く．
    if re.search(r'boxed|□|\?|？|…|\\square', eq): return False
    if re.search(r'(->|<=>|<-)\[\\text\{\}\]', eq): return False
    if len(re.findall(r'<=>|->|<-', eq)) != 1: return False
    return ' + ' in eq or re.search(r'(^|\s)\d+(/\d+)?\s*[A-Z(\[]', eq) is not None
def check_all_equations():
    eqs = sorted({b for k in ('problem','answer','solution') for b in ce_blocks(it[k]) if should_check(b)})
    for e in eqs:
        r = check_equation(clean_ce(e)); print('式', e, '→', r or 'OK'); assert r is None, e
    return eqs
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    nums = [int(x) for x in re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])]
    assert nums == list(range(1, n_items+1)), nums
    print('小問数', n_items)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
check_structure()
EQS = check_all_equations()

# ---- 個別の検証 ----

from sympy import Rational as R, symbols, solve
need = ['NaCl + NH3 + CO2 + H2O -> NaHCO3 + NH4Cl','2NaHCO3 -> Na2CO3 + H2O + CO2','CaCO3 -> CaO + CO2','CaO + H2O -> Ca(OH)2',
        '2NH4Cl + Ca(OH)2 -> CaCl2 + 2H2O + 2NH3','NH4HCO3 -> NH3 + H2O + CO2','2NaCl + CaCO3 -> Na2CO3 + CaCl2','NH3 + CO2 + H2O -> NH4+ + HCO3-']
for e in need: assert e in EQS, e
x, b = symbols('x b')
def prec(bv, K): return min(s for s in solve((5-x)*(bv-x) - K, x) if s.is_real and 0 <= s <= 5)
x1 = prec(5, 1); print(x1); assert x1 == 4
bs = solve(R(1,2)*(b - R(45,10)) - 1, b)[0]; assert bs == R(13,2)
x3 = prec(5, R(64,100)); print(x3); assert x3 == R(42,10)
b95 = solve(R(1,4)*(b - R(475,100)) - 1, b)[0]; assert b95 == R(875,100)
assert sig(float(prec(6, 1)/5*100),3) == 87.6  # 参考
xb = R(45,10); bb = R(13,2)
Na, NH4, Cl, HCO3 = 5 - xb, bb, 5, bb - xb
assert Na + NH4 == Cl + HCO3
nh4hco3 = HCO3; nh4cl = Cl - Na; assert nh4hco3 == 2 and nh4cl == xb
caoh2 = nh4cl/2; assert caoh2 == R(9,4)
co2_back = xb/2 + nh4hco3 + caoh2; assert co2_back == bb
assert nh4hco3 + nh4cl == bb
t = R(10**6,106)*5/(xb/2)*R(585,10)/10**6; print('t', float(t)); assert sig(float(t),3) == 1.23
ideal = R(2*585,10)/106; assert sig(float(ideal),3) == 1.10
assert 5*bb > 30

print('全検証OK:', ITEM)
