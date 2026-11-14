# 検証スクリプト INO-ALK-015  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/INO-ALK-015/item.json'
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

from sympy import Rational as R
assert 'Na2CO3.10H2O -> Na2CO3.H2O + 9H2O' in EQS and '2NaOH + CO2 -> Na2CO3 + H2O' in EQS and 'Na2CO3 + 2HCl -> 2NaCl + H2O + CO2' in EQS
m10, m1, m0 = 106+180, 106+18, 106
assert R(286,10)/m10 == R(1,10) and R(1,10)*m1 == R(124,10)
x = (R(543,100)/18)/(R(457,100)/106); print('x', float(x)); assert sig(float(x),2) == 7.0
m7 = 106 + 7*18; assert m7 == 232
w = R(5,100)*m10; assert w == R(143,10)
n = w/m7; c = n/R(1,2); print('c', float(c)); assert sig(float(c),3) == 0.123
nHCl = 2*c*R(10,1000); V = nHCl/R(1,10); print('V mL', float(V*1000))
cs = 2*R(1,10)*R(10,1000)/V; print('生徒の値', float(cs)); assert sig(float(cs),3) == 0.0811 and cs < R(1,10)
assert cs == R(1,10)*R(1,10)/c
reac = json.load(open('/home/claude/kurowasan3/alchemia/content/chem/reactions.json', encoding='utf-8'))
rs = reac if isinstance(reac, list) else reac.get('reactions', reac.get('items', []))
eqs = [r.get('eq') for r in rs if isinstance(r, dict)]
assert 'Na2CO3.10H2O -> Na2CO3.H2O + 9H2O' in eqs and '2NaOH + CO2 -> Na2CO3 + H2O' in eqs

print('全検証OK:', ITEM)
