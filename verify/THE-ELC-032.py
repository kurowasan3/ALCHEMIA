# 検証スクリプト THE-ELC-032  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-ELC-032/item.json'
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
assert '2H2O + 2e- -> H2 + 2OH-' in EQS and 'Cl2 + 2OH- -> Cl- + ClO- + H2O' in EQS and '2Cl- -> Cl2 + 2e-' in EQS
F = 96500
nNaOH = R(100*1000*32, 100) / 40; assert nNaOH == 800
I = nNaOH*F/3600; print('I', float(I)); assert sig(float(I)/1000, 3) == 21.4
feed = (100 - 32) + nNaOH*18/1000; print('feed', float(feed)); assert feed == R(824,10)
x = symbols('x'); s = solve((R(26,100)*x - R(468,10)) - R(20,100)*(x - R(468,10)), x)[0]
print('brine', s); assert s == 624
assert (R(26,100)*s - R(468,10))/(s - R(468,10)) == R(20,100)
assert nNaOH*R(585,10)/1000 == R(468,10)
feed40 = 32/R(40,100) - 32 + nNaOH*18/1000; print('feed40', float(feed40)); assert feed40 == R(624,10)
import os; assert os.path.exists(os.path.join(os.path.dirname(ITEM), it['problem_fig']))

print('全検証OK:', ITEM)
