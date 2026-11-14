# 検証スクリプト THE-ELC-031  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-ELC-031/item.json'
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
F = 96500
assert 'Zn + Cu^2+ -> Zn^2+ + Cu' in EQS and 'Zn + Ag2O -> ZnO + 2Ag' in EQS
assert len(EQS) >= 7, len(EQS)
dH1 = -365 - (-146); assert dH1 == -219
E1 = R(-dH1*1000, 2*F); print('E1', float(E1)); assert sig(float(E1),3) == 1.13
dH2 = -350 - (-31); assert dH2 == -319
E2 = R(-dH2*1000, 2*F); print('E2', float(E2)); assert sig(float(E2),3) == 1.65
E1x2 = R(-2*dH1*1000, 4*F); assert E1x2 == E1
ratio = R(319,219); print('ratio', float(ratio), float(R(165,113))); assert sig(float(ratio),2) == 1.5 and sig(1.65/1.13,2) == 1.5
# 元の標準生成エンタルピーとの整合（丸め前）
Zn2, Cu2, Ag1 = -153.9, 64.8, 105.6
assert abs((Zn2 - 2*Ag1) - (-365)) < 0.5 and abs((Cu2 - 2*Ag1) - (-146)) < 0.5
assert abs((-350.5 + 31.1) - (-319)) < 0.5
for k in ('problem','answer','solution'):
    assert '\\boldsymbol{約' not in it[k]

print('全検証OK:', ITEM)
