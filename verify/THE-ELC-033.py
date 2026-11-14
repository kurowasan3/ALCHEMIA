# 検証スクリプト THE-ELC-033  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-ELC-033/item.json'
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
F, T, n = 96500, 298, 2
assert len(EQS) == 2
res = {}
for name, E25, E35 in (('A', R(985,1000), R(981,1000)), ('B', R(2041,1000), R(2043,1000))):
    dEdT = (E35 - E25)/10
    W = n*F*E25; Q = -n*F*T*dEdT; dH = -(W + Q)
    eta = W/(-dH); Ep = -dH/(n*F)
    assert Ep == E25 - T*dEdT
    print(name, float(dEdT), float(W)/1000, float(Q)/1000, float(dH)/1000, float(eta), float(Ep))
    res[name] = (W, Q, dH, eta, Ep, dEdT)
W,Q,dH,eta,Ep,d = res['A']
assert sig(float(W)/1000,3)==190 and sig(float(Q)/1000,3)==23.0 and sig(float(dH)/1000,3)==-213
assert sig(float(eta)*100,3)==89.2 and sig(float(Ep),3)==1.10 and Ep > R(985,1000)
W,Q,dH,eta,Ep,d = res['B']
assert sig(float(W)/1000,3)==394 and sig(float(Q)/1000,3)==-11.5 and sig(float(dH)/1000,3)==-382
assert sig(float(eta)*100,3)==103 and sig(float(Ep),3)==1.98 and Ep < R(2041,1000)
I = T*d/R(4,100); print('I*', float(I)); assert sig(float(I),3)==1.49
# I < I* なら正味の発熱 P<0（冷える）
for Iv in (R(1,2), R(1)): assert Iv**2*R(4,100) - Iv*T*d < 0
assert R(2)**2*R(4,100) - 2*T*d > 0

print('全検証OK:', ITEM)
