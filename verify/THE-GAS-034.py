# 検証スクリプト THE-GAS-034  python3 -I この.py
ITEM = '/home/claude/chem2/outh/items/THE-GAS-034/item.json'
import sys, json, re
from math import log10, floor
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
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
def check_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        print('式', e, '->', 'OK' if r is None else r)
        assert r is None, (e, r)
    return sorted(eqs)
def sig(x, n):
    return float(f'{x:.{n}g}')
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
def close(a, b, rel=1e-9):
    assert abs(a - b) <= rel * max(abs(a), abs(b)), (a, b)
check_structure()
from sympy import Rational as Q, symbols, solve, nsimplify
print('反応式', check_equations())
R = Q(83, 10)*1000; TA, TB = 310, 300; VA, VB = Q(2), Q(1, 2)
PA, PB = Q(63, 10)*1000, Q(36, 10)*1000
cA, cB = VA/(R*TA), VB/(R*TB)
# (1)
maxA = PA*cA; print('(1) A に気体で存在できる最大量', float(maxA)); assert maxA < Q(1, 10)
# 平衡状態の候補を総当たりで調べる：液体が A のみ / B のみ / 両方 / なし
cands = {}
for liqA in (0, 1):
    for liqB in (0, 1):
        if liqA and liqB:
            ok = (PA == PB)                      # 圧力一様と矛盾
        elif liqA:
            ok = PA <= PB                        # B で凝縮しない条件
        elif liqB:
            ok = PB <= PA and Q(1, 10) >= PB*(cA+cB)
        else:
            P = Q(1, 10)/(cA+cB); ok = P <= PB and P <= PA
        cands[(liqA, liqB)] = ok
print('(1) 平衡になりうる状態', [k for k, v in cands.items() if v]); assert [k for k, v in cands.items() if v] == [(0, 1)]
vap = PB*(cA+cB); liq = Q(1, 10) - vap
print('(1) 気体の水', float(vap), '液体', float(liq)); assert sig(float(liq), 3) == 0.0965
# (2)
x = PB/PA; assert x == Q(4, 7)
x0 = Q(150, 180); print('(2) 初めの A の水蒸気圧', float(x0*PA)); assert x0*PA > PB
n = symbols('n', positive=True)
nA = solve(n/(n+Q(3, 10)) - x, n)[0]; print('(2) A の水', nA, 'B の水', Q(3, 2)-nA)
assert nA == Q(2, 5) and sig(float(Q(3, 2)-nA), 3) == 1.10
# (3)
x0 = Q(1, 2); assert x0*PA < PB      # B から A へ
moved = nA - Q(3, 10); B = 1 - moved
print('(3) 移動', moved, 'B の水', B, '熱量 kJ', float(moved*44)); assert B == Q(9, 10) and sig(float(moved*44), 2) == 4.4
# (4)
nA4 = Q(3, 10) + Q(5, 100); x4 = nA4/(nA4+Q(3, 10)); P4 = x4*PA
print('(4) x', float(x4), 'P', float(P4)); assert P4 < PB and sig(float(P4), 2) == 3400
print('(4) 必要な B の水', moved); assert moved == Q(1, 10)
# 気相の水を無視した近似の妥当性（(2)-(4)）
print('気相の水', float(vap), 'mol（液体 0.3〜1.5 mol に対して 1 % 前後）')
print('すべて一致')
