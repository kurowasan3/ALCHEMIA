# 検証スクリプト THE-SOL-038  python3 -I この.py
ITEM = '/home/claude/chem2/out/items/THE-SOL-038/item.json'
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
from sympy import symbols, solve, Eq
print('反応式', check_equations())
R = 8.3e3
k10, k40 = 1.19/22.4/1.0e5, 0.53/22.4/1.0e5   # mol/(L·Pa)（水 1.00 L）
P = symbols('P', positive=True)
res = {}
for Vg in (0.50, 5.0):
    nd = k10*2.0e5; ng = 2.0e5*Vg/(R*283); n = nd + ng
    P40 = float(solve(Eq(P*Vg/(R*313) + k40*P, n), P)[0])
    res[Vg] = (nd, ng, n, P40)
    print('空間', Vg, '溶解', nd, '気体', ng, '総量', n, 'P40', P40, '倍率', P40/2.0e5, '溶解割合', nd/n)
nd, ng, n, P40 = res[0.50]
assert sig(nd, 2) == 0.11 and sig(ng, 2) == 0.043 and sig(P40, 2) == 3.5e5
ch = 2.0e5*313/283; print('溶解無視', ch); assert sig(ch, 2) == 2.2e5 and P40 > ch
nd5, ng5, n5, P5 = res[5.0]
assert sig(P5, 2) == 2.5e5 and sig(P5/2.0e5, 2) == 1.2 and sig(P40/2.0e5, 2) == 1.7 and P5/2e5 < P40/2e5
print('40℃の溶解量', k40*P40, '放出', nd - k40*P40)
assert sig(k40*P40, 2) == 0.082 and sig(P40*0.5/(R*313), 2) == 0.067
assert round(nd5/n5*100) == 20 and round(nd/n*100) == 71
print('すべて一致')
