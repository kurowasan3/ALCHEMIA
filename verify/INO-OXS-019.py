# 検証スクリプト INO-OXS-019  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/INO-OXS-019/item.json'
import sys, json, re
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
def check_all_equations():
    n = 0
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if re.search(r'boxed|\\text|\\mathrm|\$', b): continue
            if len(re.findall(r'<=>|->|<-', b)) != 1: continue
            r = check_equation(b)
            assert r is None, (b, r)
            n += 1
    print('反応式の検査', n, '本 OK')
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
check_structure()
check_all_equations()

from sympy import symbols, solve, Rational
a = symbols('a')
assert solve(1 - a - 2*a, a)[0] == Rational(1, 3)
assert abs(30.0 / 2 - 15.0) < 1e-12
def run(frac, n=1.0):
    so2 = n*frac; h2s = n - so2
    if h2s >= 2*so2:
        S = 3*so2; return ('H2S', h2s - 2*so2, S)
    S = 1.5*h2s; return ('SO2', so2 - h2s/2, S)
r3 = run(0.40); r4 = run(0.25); print(r3, r4)
assert r3[0] == 'SO2' and abs(r3[1] - 0.10) < 1e-9 and sig(r3[2]*32, 2) == 29
assert r4[0] == 'H2S' and abs(r4[1] - 0.25) < 1e-9
for e in ['2H2S + 3O2 -> 2SO2 + 2H2O', '2H2S + SO2 -> 3S + 2H2O', '2H2S + O2 -> 2S + 2H2O', '3H2S + 3/2O2 -> 3S + 3H2O']:
    assert check_equation(e) is None, e

print('ALL OK')
