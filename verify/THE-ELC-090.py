# 検証スクリプト THE-ELC-090  python3 -I この.py
ITEM = '/home/claude/chemk/out/items/THE-ELC-090/item.json'
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

from sympy import Rational
F = 96500; Li, S, Ti, C = 6.9, 32, 48, 12
n_host = 56.0e-3 / (Ti + 2*S)
assert abs(n_host - 5.00e-4) < 1e-12
ne = 2.00e-3 * 120 * 60 / F
x = ne / n_host
print('ne', ne, 'x', x, 'ox', 4 - x, 'dm mg', ne*Li*1e3)
assert sig(x, 2) == 0.30 and sig(4 - x, 2) == 3.7 and sig(ne*Li*1e3, 3) == 1.03
t = n_host * F / 2.00e-3 / 60
print('t min', t); assert sig(t, 3) == 402
r = (Ti + 2*S) / (6*C); print('ratio', r); assert sig(r, 2) == 1.6
assert check_equation('TiS2 + Li+ + e- -> LiTiS2') is None and check_equation('Li+ + e- -> Li') is None

print('ALL OK')
