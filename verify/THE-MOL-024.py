# 検証スクリプト THE-MOL-024  python3 -I この.py
ITEM = '/home/claude/chem/out2/items/THE-MOL-024/item.json'
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
def all_equations():
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b:
                eqs.add(b)
    return sorted(eqs)
def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))
def check_structure():
    n_items = len(re.findall(r'\\item\b', it['problem']))
    n_ans = len(re.findall(r'(?:^|\s)\((\d+)\)', it['answer']))
    assert it['problem'].count('\\begin{enumerate}') <= 1
    assert n_items == n_ans, (n_items, n_ans)
    print('小問数', n_items, '答え数', n_ans)
from sympy import Rational as R, symbols, solve
check_structure()
eqs = all_equations(); print(eqs)
for e in eqs: assert check_equation(e) is None, e
assert {'2KClO3 -> 2KCl + 3O2', '4KClO3 -> 3KClO4 + KCl', 'KClO4 -> KCl + 2O2', 'Ag+ + Cl- -> AgCl'} <= set(eqs)
K, Cl, O = 39, R(71,2), 16
M3 = K + Cl + 3*O; M4 = K + Cl + 4*O; MKCl = K + Cl
assert M3 == R(245,2) and M4 == R(277,2)
n0 = R(980,100) / M3; assert n0 == R(8,100)
a, b, c = symbols('a b c')
sol = solve([R(3,2)*a - R(96,100)/32, a + b/4 - R(1,10)*R(32,1000)*10, a + b + c - n0], [a, b, c]); print(sol)
assert sol == {a: R(2,100), b: R(48,1000), c: R(12,1000)}
a, b, c = sol[a], sol[b], sol[c]
mKCl = (a + b/4)*MKCl; m4 = R(3,4)*b*M4; m3 = c*M3; tot = mKCl + m4 + m3
print(float(mKCl), float(m4), float(m3), float(tot)); assert tot == R(980,100) - R(96,100)
assert sig(float(m4/tot*100), 3) == 56.4
lo = (a + R(0)/4) / 10 / R(1,10) * 1000; hi = (a + (n0 - a)/4) / 10 / R(1,10) * 1000
print('滴定範囲', float(lo), float(hi)); assert (lo, hi) == (20, 35)
extra = R(3,2)*c + 2*R(3,4)*b; print('追加 O2', float(extra), float(extra*R(224,10)))
assert extra == R(9,100) and sig(float(extra*R(224,10)), 3) == 2.02
assert R(3,100) + extra == 3*n0/2
assert a/n0 == R(1,4) and (a + b)/n0 == R(85,100)
print('全検証OK:', ITEM)
