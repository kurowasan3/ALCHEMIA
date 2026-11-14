ITEM = '/home/claude/chem2/out/items/INO-OXS-015/item.json'
# 検証スクリプト（グループ g08）  python3 -I この.py
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
    eqs = set()
    for key in ('problem', 'answer', 'solution'):
        for b in ce_blocks(it[key]):
            if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b and not re.search(r'R', b):
                eqs.add(b)
    for e in sorted(eqs):
        r = check_equation(e)
        assert r is None, (e, r)
        print('式OK:', e)
    return eqs

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

def check_parts():
    n_items = it['problem'].count('\\item')
    nums = re.findall(r'(?:^|\s)\((\d+)\)', it['answer'])
    assert nums == [str(i) for i in range(1, n_items + 1)], (n_items, nums)
    print('小問の数OK:', n_items)


check_parts()
eqs = check_all_equations()
assert 'SO3 + H2O -> H2SO4' in eqs and 'H2SO4 + 2NaOH -> Na2SO4 + 2H2O' in eqs
from sympy import Rational as R, symbols, solve, nsimplify
H, O, S = R(1), R(16), R(32)
M_SO3, M_H2SO4 = S + 3*O, 2*H + S + 4*O
assert (M_SO3, M_H2SO4) == (80, 98)
m = 80 + 20 * M_H2SO4 / M_SO3
assert m == R(209, 2); print('(2)', float(m))
w = m / R(98, 100) - 100
print('(3)', float(w)); assert sig(float(w), 3) == 6.63
x = symbols('x')
n = R(216, 10) / 1000 * R(1, 10) / 2 * 10
xs = solve(x / M_SO3 + (1 - x) / M_H2SO4 - n, x)[0]
print('(4)', float(xs)); assert sig(float(xs) * 100, 3) == 26.0

print('全検証OK:', ITEM)
