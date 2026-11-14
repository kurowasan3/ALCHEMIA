ITEM = '/home/claude/chemk/out/items/THE-GAS-070/item.json'
# 検証スクリプト（グループ g04）  python3 -I この.py
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

eqs = set()
for key in ('problem', 'answer', 'solution'):
    for b in ce_blocks(it[key]):
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and '\\text' not in b and '\\mathrm' not in b:
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import sympy as sp
M = sp.symbols('M')
gS = 9.80/2.00e-4
P = lambda m, MM: 1.00e5 - gS*(MM+m)
Ms = sp.solve(sp.Eq(P(0.300, M)*20.0, P(0.800, M)*29.6), M)[0]; print('M =', Ms)
assert sig(float(Ms), 2) == 0.20
for MM in (0.20, float(Ms)):
    PV = P(0.300, MM)*20.0e-6
    n = PV/(8.31*300); mmax = (1.00e5 - PV/50.0e-6)/gS - MM; V4 = PV/(1.00e5 + gS*(MM+0.300))*1e6
    print(MM, n, mmax, V4)
    assert sig(n, 3) == 6.06e-4 and sig(mmax, 2) == 1.2 and sig(V4, 2) == 12

print('全検証OK:', ITEM)
