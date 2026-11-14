# 検証スクリプト（グループ g01）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
ITEM = '/home/claude/chem/out/items/THE-MOL-014/item.json'
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
        if len(re.findall(r'->|<=>|<-', b)) == 1 and '\\boxed' not in b and not re.search(r'_n|\bn[A-Z]|R', b):
            eqs.add(b)
for e in sorted(eqs):
    r = check_equation(e)
    assert r is None, (e, r)
    print('式OK:', e)

def near(a, b, rel=0.005):
    return abs(a - b) <= rel * abs(b)

def sig(x, n):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))) + (n - 1))

import sympy as sp
NA, k = sp.symbols('N_A k', positive=True)   # k=2（1 mol の粒子数が k 倍）
NA2 = k * NA; M = sp.Integer(44); M2 = k * M; n = sp.Integer(44) / M; n2 = sp.Integer(44) / M2
Vm = sp.Rational(224, 10); Vm2 = Vm * NA2 / NA
rho = 32 / Vm; rho2 = (32 * k) / Vm2
N = n * NA; N2 = n2 * NA2
ratios = [NA2 / NA, sp.Integer(1), M2 / M, n2 / n, Vm2 / Vm, sp.simplify(rho2 / rho), sp.simplify(N2 / N)]
ratios = [sp.simplify(r.subs(k, 2)) for r in ratios]
print(ratios); assert ratios == [2, 1, 2, sp.Rational(1, 2), 2, 1, 1]

print('全検証OK:', ITEM)
