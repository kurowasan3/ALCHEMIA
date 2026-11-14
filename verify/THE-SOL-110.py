# 検証スクリプト（グループ g03）  python3 -I この.py
ITEM = '/home/claude/chemk/outh/items/THE-SOL-110/item.json'
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
import numpy as np
R = 8.31e3
# データの作り方：M=2.0e5, B=4.0e-6, θ=305 K
M0, B0, th0 = 2.0e5, 4.0e-6, 305.0
c = np.array([2.0, 4.0, 6.0, 8.0])
given = {298.0: [23.8, 45.8, 65.9, 84.2], 323.0: [29.2, 63.3, 102, 146]}
for T, P in given.items():
    A2 = B0*(1 - th0/T)
    true = c*R*T*(1/M0 + A2*c)
    print(T, np.round(true, 2))
    for p, t in zip(P, true): assert sig(t, 3) == p
res = {}
for T, P in given.items():
    y = np.array(P)/c
    s2 = (y[3]-y[0])/6.0; i2 = y[0] - s2*2.0
    sl, ic = np.polyfit(c, y, 1)
    res[T] = dict(slope=s2, icpt=i2, M=R*T/i2, A2=s2/(R*T), Mfit=R*T/ic, A2fit=sl/(R*T))
    print(T, res[T])
r3 = res[323.0]; r2 = res[298.0]
assert round(r3['icpt'], 1) == 13.4 and round(r3['slope'], 3) == 0.608
assert sig(r3['M'], 2) == 2.0e5 and sig(r2['M'], 2) == 2.0e5 and sig(r3['Mfit'], 2) == 2.0e5
assert sig(r3['A2'], 2) == 2.3e-7 and sig(r2['A2'], 2) == -9.3e-8
assert sig(r3['A2fit'], 2) == 2.3e-7 and sig(r2['A2fit'], 2) == -9.3e-8
th = sp.symbols('th')
sol = sp.nsolve((1 - th/298)/(1 - th/323) - r2['A2']/r3['A2'], th, 300)
print('θ', sol); assert sig(float(sol), 3) == 305
solfit = sp.nsolve((1 - th/298)/(1 - th/323) - r2['A2fit']/r3['A2fit'], th, 300)
print('θ(全点の回帰)', solfit); assert sig(float(solfit), 3) == 305
Mapp = 8.00*R*323/146; print('M見かけ', Mapp); assert sig(Mapp, 2) == 1.5e5
print('OK')
