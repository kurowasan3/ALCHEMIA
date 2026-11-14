# 検証スクリプト（グループ g07）  python3 -I この.py
import sys, json, re
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from math import log10, sqrt
from sympy import symbols, nsolve, Rational as Q
ITEM = '/home/claude/chem2/outh/items/THE-ION-050/item.json'
it = json.load(open(ITEM, encoding='utf-8'))

def ce_blocks(s):
    out, i = [], 0
    while True:
        i = s.find('\\ce{', i)
        if i < 0: return out
        j = i + 4; d = 1
        while d:
            if s[j] == '{': d += 1
            elif s[j] == '}': d -= 1
            j += 1
        out.append(s[i+4:j-1]); i = j

eqs = {b for k in ('problem', 'answer', 'solution') for b in ce_blocks(it[k]) if len(re.findall(r'->|<=>', b)) == 1}
for e in sorted(eqs):
    r = check_equation(e); assert r is None, (e, r); print('式OK:', e)

Kw = 1e-14; Kb = 2.0e-5
K3, K2, KM = 1e-38, 1e-15, 1e-11
cFe, cMg, lim = 1e-2, 1e-1, 1e-5
pH = lambda oh: 14 + log10(oh)
OH = {'Fe3 start': (K3/cFe)**(1/3), 'Fe3 full': (K3/lim)**(1/3), 'Fe2 start': (K2/cFe)**0.5, 'Fe2 full': (K2/lim)**0.5, 'Mg start': (KM/cMg)**0.5}
for k, v in OH.items(): print(k, '%.3g' % v, 'pH %.2f' % pH(v))
assert [round(pH(OH[k]), 1) for k in OH] == [2.0, 3.0, 7.5, 9.0, 9.0]
assert not (OH['Fe2 full'] < OH['Mg start'])         # 窓が空
r_lo, r_hi = OH['Fe3 full']/Kb, OH['Mg start']/Kb
print('r range', r_lo, r_hi); assert abs(r_lo - 5e-7) < 1e-12 and abs(r_hi - 0.5) < 1e-12
oh_nh3 = sqrt(Kb*1.0); print('NH3 only OH', oh_nh3, 'IP', cMg*oh_nh3**2); assert cMg*oh_nh3**2 > KM
# (4)
nh3, nh4 = 0.20 - 0.030, 1.0 + 0.030
oh = Kb*nh3/nh4; print('(4) OH', oh, 'pH', pH(oh), 'pH(log3.3=0.52)', 14 - (6 - 0.52))
assert abs(oh - 3.30e-6) < 1e-8 and round(pH(oh), 2) == 8.52
print('Mg IP', cMg*oh**2, 'Fe3 left', K3/oh**3); assert cMg*oh**2 < KM and K3/oh**3 < lim
# (5) Fe2+ を残量 x まで含めて厳密に解く
x = symbols('x')
f = x*(Kb*(Q('0.20') - 2*(Q('0.010') - x))/(1 + 2*(Q('0.010') - x)))**2 - K2
xs = float(nsolve(f, x, 8e-5)); print('(5) Fe2+ left', xs)
oh5 = Kb*0.18/1.02; approx = K2/oh5**2; print('近似', approx)
assert round(xs, 6) == 8.0e-5 and round(approx, 6) == 8.0e-5 and xs > lim
print('すべて OK')
