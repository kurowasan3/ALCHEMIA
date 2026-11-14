# THE-MOL-041 検証：密閉容器での固体単体の燃焼と圧力
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, N, symbols, solve
eqs = ['S + O2 -> SO2', 'C + O2 -> CO2', '2C + O2 -> 2CO', '4P + 5O2 -> P4O10',
       '2Mg + O2 -> 2MgO', '4Al + 3O2 -> 2Al2O3', '3C + 2O2 -> 2CO + CO2']
for e in eqs:
    r = check_equation(e); print(e, '->', r or 'OK'); assert r is None
# (1) 0.48 g の単体 -> 酸化物 0.96 g，圧力不変
O2 = (Q('0.96') - Q('0.48')) / 32
M = Q('0.96') / O2
print('O2 mol', O2, 'M =', M); assert M == 64
# 酸化物 1 分子に O2 1 分子分（O 原子 2 個）→ 残り 64-32=32 が単体の原子分
assert M - 32 == 32
# (2) C 0.36 g = 0.030 mol, O2 0.020 mol，すべて消費
x, y = symbols('x y')
s = solve([x + y - Q('0.030'), x + 2*y - 2*Q('0.020')], [x, y])
print('CO', s[x], 'CO2', s[y], 'P ratio', (s[x] + s[y]) / Q('0.020'))
assert s[x] == Q('0.020') and s[y] == Q('0.010') and (s[x]+s[y])/Q('0.020') == Q(3, 2)
# 質量保存
mass_before = Q('0.36') + Q('0.020')*32
mass_after = s[x]*28 + s[y]*44
assert mass_before == mass_after; print('mass', mass_before, mass_after)
# (3) 気体の物質量変化（消費 O2 1 mol あたり生成気体 mol）
print('C:1->1 不変, S:1->1 不変, P: 5->0 減少, Mg: 1->0 減少, Al: 3->0 減少')
print('OK')
