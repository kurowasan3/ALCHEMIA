# THE-EQL-126 検証：エタンの脱水素平衡（定圧と定積）
from sympy import symbols, solve, Rational, nsimplify
R = 8.31e3; T = 1000; P = 1.00e5
V1 = 108.0
ntot = P*V1/(R*T); x = ntot - 1.00
print('ntot', ntot, 'x', x)
x = round(x, 3)
assert x == 0.300
Kc = x**2/((1 - x)*V1)
print('Kc', Kc)
# 定積 50.0 L で C2H4 0.150 mol
n0 = 0.150 + 0.150**2/(Kc*50.0)
print('n0', n0)
assert round(n0, 2) == 0.53
# 定圧 2.00e5 Pa での反応率
Kp = x**2/(1 - x**2)*P
y = symbols('y', positive=True)
ys = [s for s in solve(y**2/(1 - y**2)*2.00e5 - Kp, y) if 0 < s < 1][0]
print('Kp', Kp, 'conversion at 2e5', ys)
assert round(float(ys), 3) == 0.217
print('conversion exp2', 0.150/n0)
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
assert check_equation('C2H6 <=> C2H4 + H2') is None
print('OK')
