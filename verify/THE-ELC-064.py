# THE-ELC-064 検証：AgCl|Cl- と AgBr|Br- の電池，沈殿の消失
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, symbols, solve
for e in ['AgCl + e- -> Ag + Cl-', 'Ag + Br- -> AgBr + e-']:
    assert check_equation(e) is None
K1, K2 = Q('1.8e-10'), Q('5.4e-13'); c0, b0 = Q('1.0e-3'), Q('0.10')
print(float(K1/c0), float(K2/b0))
x = symbols('x', positive=True)
xs = [s for s in solve(K1*(b0 - x) - K2*(c0 + x), x)][0]
print(float(xs), float(c0 + xs), float(b0 - xs))
assert round(float(c0 + xs), 2) == 0.10 and abs(float(b0 - xs) - 3.0e-4) < 5e-6
# AgCl 0.0500 mol のみ
n = Q('0.0500'); cl, br = c0 + n, b0 - n
agA_sat, agB = K1/cl, K2/br
print(float(agA_sat), float(agB)); assert agA_sat > agB
assert float(agA_sat*1) < 1e-8          # 溶存 Ag+ の量（1 L）は無視できる
print(float(cl), float(br), float(agB), float(agB*cl))
assert round(float(cl), 3) == 0.051 and float(br) == 0.05 and round(float(agB)*1e11, 1) == 1.1 and agB*cl < K1
print(float(n*108)); assert round(float(n*108), 2) == 5.40
print('OK')
