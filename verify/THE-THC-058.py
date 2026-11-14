import sys, sympy as sp
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
Q1 = 4*394 + 5*286 - 126; Q2 = 4*111 + 5*286 - 126
assert Q1 == 2880 and Q2 == 1748
x = sp.symbols('x'); xs = sp.solve(Q1*x + Q2*(2 - x) - 5190, x)[0]
nCO = 4*(2 - xs); nO2 = sp.Rational(13,2)*xs + sp.Rational(9,2)*(2 - xs)
print('x', float(xs), 'CO', float(nCO), 'O2', float(nO2))
assert abs(float(nCO) - 570/283) < 1e-9
assert round(float(nCO), 1) == 2.0 and round(float(nO2), 1) == 12.0
for e in ['C4H10 + 13/2O2 -> 4CO2 + 5H2O', 'C4H10 + 9/2O2 -> 4CO + 5H2O', 'CO + 1/2O2 -> CO2']:
    assert check_equation(e) is None, e
print('OK')
