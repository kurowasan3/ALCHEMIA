import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, N
for e in ['2Cu^2+ + 4I- -> 2CuI + I2', 'I2 + 2S2O3^2- -> 2I- + S4O6^2-', '2Fe^3+ + 2I- -> 2Fe^2+ + I2']:
    assert check_equation(e) is None, e
nS = R(1, 10)*R(196, 10)/1000
nCu = nS                      # Cu2+ 2 -> I2 1 -> S2O3 2
M = R(635, 10) + 32 + 64 + 5*R(180, 10)
assert M == R(4990, 20)
pur = nCu*M/R(1, 2)*100
kI = 2*nCu*(39 + 127)
print('M', N(M), 'purity', N(pur, 6), 'KI g', N(kI, 5))
assert round(float(pur), 1) == 97.8 and round(float(kI), 2) == 0.65
print('OK')
