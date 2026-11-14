import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, symbols, solve, N
for e in ['S2O3^2- + 4Br2 + 5H2O -> 2SO4^2- + 8Br- + 10H+', '2S2O3^2- + I2 -> S4O6^2- + 2I-', 'Br2 + 2I- -> I2 + 2Br-',
          'I2 + 2OH- -> IO- + I- + H2O', 'S2O3^2- + 4IO- + 2OH- -> 2SO4^2- + 4I- + H2O', '2S2O3^2- -> S4O6^2- + 2e-']:
    assert check_equation(e) is None, e
c = R(1, 10)
nBr2 = c*R(24, 1000)/2
conc = nBr2/R(10, 1000)
wrong = c*R(8, 1000)/2/R(10, 1000)
e_tot = 2*nBr2
Vmax = e_tot/1/c*1000; Vmin = e_tot/8/c*1000
f = symbols('f')
fv = solve(c*R(8, 1000)*(1 + 7*f) - e_tot, f)[0]
print(conc, wrong, Vmin, Vmax, fv, N(fv, 4))
assert conc == R(12, 100) and wrong == R(4, 100) and Vmin == 3 and Vmax == 24 and fv == R(2, 7)
print('OK')
