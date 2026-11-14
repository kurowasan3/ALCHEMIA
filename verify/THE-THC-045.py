"""THE-THC-045 検証"""
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, solve, Symbol
CO2, H2Ol, H2Og, CH4, C3H8, NH3 = -394, -286, -242, -75, -105, -46
dC = CO2
dCH4 = CO2 + 2 * H2Ol - CH4
dC3 = 3 * CO2 + 4 * H2Ol - C3H8
dNH3 = Q(6 * H2Ol - 4 * NH3, 4)
assert (dCH4, dC3, dNH3) == (-891, -2221, -383)
per = {'C': -dC, 'CH4': -dCH4, 'C3H8': Q(-dC3, 3)}
order = sorted(per, key=lambda k: per[k])
assert order == ['C', 'C3H8', 'CH4']
x = Symbol('x'); xs = solve(394 + 383 * x - Q(2221, 3), x)[0]
assert f'{float(xs):.2g}' == '0.9'
dCH4g = CO2 + 2 * H2Og - CH4
assert dCH4g == -803
r1, r2 = 891 / 394, 803 / 394
assert f'{r1:.2g}' == '2.3' and f'{r2:.2g}' == '2'
for e in ['C3H8 + 5O2 -> 3CO2 + 4H2O', '4NH3 + 3O2 -> 2N2 + 6H2O', 'CH4 + 2O2 -> CO2 + 2H2O', 'C + O2 -> CO2']:
    assert check_equation(e) is None
print('THE-THC-045 OK', float(xs), r1, r2)
