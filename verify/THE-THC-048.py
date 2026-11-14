"""THE-THC-048 検証"""
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import symbols, simplify, Rational as Q
for e in ['CH4 + H2O -> CO + 3H2', 'CO + H2O -> CO2 + H2', 'N2 + 3H2 -> 2NH3',
          '3CH4 + 6H2O + 4N2 -> 3CO2 + 8NH3', 'CH4 + 2O2 -> CO2 + 2H2O', '4NH3 + 3O2 -> 2N2 + 6H2O']:
    assert check_equation(e) is None, e
def calc(co2, h2o, ch4, nh3, eta):
    dH = (3 * co2 + 8 * nh3) - (3 * ch4 + 6 * h2o)        # 8 NH3 あたり
    qCH4 = -(co2 + 2 * h2o - ch4)
    qNH3 = -(6 * h2o - 4 * nh3) / 4
    n = Q(3, 8) + dH / 8 / (qCH4 * eta)
    return dH, qCH4, qNH3, n, n / qNH3 / (1 / qCH4)
dH, qCH4, qNH3, n1, r1 = calc(-394, -242, -75, -46, 1)
_, _, _, n2, r2 = calc(-394, -242, -75, -46, Q(1, 2))
assert dH == 127 and qCH4 == 803 and qNH3 == 317
print(float(dH / 8), float(n1), float(n2), float(r1), float(r2))
assert f'{float(dH/8):.3g}' == '15.9' and f'{float(n1):.3g}' == '0.395' and f'{float(n2):.3g}' == '0.415'
assert f'{float(r1):.3g}' == '1' and f'{float(r2):.3g}' == '1.05'
cap = (dH / 8 / (qCH4 * Q(1, 2))) / qNH3 * qCH4
assert f'{float(cap):.2g}' == '0.1', float(cap)
a, b, c, dd = symbols('a b c d')
assert simplify(calc(a, b, c, dd, 1)[4] - 1) == 0      # η=1 では恒等的に 1
print('THE-THC-048 OK', float(cap))
