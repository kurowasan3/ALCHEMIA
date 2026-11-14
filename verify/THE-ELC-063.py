# THE-ELC-063 検証：電解後の Cu2+ のヨウ素滴定と電流効率
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
for e in ['2Cu^2+ + 4I- -> 2CuI + I2', 'I2 + 2S2O3^2- -> 2I- + S4O6^2-', 'Cu^2+ + 2e- -> Cu', '2H2O -> O2 + 4H+ + 4e-', '2H+ + 2e- -> H2']:
    assert check_equation(e) is None
s = Q('0.100')*Q('22.75'); c = s/Q('5.00'); remain = c*2; dep = 1 - remain
e = Q('5.00')*3860/Q('9.65e4'); theo = e/2; eff = dep/theo*100
h2 = (e - 2*dep)/2; o2 = e/4
print(float(c), float(dep), float(e), float(eff), float(h2), float(o2))
assert c == Q('0.455') and dep == Q('0.09') and e == Q('0.2') and eff == 90 and h2 == Q('0.01') and o2 == Q('0.05')
print('OK')
