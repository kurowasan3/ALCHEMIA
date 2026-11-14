# THE-ACB-044 verification: haloacetic acids, pH without approximation
from sympy import symbols, solve, log, sqrt, Rational, N
x = symbols('x', positive=True)
c = Rational(1, 10)
res = {}
for name, Ka in (('ClCH2COOH', Rational(14, 10000)), ('Cl2CHCOOH', Rational(45, 1000))):
    s = [r for r in solve(x**2 / (c - x) - Ka, x) if r < c][0]
    pH = N(-log(s, 10)); pHa = N(-log(sqrt(Ka * c), 10))
    print(name, N(s), pH, 'approx', pHa, 'alpha', N(s / c))
    res[name] = (round(float(s), 3 if s < 0.02 else 3), round(float(pH), 2), round(float(pHa), 2))
assert res['ClCH2COOH'][1:] == (1.95, 1.93)
assert res['Cl2CHCOOH'][1:] == (1.32, 1.17)
# rounded helper values given in the problem
import math
assert abs(math.sqrt(5.6196) - 2.371) < 1e-3 and abs(math.sqrt(2.0025) - 1.415) < 1e-3
assert abs(math.log10(1.12) - 0.049) < 1e-3 and abs(math.log10(4.83) - 0.684) < 1e-3
assert abs(math.log10(1.18) - 0.072) < 1e-3 and abs(math.log10(6.71) - 0.827) < 1e-3
# ordering by electronegativity / literature pKa
pKa = {'F': 2.59, 'Cl': 2.87, 'Br': 2.90, 'I': 3.18}
assert sorted(pKa, key=pKa.get) == ['F', 'Cl', 'Br', 'I']
assert 2.83 < 4.08 < 4.87
print('OK')
