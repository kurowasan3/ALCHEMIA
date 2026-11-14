# ORG-ALC-021 検証：二クロム酸カリウムによるアルコールの酸化（電子数の収支）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as R, symbols, solve
eqs = ['Cr2O7^2- + 14H+ + 6e- -> 2Cr^3+ + 7H2O',
       'CH3CH(OH)CH2CH3 -> CH3COCH2CH3 + 2H+ + 2e-',
       '3CH3CH(OH)CH2CH3 + Cr2O7^2- + 8H+ -> 3CH3COCH2CH3 + 2Cr^3+ + 7H2O',
       'CH3CH2CH2CH2OH + H2O -> CH3CH2CH2COOH + 4H+ + 4e-',
       '3CH3CH2CH2CH2OH + 2Cr2O7^2- + 16H+ -> 3CH3CH2CH2COOH + 4Cr^3+ + 11H2O']
for e in eqs:
    assert check_equation(e) is None, (e, check_equation(e))
# (3) 各 0.30 mol：失う電子 第一級→酸 4, 第二級 2, 第三級 0
n = R(30, 100)
cons = {k: n*e/6 for k, e in {'a': 4, 'b': 2, 'c': 4, 'd': 0}.items()}
assert cons == {'a': R(1, 5), 'b': R(1, 10), 'c': R(1, 5), 'd': 0}
# (4)
x, y = symbols('x y')
s = solve([x + y - R(300, 1000), 4*x - R(120, 1000)*6], [x, y])
assert s[x] == R(180, 1000) and s[y] == R(120, 1000)
# (5) 0.0500 mol × 6 = 0.300 mol = 0.150 mol × 2 → 2 電子の酸化
assert R(500, 10000)*6 == R(150, 1000)*2
print('ORG-ALC-021 OK', cons, s)
