# THE-ELC-062 検証：AgCl を含む銀の濃淡電池
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, symbols, solve
for e in ['Ag + Cl- -> AgCl + e-', 'AgCl + e- -> Ag + Cl-', 'Ag+ + e- -> Ag', 'Ag -> Ag+ + e-']:
    assert check_equation(e) is None, e
Ksp = Q('1.8e-10'); A0, B0 = Q('0.10'), Q('0.0010')
print(float(Ksp/A0), float(Ksp/B0))
x = symbols('x'); xs = solve(A0 - x - (B0 + x), x)[0]
print(float(xs), float(A0 - xs), float(xs*Q('9.65e4')))
assert xs == Q('0.0495') and A0 - xs == Q('0.0505') and round(float(xs*96500), -2) == 4800
print('OK')
