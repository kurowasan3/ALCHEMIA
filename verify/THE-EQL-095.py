"""THE-EQL-095 検証"""
import sys
sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q, N
R, T, Kp = Q(831, 100) * 10**3, 300, Q(150000)
V0 = 1 * R * T / 10**5
P_on = Kp / Q(6, 10)
pN2 = 600000 - Kp
nH2g = Q(4, 10) * Kp / pN2
absorbed = Q(6, 10) - nH2g
V = Q(4, 10) * R * T / pN2
print(N(V0), N(P_on), N(nH2g), N(absorbed), N(V))
assert f'{float(V0):.2g}' == '25' and f'{float(P_on):.2g}' == '2.5e+05'
assert f'{float(absorbed):.2g}' == '0.47' and f'{float(V):.2g}' == '2.2'
assert check_equation('MgH2 <=> Mg + H2') is None  # 形の確認（M を Mg で代用）
print('THE-EQL-095 OK')
