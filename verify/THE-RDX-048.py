# THE-RDX-048 検証：H2O2 から O2（MnO2 触媒 と KMnO4 酸化剤）
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
for e in ['2H2O2 -> 2H2O + O2', '2MnO4- + 5H2O2 + 6H+ -> 2Mn^2+ + 5O2 + 8H2O']:
    assert check_equation(e) is None, e
n = 0.500*20.0e-3
O2a = n/2; O2b = n
print('O2 a,b mol', O2a, O2b, 'L', O2a*22.4, O2b*22.4)
assert round(O2a*22.4, 3) == 0.112 and round(O2b*22.4, 3) == 0.224
nMn = n*2/5
print('KMnO4 mL', nMn/0.100*1000)
assert abs(nMn/0.100*1000 - 40.0) < 1e-9
print('OK')
