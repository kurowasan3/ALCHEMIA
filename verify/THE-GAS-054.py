# THE-GAS-054 検証：密閉容器での水の電気分解と水蒸気圧
import sys; sys.path.insert(0, '/home/claude/kurowasan3/alchemia/tools')
from chemlib import check_equation
from sympy import Rational as Q
assert check_equation('2H2O + 2e- -> H2 + 2OH-') is None and check_equation('2H2O -> O2 + 4H+ + 4e-') is None
R = Q('8.31e3'); pv27, pv47 = Q('3.6e3'), Q('1.06e4'); P0 = Q('1.013e5')
air = P0 - pv27; e = Q('0.500')*1930/Q('9.65e4'); ngas = e/2 + e/4
dp = ngas*R*300/Q('0.500'); p27 = P0 + dp
p47 = (air + dp)*Q(320, 300) + pv47; wrong = p27*Q(320, 300)
print(float(air), float(e), float(ngas), float(dp), float(p27), float(p47), float(wrong))
assert e == Q(1, 100) and ngas == Q('0.0075') and round(float(p27), -3) == 139000 and round(float(p47), -3) == 155000 and wrong < p47
print('OK')
